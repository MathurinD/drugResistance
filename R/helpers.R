############################## helpers.R ###################################################

#' Convert readable concentration format to numerals
#'
#' Convert from the readable "XX[pnum]M" format to a numeral M format
#' @param concentrations A vector of readable concentrations
#' @return A numeric vector corresponding to the converted values
#' @export
concentration_values = function(concentrations) {
    # NOTE: does not manage non valid formats correctly
    scaling = list("nM"=1e-9, "µM"=1e-6, "uM"=1e-6, "mM"=1e-3, "k"=1e3, "M"=1e6)
    numbers = gsub("[^0-9.]*", "", concentrations)
    modifiers = gsub("[0-9.]*", "", concentrations)
    return( as.numeric(numbers) * unlist(sapply(modifiers, function(xx){ if (xx %in% names(scaling)){scaling[[xx]]} else{1} })) )
}

#' The sigmoid function
#'
#' The sigmoid function: a / (1 + exp(-b*(x-c)))
#' @param params A vector of size 2 for the 'b' and 'c' parameters
#' @param xx A vector of numerals to which the function should be applied
#' @return A vector of numerals
#' @export
sigmoid <- function(params, xx) {
  return( 1 / (1 + exp(-params[1] * (xx - params[2]))) )
}

#' Fit a sigmoid
#'
#' Helper function to fit a sigmoid with a=1 to data using the LM algorithm
#' @param fit_data A tibble or data.frame with columns 'Concentration_value' and 'Viability' corresponding defining the points to fit
#' @param slopes Initial values to start the fit of slopes
#' @param ic50s Initial values to start the fit of ic50s
#' @return A vector of size 3 corresponding to parameters b and c, and the residual of the fit
#' @seealso sigmoid
#' @export
fit_sigmoid <- function(fit_data, slopes=c(-0.01, -0.03, -0.1), ic50s=c(-4, -5, -6)) {
    fit_data = fit_data %>% mutate(Viability = sapply(Viability, max, -0.1, na.rm=TRUE)) # Threshold negative viability for the fit as it seems to ruin it
    fit_data$Concentration_value %>% range(na.rm=TRUE) %>% log10 -> crange
    # Try fitting with each combination of slopes and IC50s as starting points
    init_fits = lapply(slopes, function(xx) {
        lapply(ic50s, function(yy) {
            # IC50s tend to be between 1 and 100 uM
            try({ nlsLM( Viability ~ 1/(1+exp(-(log10(Concentration_value) - ic50)*slope)), start=list(slope=xx, ic50=yy), data=fit_data, upper=c(slope=-0.00001, ic50=10), lower=c(slope=-100, ic50=-15), control=nls.lm.control(epsfcn=1e-5)) }, silent=TRUE)
        })
    }) %>% flatten
    # Remove fits that failed
    fits = list()
    for (ff in init_fits) {
        if (length(ff) > 1) {
            fits[[length(fits)+1]] = ff
        }
    }
    if (length(fits) > 0) {
        fits = lapply(fits, function(ff) { list(coef=coef(ff), residual=sum(resid(ff)^2), model=ff) })
        bfid = order(sapply(fits, function(ff) { ff$residual }))[1]
        return( fits[[bfid]] )
    } else {
        return(list(coef=c(-0.01, 0), residual=NaN, conf=matrix(c(-10, -10, 10, 10), nrow=2), model=NULL))
    }
}

#' Find the concentration for a given inhibition
#'
#' @param nls_model Model fitted by fit_drug_sensitivity
#' @param target_effect Percent viability that should be reached
#' @param explore Starting range to look for the concentration, in log10
#' @param precision How close to the target effect can the algorithm stop
#' @value log10(target_concentration)
#' @description Find the target concentration by predicting the effect on a range of concentration and shrinking from there
#' @rdname findIC
find_effective_concentration <- function(nls_model, target_effect=0.5, explore=c(-9:1), precision=0.01) {
    best=10
    explore = sort(explore)
    while(best-target_effect > precision) {
        prediction = predict(plots$fits$Venetoclax$model, list(Concentration_value=10^explore))
        if (prediction[1] < target_effect) {stop('First value of the explore range is already more efficient than the target effect')}
        if (prediction[length(prediction)] > target_effect) {stop('Last value of the explore range does not reach the target effect')}
        ii = 2
        while (prediction[ii] > target_effect) {ii = ii+1}
        best = prediction[ii-1]
        explore = seq(explore[ii-1], explore[ii], length.out=10)
    }
    return(explore[1])
}
#' @rdname findIC
findIC20 <- function(nls_model, explore=c(-9:1), precision=0.01) { find_effective_concentration(nls_model, target_effect=0.8, explore, precision)}
