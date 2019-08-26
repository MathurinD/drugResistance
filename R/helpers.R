############################## helpers.R ###################################################

#' Convert readable concentration format to numerals
#'
#' Convert from the readable "XX[pnum]M" format to a numeral M format
#' @param concentrations A vector of readable concentrations
#' @return A numeric vector corresponding to the converted values
#' @export
concentration_values = function(concentrations) {
    # NOTE: does not manage non valid formats correctly
    scaling = list("nM"=1e-9, "uM"=1e-6, "mM"=1e-3, "k"=1e3, "M"=1e6)
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
#' @return A vector of size 3 corresponding to parameters b and c, and the residual of the fit
#' @seealso sigmoid
#' @export
fit_sigmoid <- function(fit_data) {
    fits = list()
    fit_data = fit_data %>% mutate(Viability = sapply(Viability, max, -0.1, na.rm=TRUE)) # Threshold negative viability for the fit as it seems to ruin it
    # IC50s tend to be between 1 and 100 uM
    try({ fits[[length(fits)+1]] = nlsLM( Viability ~ 1/(1+exp(-(log10(Concentration_value) - ic50)*slope)), start=list(slope=-0.01, ic50=-6), data=fit_data) }, silent=TRUE)
    try({ fits[[length(fits)+1]] = nlsLM( Viability ~ 1/(1+exp(-(log10(Concentration_value) - ic50)*slope)), start=list(slope=-0.03, ic50=-6), data=fit_data) }, silent=TRUE)
    try({ fits[[length(fits)+1]] = nlsLM( Viability ~ 1/(1+exp(-(log10(Concentration_value) - ic50)*slope)), start=list(slope=-0.1, ic50=-6), data=fit_data) }, silent=TRUE)
    try({ fits[[length(fits)+1]] = nlsLM( Viability ~ 1/(1+exp(-(log10(Concentration_value) - ic50)*slope)), start=list(slope=-0.1, ic50=-4), data=fit_data) }, silent=TRUE)
    try({ fits[[length(fits)+1]] = nlsLM( Viability ~ 1/(1+exp(-(log10(Concentration_value) - ic50)*slope)), start=list(slope=-0.03, ic50=-4), data=fit_data) }, silent=TRUE)
    try({ fits[[length(fits)+1]] = nlsLM( Viability ~ 1/(1+exp(-(log10(Concentration_value) - ic50)*slope)), start=list(slope=-0.01, ic50=-4), data=fit_data) }, silent=TRUE)
    if (length(fits) > 0) {
        fits = lapply(fits, function(ff) { list(coef=coef(ff), residual=sum(resid(ff)^2), model=ff) })
        bfid = order(sapply(fits, function(ff) { ff$residual }))[1]
        return( fits[[bfid]] )
    } else {
        return(list(coef=c(-0.01, 0), residual=NaN, conf=matrix(c(-10, -10, 10, 10), nrow=2), model=NULL))
    }
}
