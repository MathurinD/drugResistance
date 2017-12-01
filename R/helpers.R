#' Convert readable concentration format to numerals
#'
#' Convert from the readable "XX[pnµm]M" format to a numeral M format
#' @param concentrations A vector of readable concentrations
#' @return A numeric vector corresponding to the converted values
#' @export
# NOTE: does not manage non valid formats correctly
concentration_values = function(concentrations) {
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
#' @value A vector of numerals
#' @export
sigmoid <- function(params, xx) {
  return( 1 / (1 + exp(-params[1] * (xx - params[2]))) )
}

#' Fit a sigmoid
#'
#' Helper function to fit a sigmoid with a=1 to data using the LM algorithm
#' @param fit_data A tibble or data.frame with columns 'Concentration_value' and 'Viability' corresponding defining the points to fit
#' @return A vector of size 3 corresponding to parameters b and c, and the residual of the fit
#' @export
#' @seealso sigmoid
fit_sigmoid <- function(fit_data) {
    fits = list()
    try({ fits[[length(fits)+1]] = nlsLM( Viability ~ 1/(1+exp(-(log10(Concentration_value) - cc)*bb)), start=list(bb=-0.01, cc=-6), data=fit_data) }, silent=TRUE)
    try({ fits[[length(fits)+1]] = nlsLM( Viability ~ 1/(1+exp(-(log10(Concentration_value) - cc)*bb)), start=list(bb=-0.03, cc=-6), data=fit_data) }, silent=TRUE)
    try({ fits[[length(fits)+1]] = nlsLM( Viability ~ 1/(1+exp(-(log10(Concentration_value) - cc)*bb)), start=list(bb=-0.1, cc=-6), data=fit_data) }, silent=TRUE)
    if (length(fits) > 0) {
        fits = sapply(fits, function(ff) { c(coef(ff), residual=sum(resid(ff)^2)) })
        return(fits[c("bb", "cc"), order(fits["residual",])[1]])
    } else {
        return(c(-0.01, 0, NaN))
    }
}
