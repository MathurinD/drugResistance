% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/synergy_calculations.R
\name{get_synergy_table}
\alias{get_synergy_table}
\title{Get a valid input for ReshapeData

Get a synergy table that can be used by the function ReshapeData synergyfinder.}
\usage{
get_synergy_table(pdata, restrict = FALSE)
}
\arguments{
\item{pdata}{A tibble as outputed from process_growth_curves. Must have columns 'Treatment' with format drugRow_rowConcentration+drugCol_colConcentation and 'Viability'}

\item{scale}{Whether the data should be strictly restricted between 0 and 100. Apply a cutoff, not a scaling.}
}
\value{
A tibble with column names as required by ReshapeData
}
\description{
Get a valid input for ReshapeData

Get a synergy table that can be used by the function ReshapeData synergyfinder.
}