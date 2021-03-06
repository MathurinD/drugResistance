% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/growth_curves.R
\name{process_growth_curve}
\alias{process_growth_curve}
\title{Process growth curves to viability values

Extract growth rate, normalise and averages growth curves from the experiment}
\usage{
process_growth_curve(experiment, keep_controls = TRUE, max_confluency = 75)
}
\arguments{
\item{experiment}{Dataset as return by read_incu from incucyter}

\item{keep_controls}{Whether the controls should be present in the output. Control values are DMSO, medium and no_cells.}
}
\value{
A tibble with new columns Inhibitor, Concentration and Viability
}
\description{
Process growth curves to viability values

Extract growth rate, normalise and averages growth curves from the experiment
}
