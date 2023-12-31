% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/do_km_in_scenarios.R
\name{do_km_in_scenarios}
\alias{do_km_in_scenarios}
\title{Do k-means clustering in scenarios}
\usage{
do_km_in_scenarios(scenario_names, scenarios, n_repeats, nrows)
}
\arguments{
\item{scenario_names}{Names of scenarios (Output of the make_scenarios function)}

\item{scenarios}{List where each element is a parameter list (Output of the make_scenarios function)}

\item{n_repeats}{Number of repeats}

\item{nrows}{Number of rows (for plotting)}
}
\value{
A plot with subplots of the silhouette scores resulting from k-means clustering in each scenario.
}
\description{
Do k-means clustering in scenarios
}
\examples{
do_km_in_scenarios(scenario_names = scenario_names,
scenarios = scenario_params_list,
n_repeats = 100,
nrows = 3)

}
