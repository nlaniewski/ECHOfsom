% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions_fsom_plots.R
\name{pheat}
\alias{pheat}
\title{Pheatmap of FlowSOM Clusters/Nodes}
\usage{
pheat(
  dat.mat,
  color.type = c("divergent", "sequential"),
  break.vec = seq(0, 1, by = 0.05),
  ...
)
}
\arguments{
\item{dat.mat}{numeric matrix}

\item{color.type}{c("divergent","sequential"); heatmap colorscale}

\item{break.vec}{a numeric vector for generating color breaks}

\item{...}{additional arguments for pheatmap::pheatmap}
}
\value{
a pheatmap
}
\description{
Pheatmap of FlowSOM Clusters/Nodes
}
