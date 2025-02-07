\name{CGclustAnn}
\alias{CGclustAnn}
\title{Clustering of CG probesets based on Genomic Position, RefGene Name, or RefGene Annotation}
\description{
The \code{CGclustAnn} function performs the clustering of consecutive CG probesets according to their genomic coordinates. The analysis may be performed considering the CG probeset position within gene or gene regions .
}
\usage{
CGclustAnn(cgAnnotation, by = "genomePosition", nCGclusters = 3, CGclusterSize = 200)
}
\arguments{
  \item{cgAnnotation}{
    A data frame containing the CG probeset annotations, including  IDs, genomic coordinates, and associated gene Names, RefGene  Accessions and  RefGene  Regions.
  }
  \item{by}{
    A character string specifying the method  of CG clustering. Options are \code{"genomePosition"} (default), \code{"RefGene_Name"}, or \code{"RefGene_Accession"}.
  }
  \item{nCGclusters}{
    The minimum number of CG probesets for each CG cluster.
  }
  \item{CGclusterSize}{
    The maximum base pair distance between adjacent CG probesets.
  }
}
\value{
The resulting data frame is composed of the CG probesets annotation details with additional columns including:
\itemize{
  \item \code{clusterName}: A unique identifier for each CG cluster.
  \item \code{nCG}: The number of CG probesets in the CG cluster.
  \item \code{size}: The size (bp) of the CG cluster.
  \item \code{chrcluster}: The chromosome location of CG cluster.
  \item \code{startcluster}: The starting Chr position of CG cluster.
  \item \code{stopcluster}: The ending Chr position of CG cluster.
  \item \code{clusterType}: The type of executed CG clustering.
}
}
\examples{
# Example usage of CGclustAnn function
# To testing the function run the code using cgAnnotation test dataset available in EpiMethEx2.0  package as follows:
library(EpiMethEx2.0)
data(cgAnnotation)
CGclustAnn <- CGclustAnn(cgAnnotation, by = "genomePosition", nCGclusters = 2, CGclusterSize = 200)

head(CGclustAnn)
}
