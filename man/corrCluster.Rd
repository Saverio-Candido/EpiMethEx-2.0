\name{corrCluster}
\alias{corrCluster}
\title{Clustering of CG probesets according to correlation median levels}
\description{
The \code{corrCluster} function performs the clustering of consecutive CG probesets according to correlation values between methDNA levels of each CG probeset and expression of relative genes.
CG probesets are clusterized and annotated according to the following groups: strongPosCorr: (r > 0.7), moderatePosCorr (0.3 ≤ r ≤ 0.7), weakPosCorr (0 ≤ r < 0.3), strongNegCorr (r < -0.7), moderateNegCorr (-0.7 ≤ r ≤ -0.3), and weakNegCorr (-0.3 < r ≤ 0)

}
\usage{
corrCluster(cgAnnotation, CGclusterAnn, corrMatrix, by = "RefGene_Name", nCGclusters = 1)
}
\arguments{
  \item{cgAnnotation}{
    A data frame containing the CG probeset annotations, including  IDs, genomic coordinates, and associated gene Names, RefGene  Accessions and  RefGene  Regions.
  }
  \item{CGclusterAnn}{
    A data frame of RefGene_Name or RefGene_Accession CG cluster Annotation.
  }
  \item{corrMatrix}{
    A data frame containing the correlation r and significance (p value) values computed between each CG probeset and relative genes (filter corrMatrix according to p values).
  }
  \item{by}{
    A character string specifying the method of CG clustering. Options are \code{"RefGene_Name"} (default), or \code{"RefGene_Accession"}.
  }
  \item{nCGclusters}{
    The minimum number of CG probesets for each CG cluster.
  }
  }
\value{
The resulting data frame is composed of genomicPosition or RefGene_Accession CG cluster Annotation details with additional columns including:
\itemize{
  \item \code{cor}: correlation median values for each CG probeset.
  \item \code{p}: p value relative to correlation analysis for each CG probeset.
  \item \code{corrClusterName}: A unique identifier for each correlation CG cluster.
  \item \code{ClusterCoverage}: The ratio between the number of probesets of each correlation CG cluster and relative CG cluster Annotation.
  \item \code{cgN}: The number of CG probesets in the correlation CG cluster.
  \item \code{CGsize}: The size (bp) of the correlation CG cluster.
  \item \code{CGstartcluster}: The starting Chr position of correlation CG cluster.
  \item \code{CGstopcluster}: The ending Chr position of correlation CG cluster.
  \item \code{Corr.median.cluster}: cumulative median of correlation levels relative to CG probesets belonging to correlation CG cluster.
  \item \code{Corr.IQR.cluster}: Interquartile of the median of correlation levels relative to CG probesets belonging to correlation CG cluster.
  \item \code{refGeneIso}: (available only for RefGene_Name clustering option) List of RefGene Accession and region related to correlation CG cluster.
}
}
\examples{
# Example usage of corrCluster function
# To testing the function run the code using cgAnnotation, CGclusterAnn (select CGclusterAnn_Name for "RefGene_Name" clustering method and CGclusterAnn_Accession for "RefGene_Accession" clustering method and rename as CGclusterAnn), corrMatrix test datasets available in EpiMethEx2.0 package. The appropiate method of CG clustering and the minimum number of CG probesets for each CG cluster have be set:

library(EpiMethEx2.0)
data(cgAnnotation)
data(CGclusterAnn_Name)
CGclusterAnn <- CGclusterAnn_Name # rename CGclustAnn test dataset
data(corrMatrix)
corrCluster <- corrCluster(cgAnnotation, CGclusterAnn, corrMatrix, by = "RefGene_Name", nCGclusters = 1)


head(corrCluster)
}
