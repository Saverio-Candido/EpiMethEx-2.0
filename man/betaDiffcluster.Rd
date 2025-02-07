\name{betaDiffcluster}
\alias{betaDiffcluster}
\title{Clustering of CG probesets according to Beta difference median levels}
\description{
The \code{betaDiffcluster} function performs the clustering of consecutive CG probesets according to Beta difference values computed between the comparison groups.
CG probesets are clusterized and annotated according to the following groups: Strong methylation  (Beta diff > 0.5), Weak methylation (≥ 0 Beta diff  ≤ 0.5), Strong de-methylation (Beta diff < -0.5), Weak de-methylation ( ≥ -0.5 Beta diff ≤ 0).
}
\usage{
betaDiffcluster(cgAnnotation, CGclusterAnn, methDNAmatrix, betaDiffmatrix, by = "genomePosition", nCGclusters = 1)
}
\arguments{
  \item{cgAnnotation}{
    A data frame containing the CG probeset annotations, including  IDs, genomic coordinates, and associated gene Names, RefGene  Accessions and  RefGene  Regions.
  }
  \item{CGclusterAnn}{
    A data frame of genomicPosition or RefGene_Accession CG cluster Annotation.
  }
  \item{methDNAmatrix}{
    A data frame containing the median methDNA values of CG probesets computed for all samples.
  }
   \item{betaDiffmatrix}{
    A data frame containing the Beta difference and significance (p value) values of CG probesets computed between the comparison groups (e.g. tumor vs normal) (filter betaDiffmatrix according to p values).
  }
  \item{by}{
    A character string specifying the method of CG clustering. Options are \code{"genomePosition"} (default), or \code{"RefGene_Accession"}.
  }
  \item{nCGclusters}{
    The minimum number of CG probesets for each CG cluster.
  }
  }
\value{
The resulting data frame is composed of genomicPosition or RefGene_Accession CG cluster Annotation details with additional columns including:
\itemize{
  \item \code{medianvalue}: methDNA median value of each CG probeset.
  \item \code{beta_diff}: beta difference value of each CG probeset computed between the comparison groups.
  \item \code{pValue}: significance of beta difference analysis for each CG probeset.
  \item \code{cgclusterName}: A unique identifier for each betaDiff CG cluster.
  \item \code{ClusterCoverage}: The ratio between the number of probesets of each betaDiff CG cluster and relative CG cluster Annotation.
  \item \code{cgN}: The number of CG probesets in the betaDiff CG cluster.
  \item \code{CGsize}: The size (bp) of the betaDiff CG cluster.
  \item \code{CGstartcluster}: The starting Chr position of betaDiff CG cluster.
  \item \code{CGstopcluster}: The ending Chr position of betaDiff CG cluster.
  \item \code{BetaDiff.median.cluster}: cumulative median of Beta difference levels relative to CG probesets belonging to betaDiff CG cluster.
  \item \code{CG.IQR.cluster}: Interquartile of the Beta difference levels relative to CG probesets belonging to betaDiff CG cluster.
  \item \code{gene}: (available only for genomePosition clustering option) List of genes (refGene Name) related to betaDiff CG cluster. The "." character indicates that no genes are associated to the cluster.
}
}
\examples{
# Example usage of betaDiffcluster function
# To testing the function run the code using cgAnnotation, CGclusterAnn (select CGclusterAnn_Pos for "genomePosition" clustering method and CGclusterAnn_Accession for "RefGene_Accession" clustering method and rename as CGclusterAnn), methDNAmatrix, and betaDiffmatrix test datasets available in EpiMethEx2.0 package. the appropiate method of CG clustering and the minimum number of CG probesets for each CG clusterhave be set:
library(EpiMethEx2.0)
data(cgAnnotation)
data(CGclusterAnn_Pos)
CGclusterAnn <- CGclusterAnn_Pos # rename CGclustAnn test dataset
data(methDNAmatrix)
data(betaDiffmatrix)
betaDiffcluster <- betaDiffcluster(cgAnnotation, CGclusterAnn, methDNAmatrix, betaDiffmatrix, by = "genomePosition", nCGclusters = 1)

head(betaDiffcluster)
}
