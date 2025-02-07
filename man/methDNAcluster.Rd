\name{methDNAcluster}
\alias{methDNAcluster}
\title{Clustering of CG probesets according to methDNA median levels}
\description{
The \code{methDNAcluster} function performs the clustering of consecutive CG probesets according to methDNA median levels ( hypo- = median Beta values < 0.2; partially = 0.2 ≤ median Beta values ≤ 0.6; hyper-methylated = median Beta values > 0.6) and by genomic Position or RefGene Accession CG cluster Annotation.
}
\usage{
methDNAcluster(cgAnnotation, CGclusterAnn, methDNAmatrix, by = "genomePosition", nCGclusters = 1)
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
  \item \code{cgclusterName}: A unique identifier for each methDNA CG cluster.
  \item \code{ClusterCoverage}: The ratio between the number of probesets of each methDNA CG cluster and relative CG cluster Annotation.
  \item \code{cgN}: The number of CG probesets in the methDNA CG cluster.
  \item \code{CGsize}: The size (bp) of the methDNA CG cluster.
  \item \code{CGstartcluster}: The starting Chr position of methDNA CG cluster.
  \item \code{CGstopcluster}: The ending Chr position of methDNA CG cluster.
  \item \code{CG.median.cluster}: cumulative median of methDNA levels relative to CG probesets belonging to methDNA CG cluster.
  \item \code{CG.IQR.cluster}: Interquartile of the median of methDNA levels relative to CG probesets belonging to methDNA CG cluster.
  \item \code{gene}: (available only for genomePosition clustering option) List of genes (refGene Name) related to methDNA CG cluster. The "." character indicates that no genes are associated to the cluster.
}
}
\examples{
# Example usage of methDNAcluster function
# To testing the function run the code using cgAnnotation, CGclusterAnn (select CGclusterAnn_Pos for "genomePosition" clustering method and CGclusterAnn_Accession for "RefGene_Accession" clustering method and rename as CGclusterAnn), methDNAmatrix test datasets available in EpiMethEx2.0 package, and select appropiate method of CG clustering as follows:
library(EpiMethEx2.0)
data(cgAnnotation)
data(CGclusterAnn_Pos)
CGclusterAnn <- CGclusterAnn_Pos # rename CGclustAnn test dataset
data(methDNAmatrix)

methDNAcluster <- methDNAcluster(cgAnnotation, CGclusterAnn, methDNAmatrix, by = "genomePosition", nCGclusters = 1)

head(methDNAcluster)
}
