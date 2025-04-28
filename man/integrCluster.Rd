\name{integrCluster}
\alias{integrCluster}
\title{Integration of methDNA, Corr, and betaDiff CG clusters}
\description{
The \code{integrCluster} function performs the merging of CG clusters obtained from methDNA, correlation, and Beta difference clustering analysys according to RefGene_Accession.
The function allows to select consecutive CG probesets showing concordance for methDNA, correlation, and Beta difference values, generating the new integrated CG clusters.
}
\usage{
integrCluster(methDNAmatrix, betaDiffmatrix, corrMatrix, diffExprMatrix, CGclusterAnn, methDNAclusters, betaDiffclusters, corrClusters)
}
\arguments{
   \item{methDNAmatrix}{
    A data frame containing the median methDNA values of CG probesets computed for all samples.
  }
   \item{betaDiffmatrix}{
    A data frame containing the Beta difference and significance (p value) values of CG probesets computed between the comparison groups (e.g. tumor vs normal) (filter betaDiffmatrix according to p values).
  }
   \item{corrMatrix}{
    A data frame containing the correlation r and significance (p value) values computed between each CG probeset and relative genes (filter corrMatrix according to p values).
  }
  \item{diffExprMatrix}{
    A data frame containing differential analysis of gene expression between each tumor type and normalPool. The mean of expression levels of comparison groups, Fold Change, and p value are reported.
  }
   \item{CGclusterAnn}{
    A data frame of RefGene_Accession CG cluster Annotation.
  }
  \item{methDNAclusters}{
    A data frame containing methDNA CG clusters obtained from clustering procedure according to RefGene_Accession CG cluster Annotation (nCGcluster ≥ 2, CG probeset distance ≤ 200 bp)
  }
  \item{betaDiffclusters}{
    A data frame containing betaDiff CG clusters obtained from clustering procedure according to RefGene_Accession CG cluster Annotation (nCGcluster ≥ 2, CG probeset distance ≤ 200 bp)
  }
  \item{corrClusters}{
    A data frame containing Corr CG clusters obtained from clustering procedure according to RefGene_Accession CG cluster Annotation (nCGcluster ≥ 2, CG probeset distance ≤ 200 bp)
  }
  }
\value{
The resulting data frame is composed of RefGene_Accession CG cluster Annotation details with additional columns including:
\itemize{
  \item \code{integratedCluster}: A unique identifier for each integrated CG cluster.
  \item \code{ClusterGenCoordinate}: Genomic coordinates of integrated CG cluster.
  \item \code{CG.median.methDNA}: median methDNA levels of each CG probeset belonging to integrated CG cluster.
  \item \code{CG.median.corr}: median corr levels of each CG probeset belonging to integrated CG cluster.
  \item \code{CG.coor.pValue}: p value relative to correlation analysis of each CG probeset belonging to integrated CG cluster.
  \item \code{CG.median.betaDiff}: median betaDiff levels of each CG probeset belonging to integrated CG cluster.
  \item \code{CG.betaDiff.pValue}: p value relative to Beta difference analysis of each CG probeset belonging to integrated CG cluster.
  \item \code{cluster.median.methDNA}: median methDNA levels of each integrated CG cluster.
  \item \code{cluster.median.corr}: median corr levels of each integrated CG cluster.
  \item \code{cluster.median.betaDiff}: median betaDiff levels of each integrated CG cluster.
  \item \code{mean.expr.Ref}: mean expression levels of gene relative to integrated CG cluster computed for reference group.
  \item \code{logFC_Ref.vs.CTRL}: LogFC of gene expression differential analysis between reference group and control group.
  \item \code{pValue_LogFC.expr}: p value of gene expression differential analysis between reference group and control group.
}
}
\examples{
# Example usage of integrCluster function
# To testing the function run the code selecting CGclusterAnn_Accession (rename it as CGclusterAnn), methDNAmatrix, betaDiffmatrix, corrMatrix, diffExprMatrix, methDNAclusters, betaDiffclusters, and corrClusters test datasets available in EpiMethEx2.0 package:

library(EpiMethEx2.0)

data("methDNAmatrix")
data("corrMatrix")
data("betaDiffmatrix")
data("diffExprMatrix")

data("CGclusterAnn_Accession")
CGclusterAnn <- CGclusterAnn_Accession

data("methDNAclusters")
data("betaDiffclusters")
data("corrClusters")

integrCluster <- integrCluster(methDNAmatrix, betaDiffmatrix, corrMatrix, diffExprMatrix, CGclusterAnn, methDNAclusters, betaDiffclusters, corrClusters)

head(integrCluster)
}
