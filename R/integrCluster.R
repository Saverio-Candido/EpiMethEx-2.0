integrCluster <- function(methDNAmatrix, betaDiffmatrix, corrMatrix, diffExprMatrix, CGclusterAnn, methDNAclusters, betaDiffclusters, corrClusters){

  if(unique(CGclusterAnn$clusterType) == "RefGene_Accession"){"Valid CGclusterAnn"}else { stop("Invalid CGclusterAnn")}
  CGclusterAnn$clusterType <- NULL


  dfmeth <- methDNAclusters[!is.na(methDNAclusters$cgclusterName), ]
  dfcorr <- corrClusters[!is.na(corrClusters$corrClusterName), ]
  dfBetaDif <- betaDiffclusters[!is.na(betaDiffclusters$cgclusterName), ]
  AccessionList <- unique(CGclusterAnn$clusterName)

    count = 0
  clusterlist <- list()
  for(i in AccessionList){
    count <- count + 1
    n = 1
    methCluster <- unique(dfmeth[dfmeth$clusterName %in% i, c("cgclusterName")])
    corrCluster <- unique(dfcorr[dfcorr$clusterName %in% i, c("corrClusterName")])
    betadiffCluster <- unique(dfBetaDif[dfBetaDif$clusterName %in% i, c("cgclusterName")])

    if(length(na.omit(methCluster))>0 & length(na.omit(corrCluster))>0 & length(na.omit(betadiffCluster))>0){
      for (ii in methCluster){
        methVector <- dfmeth[dfmeth$cgclusterName %in% ii, c("ID")]
        for (iii in corrCluster){
          corrVector <- dfcorr[dfcorr$corrClusterName %in% iii,c("ID")]
          for (iiii in betadiffCluster){
            betadiffVector <- dfBetaDif[dfBetaDif$cgclusterName %in% iiii,c("ID")]
            interVect <- Reduce(intersect,list(methVector, corrVector, betadiffVector))
            if (length(interVect) >= 2){
              methType <- (strsplit(ii, "_"))
              methName <- sapply(methType, function(x) x[4])
              corrType <- (strsplit(iii, "_"))
              corrName <- sapply(corrType, function(x) x[2])
              betadiffType2 <- (strsplit(iiii, "_"))
              betadiffName2  <- sapply(betadiffType2, function(x) x[4:5])
              betadiffName  <- paste(betadiffName2, collapse = "_")

              SubsetAnn <- CGclusterAnn[CGclusterAnn$clusterName %in% i & CGclusterAnn$ID %in% interVect, ]

              clusterName <- paste0(unique(SubsetAnn$UCSC_RefGene_Name), "_",unique(SubsetAnn$UCSC_RefGene_Accession), "_", unique(SubsetAnn$UCSC_RefGene_Group), "_",unique(SubsetAnn$chr), ":",min(SubsetAnn$CGposition),"-", max(SubsetAnn$CGposition))
              VectName <- paste0(clusterName,"_", methName, "Meth_", corrName, "_", betadiffName,"_", n)
              SubsetAnn$integratedCluster <- VectName
              SubsetAnn$ClusterGenCoordinate <- paste0(unique(SubsetAnn$chr), ":",min(SubsetAnn$CGposition),"-", max(SubsetAnn$CGposition))

              Subsetmeth <- dfmeth[dfmeth$clusterName %in% i & dfmeth$ID %in% interVect,c("ID", "medianvalue")]
              colnames(Subsetmeth)[colnames(Subsetmeth) == "medianvalue"] <- "CG.median.methDNA"
              merged1 <- merge(SubsetAnn, Subsetmeth, by = "ID")
              Subsetcorr <- dfcorr[dfcorr$clusterName %in% i & dfcorr$ID %in% interVect, c("ID", "cor", "p") ]
              colnames(Subsetcorr)[colnames(Subsetcorr) == "cor"] <- "CG.median.corr"
              colnames(Subsetcorr)[colnames(Subsetcorr) == "p"] <- "CG.corr.pValue"
              merged2 <- merge(merged1, Subsetcorr, by = "ID")
              Subsetbetadiff <- dfBetaDif[dfBetaDif$clusterName %in% i & dfBetaDif$ID %in% interVect, c("ID", "beta_diff", "pValue") ]
              colnames(Subsetbetadiff)[colnames(Subsetbetadiff) == "beta_diff"] <- "CG.median.betaDiff"
              colnames(Subsetbetadiff)[colnames(Subsetbetadiff) == "pValue"] <- "CG.betaDiff.pValue"
              merged3 <- merge(merged2, Subsetbetadiff, by = "ID")
              clusterlist[[count]] <- merged3
            } # if length intersection vectors > 2
            n = n + 1
          } # for iiii
        } # for iii
      } # for ii
    } # if test number of clusters

    print (paste0("Processing: ", count, " of ", length(AccessionList)))
  }  #for i
  result <- do.call(rbind, clusterlist)

  result <- result %>%
    group_by(integratedCluster) %>%
    mutate(cluster.median.methDNA = median(CG.median.methDNA), cluster.median.corr = median(CG.median.corr), cluster.median.betaDiff = median(CG.median.betaDiff))

  diffExprMatrix <- diffExprMatrix[ , c(1, 2, 4, 5)]
  colnames(diffExprMatrix)[colnames(diffExprMatrix) == "MeanTumor"] <- "mean.expr.tumor"
  colnames(diffExprMatrix)[colnames(diffExprMatrix) == "logFC"] <- "logFC.expr.tumor.vs.normal"
  colnames(diffExprMatrix)[colnames(diffExprMatrix) == "pValue"] <- "pValue.logFC.expr"

  resultdef2 <- merge(result, diffExprMatrix, by = "UCSC_RefGene_Name", all.x = TRUE)
  resultdef <- merge(CGclusterAnn, resultdef2[ ,c(2, 16, 23:34)], by =c("clusterName", "ID"), all.x = TRUE)
  resultdef <- resultdef[ ,c(2:16, 1, 17:33)]
    }
