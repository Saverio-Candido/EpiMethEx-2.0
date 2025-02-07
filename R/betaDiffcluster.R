betaDiffcluster <- function(cgAnnotation, CGclusterAnn, methDNAmatrix, betaDiffmatrix, by = "genomePosition", nCGclusters = 1){

  betaDiffclust_Pos <- function(cgAnnotation, CGclusterAnn, methDNAmatrix, betaDiffmatrix, by = "genomePosition", nCGcluster = nCGclusters){

    if(unique(CGclusterAnn$clusterType) == "genomePosition"){"Valid CGclusterAnn"}else { stop("Invalid CGclusterAnn")}
    CGclusterAnn$clusterType <- NULL

    chrlist <- unique(CGclusterAnn$clusterName)

    chrmethclu <- list()
    for (j in chrlist) {

      df <- CGclusterAnn[CGclusterAnn$clusterName == j, c(1:3, 12:17)]
      df<- df[order(df$CGposition),]
      df$progr <- 1:nrow(df)
      CGlist <- df$ID
      subBetadiff <- betaDiffmatrix[betaDiffmatrix$ID %in% CGlist,]
      if (sum(!is.na(subBetadiff$beta_diff))>=2){
        cc <- merge(df, subBetadiff[, c(1,2)], by ="ID", all = TRUE)

        cc$methDiffType <- ifelse(cc$beta_diff > 0.5, "Strong_meth",
                                  ifelse(cc$beta_diff > 0 & cc$beta_diff <= 0.5, "Weak_meth",
                                         ifelse(cc$beta_diff < -0.5, "Strong_de-meth",
                                                ifelse(cc$beta_diff >= -0.5 & cc$beta_diff < 0, "Weak_de-meth",""))))

        qq <- cc[complete.cases(cc$beta_diff), ]
        methDiffType <- unique(qq$methDiffType)

        resultlist <- list()
        for (meth in methDiffType){

          if (meth == "Strong_meth"){
            uu <- qq[qq$methDiffType == "Strong_meth", ]
          } else if(meth == "Weak_meth"){
            uu <- qq[qq$methDiffType == "Weak_meth", ]
          }  else if(meth == "Weak_de-meth"){
            uu <- qq[qq$methDiffType == "Weak_de-meth", ]
          } else if(meth == "Strong_de-meth") {
            uu <- qq[qq$methDiffType == "Strong_de-meth", ]
          } else {
            uu <- data.frame()
          }

          if(nrow(uu) > 1){
            uu <- uu[order(uu$CGposition),]

            datalist <- list()
            y <- 1
            z <- nrow(uu)
            i <- 1
            for(i in y:z) {

              if (i == z) { t <- 0
              w <- 0  } else { t <- 1
              w <- 2   }
              if (uu[i+t, c("progr")] - uu[i, c("progr")] >= w) {
                if(y!=i){
                  aa <- uu[y:i, ]
                  aa$cgclusterName <- paste0(sapply(strsplit(uu[i, 4], "_"), function(x) paste0(x[1],"_", x[2])), "_", c(meth), "_", uu[1,2],":",uu[y,3], "-", uu[i,3])
                  aa$ClusterCoverage <- c(nrow(uu[i:y, ]))/aa$nCG
                  aa$cgN <- nrow(uu[i:y, ])
                  aa$CGsize <- uu[i, 3]-uu[y, 3]+1
                  aa$CGstartcluster <- uu[y, 3]
                  aa$CGstopcluster <- uu[i, 3]
                  annloop <- cgAnnotation[cgAnnotation$ID %in% aa$ID, ]
                  aa <- aa %>% mutate(BetaDiff.median.cluster = median(beta_diff), CG.IQR.cluster = IQR(beta_diff),
                                      gene = paste0(unique(annloop$UCSC_RefGene_Name), collapse = ", "))

                  if (nrow(aa) > nCGcluster) {datalist[[i]]<- aa}
                }
                y <- i+1
              }
            }
            result= do.call(rbind, datalist)
            resultlist[[meth]] <- result
          }
        }
        resultdef= do.call(rbind, resultlist)
        chrmethclu[[j]] <- resultdef
        print(j)
      }
    }
    resultdefdef = do.call(rbind, chrmethclu)
    resultdefdef$progr <- NULL
    finalDF <- left_join(CGclusterAnn[ , c(1:3, 12:17)], methDNAmatrix, by = "ID")
    finalDF2 <- left_join(finalDF, betaDiffmatrix, by = "ID")
    finalDF3 <- left_join(finalDF2, resultdefdef[ , c(1, 12:20)], by = "ID")
  }
  betaDiffclust_Acc <- function(cgAnnotation, CGclusterAnn, methDNAmatrix, betaDiffmatrix, by = "RefGene_Accession", nCGcluster = nCGclusters){

    if(unique(CGclusterAnn$clusterType) == "RefGene_Accession"){"Valid CGclusterAnn"}else { stop("Invalid CGclusterAnn")}
    CGclusterAnn$clusterType <- NULL

    chrlist <- unique(CGclusterAnn$clusterName)

    chrmethclu <- list()
    for (j in chrlist) {

      df <- CGclusterAnn[CGclusterAnn$clusterName == j, ]
      df <- df[order(df$CGposition),]
      df$progr <- 1:nrow(df)
      CGlist <- df$ID
      subBetadiff <- betaDiffmatrix[betaDiffmatrix$ID %in% CGlist,]
      if (sum(!is.na(subBetadiff$beta_diff))>=2){
        cc <- merge(df, subBetadiff[, c(1,2)], by ="ID", all = TRUE)

        cc$methDiffType <- ifelse(cc$beta_diff >= 0.5, "Strong_meth",
                                  ifelse(cc$beta_diff >= 0.1 & cc$beta_diff < 0.5, "Weak_meth",
                                         ifelse(cc$beta_diff <= -0.5, "Strong_de-meth",
                                                ifelse(cc$beta_diff > -0.5 & cc$beta_diff <= -0.1, "Weak_de-meth",""))))

        qq <- cc[complete.cases(cc$beta_diff), ]
        methDiffType <- unique(qq$methDiffType)

        resultlist <- list()
        for (meth in methDiffType){
          meth <- methDiffType[1]
          if (meth == "Strong_meth"){
            uu <- qq[qq$methDiffType == "Strong_meth", ]
          } else if(meth == "Weak_meth"){
            uu <- qq[qq$methDiffType == "Weak_meth", ]
          }  else if(meth == "Weak_de-meth"){
            uu <- qq[qq$methDiffType == "Weak_de-meth", ]
          } else if(meth == "Strong_de-meth") {
            uu <- qq[qq$methDiffType == "Strong_de-meth", ]
          } else {
            uu <- data.frame()
          }

          if(nrow(uu) > 1){
            uu <- uu[order(uu$CGposition),]

            datalist <- list()
            y <- 1
            z <- nrow(uu)
            i <- 1
            for(i in y:z) {

              if (i == z) { t <- 0
              w <- 0  } else { t <- 1
              w <- 2   }
              if (uu[i+t, c("progr")] - uu[i, c("progr")] >= w) {
                if(y!=i){
                  aa <- uu[y:i, ]
                  aa$cgclusterName <- paste0(uu[i,c("merge")], "_", uu[i,c("methDiffType")], "_", uu[i,c("chr")], ":", uu[y, c("CGposition")], "-", uu[i, c("CGposition")])
                  aa$ClusterCoverage <- c(nrow(uu[y:i, ]))/aa$nCG
                  aa$cgN <- nrow(uu[y:i, ])
                  aa$CGsize <- uu[i, c("CGposition")]-uu[y, c("CGposition")]+1
                  aa$CGstartcluster <- uu[y, c("CGposition")]
                  aa$CGstopcluster <- uu[i, c("CGposition")]
                  aa <- aa %>% mutate(BetaDiff.median.cluster = median(beta_diff),
                                      CG.IQR.cluster = IQR(beta_diff))

                  if (nrow(aa) > 1) {datalist[[i]]<- aa}
                }
                y <- i+1
              }
            }
            result= do.call(rbind, datalist)
            resultlist[[meth]] <- result
          }
        }
        resultdef= do.call(rbind, resultlist)
        chrmethclu[[j]] <- resultdef
        print(j)
      }
    }
    resultdefdef = do.call(rbind, chrmethclu)
    resultdefdef$progr <- NULL
    finalDF <- left_join(CGclusterAnn[ ,c(1:6, 15:21)],methDNAmatrix, by = "ID")
    finalDF$mergeIDgene <- paste0(finalDF$ID, "_", finalDF$merge)
    resultdefdef$mergeIDgene <- paste0(resultdefdef$ID, "_", resultdefdef$merge)

    finalDF2 <- left_join(finalDF, betaDiffmatrix, by = "ID")
    finalDF3 <- left_join(finalDF2, resultdefdef[ , c(32, 24:31)], by = ("mergeIDgene"))
    finalDF3 <- finalDF3[ ,-c(7,15)]
  }

  if(by == "genomePosition"){
    results <- betaDiffclust_Pos(cgAnnotation, CGclusterAnn, methDNAmatrix, betaDiffmatrix, by = "genomePosition", nCGcluster = nCGclusters)
  } else if(by == "RefGene_Accession"){
    results <- betaDiffclust_Acc(, CGclusterAnn, methDNAmatrix, betaDiffmatrix, by = "RefGene_Accession", nCGcluster = nCGclusters)
  } else { stop("Invalid input for by = ")}
}
