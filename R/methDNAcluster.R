methDNAcluster <- function(cgAnnotation, CGclusterAnn, methDNAmatrix, by = "genomePosition", nCGclusters = 1){

  methDNAcluster_Pos <- function(cgAnnotation, CGclusterAnn, methDNAmatrix, by = "genomePosition", nCGcluster = nCGclusters){

    if(unique(CGclusterAnn$clusterType) == "genomePosition"){"Valid CGclusterAnn"}else { stop("Invalid CGclusterAnn")}
    CGclusterAnn$clusterType <- NULL

    chrlist <- unique(CGclusterAnn$clusterName)

    chrmethclu <- list()
    for (j in chrlist) {

      df <- CGclusterAnn[CGclusterAnn$clusterName == j, c(1:3, 12:17)]
      df$progr <- 1:nrow(df)
      CGlist <- df$ID
      medianValue <- methDNAmatrix[methDNAmatrix$ID %in% CGlist,]
      if (sum(!is.na(medianValue$medianvalue))>=2){
        cc <- merge(df, medianValue, by ="ID", all = TRUE)
        cc$methType <- ifelse(cc$medianvalue > 0.6, "hyper",
                              ifelse(cc$medianvalue >= 0.2 & cc$medianvalue <= 0.6, "partially",
                                     ifelse(cc$medianvalue > 0 & cc$medianvalue < 0.2, "hypo","")))
        qq <- cc[complete.cases(cc$medianvalue), ]
        methType <- unique(qq$methType)

        resultlist <- list()
        for (meth in methType){

          if (meth == "hyper"){
            uu <- qq[qq$methType == "hyper", ]
          } else if(meth == "hypo"){
            uu <- qq[qq$methType == "hypo", ]
          } else if(meth == "partially") {
            uu <- qq[qq$methType == "partially", ]
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
                  annloop <- cgAnnotation[cgAnnotation$ID %in% aa$ID, c(1,4)]
                  aa <- aa %>% mutate(CG.median.cluster = median(medianvalue), CG.IQR.cluster = IQR(medianvalue),
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
        print(paste0("Processing: ", j))
      }
    }
    resultdefdef = do.call(rbind, chrmethclu)
    resultdefdef$progr <- NULL
    finalDF <- left_join(CGclusterAnn[ , c(1:3, 12:17)],methDNAmatrix, by = "ID")
    finalDF2 <- left_join(finalDF, resultdefdef[ , c(1, 12:20)], by = "ID")
  }
  methDNAcluster_RefGeneAcc <- function(CGclusterAnn, methDNAmatrix, by = "RefGene_Accession", nCGcluster = nCGclusters){

    if(unique(CGclusterAnn$clusterType) == "RefGene_Accession"){"Valid CGclusterAnn"}else { stop("Invalid CGclusterAnn")}
    CGclusterAnn$clusterType <- NULL

    chrlist <- unique(CGclusterAnn$clusterName)
    chrmethclu <- list()

    for (j in chrlist) {
      df <- CGclusterAnn[CGclusterAnn$clusterName == j, ]
      CGlist <- df$ID
      medianValue <-  methDNAmatrix[ methDNAmatrix$ID %in% CGlist,]
      if (sum(!is.na(medianValue$medianvalue))>=2){
        cc <- merge(df, medianValue, by ="ID", all = TRUE)
        cc <-  cc[order(cc$CGposition), ]
        cc$methType <- ifelse(cc$medianvalue > 0.6, "hyper",
                              ifelse(cc$medianvalue >= 0.2 & cc$medianvalue <= 0.6, "partially",
                                     ifelse(cc$medianvalue > 0 & cc$medianvalue < 0.2, "hypo","")))
        cc<- cbind(prog = seq_len(nrow(cc)), cc)
        qq <- subset(cc, !is.na(medianValue))
        qq <- cc[complete.cases(cc$medianvalue), ]
        methType <- unique(qq$methType)
        resultlist <- list()

        for (meth in methType){
          uu <- qq[qq$methType %in% meth, ]
          if(nrow(uu) > 1){
            datalist <- list()
            y <- 1
            z <- nrow(uu)
            i <- 1
            for(i in y:z) {
              if (i == z) { t <- 0
              w <- 0  } else { t <- 1
              w <- 2   }
              if (uu[i+t, c("prog")] - uu[i, c("prog")] >= w) {
                aa <- uu[y:i, ]
                aa$cgclusterName <- paste0(uu[i,c("merge")], "_", uu[i,c("methType")], "_", uu[i,c("chr")], ":", uu[y, c("CGposition")], "-", uu[i, c("CGposition")])
                aa$ClusterCoverage <- c(nrow(uu[y:i, ]))/aa$nCG
                aa$cgN <- nrow(uu[y:i, ])
                aa$CGsize <- uu[i, c("CGposition")]-uu[y, c("CGposition")]+1
                aa$CGstartcluster <- uu[y, c("CGposition")]
                aa$CGstopcluster <- uu[i, c("CGposition")]
                aa <- aa %>% mutate(CG.median.cluster = median(medianvalue),
                                    CG.IQR.cluster = IQR(medianvalue))
                if (nrow(aa) > nCGcluster) {
                  datalist[[i]]<- aa
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
        print(paste0("Processing: ", j))
      }
    }
    resultdefdef = do.call(rbind, chrmethclu)
    resultdefdef$prog <- NULL

    finalDF <- left_join(CGclusterAnn[ , c(1:6, 15:21)], methDNAmatrix, by = "ID")
    finalDF$mergeIDgene <- paste0(finalDF$ID, "_", finalDF$merge)
    resultdefdef$mergeIDgene <- paste0(resultdefdef$ID, "_", resultdefdef$merge)
    finalDF2 <- left_join(finalDF, resultdefdef[ , c(32, 24:31)], by = ("mergeIDgene"))
    finalDF2 <- finalDF2[ ,-c(7, 15)]
  }

  if(by == "genomePosition"){
    results <- methDNAcluster_Pos(cgAnnotation, CGclusterAnn, methDNAmatrix, by = "genomePosition", nCGcluster = nCGclusters)
  } else if(by == "RefGene_Accession"){
    results <-  methDNAcluster_RefGeneAcc(CGclusterAnn, methDNAmatrix, by = "RefGene_Accession", nCGcluster = nCGclusters)
  } else { stop("Invalid input for by = ")}
}
