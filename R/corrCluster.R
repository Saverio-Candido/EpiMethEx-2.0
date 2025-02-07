corrCluster <- function(cgAnnotation, CGclusterAnn, corrMatrix, by = "RefGene_Name", nCGclusters = 1){

  corrCluster_Name <-  function(cgAnnotation, CGclusteAnn, corrMatrix, by = "RefGene_Name", nCGcluster = 1){

    if(unique(CGclusterAnn$clusterType) == "RefGene_Name"){"Valid CGclusterAnn"}else { stop("Invalid CGclusterAnn")}
    CGclusterAnn$clusterType <- NULL

    allgeneCG <- distinct(cgAnnotation[cgAnnotation$UCSC_RefGene_Name != ".", c(1,4)])
    geneCG <- subset(allgeneCG, UCSC_RefGene_Name %in% unique(corrMatrix$UCSC_RefGene_Name))
    RefGeneclustlist <- subset(CGclusteAnn, ID %in% unique(corrMatrix$ID))
    RefGeneclustlist <- unique(RefGeneclustlist$clusterName)

    # loop 1 by RefGenecluster
    RefGeneClusterResult <- list()
    for (j in RefGeneclustlist) {
      annmatrix <- CGclusteAnn[CGclusteAnn$clusterName == j, ]
      genelist <- geneCG[geneCG$ID %in% annmatrix$ID, ]
      genelist <- unique(genelist$UCSC_RefGene_Name)

      # loop2: selection by gene within loop1
      resultbyGene <- list()
      for (ge in genelist){
        genematrix <- corrMatrix[corrMatrix$UCSC_RefGene_Name %in% ge, ]
        subresult <- merge(annmatrix, genematrix, by = c("ID","UCSC_RefGene_Name"),  all.x = TRUE)
        subresult <-  subresult[order(subresult$CGposition), ]
        subresult$corrType <- ifelse(subresult$cor > 0.7, "strongPosCorr",
                                     ifelse(subresult$cor >= 0.3 & subresult$cor <= 0.7, "moderatePosCorr",
                                            ifelse(subresult$cor > 0 &subresult$cor < 0.3, "weakPosCorr",
                                                   ifelse(subresult$cor < -0.7, "strongNegCorr",
                                                          ifelse(subresult$cor >= - 0.7  & subresult$cor <= -0.3, "moderateNegCorr",
                                                                 ifelse(subresult$cor > -0.3 & subresult$cor < 0, "weakNegCorr", ""))))))

        subresult<- cbind(prog = seq_len(nrow(subresult)), subresult)
        qq <- subset(subresult, !is.na(cor))

        corrType <- unique(qq$corrType)

        # loop 3: selection by corcordance corrValues within loop 1 and 2 selection
        resultbyCorr <- list()
        for (cor in corrType){
          uu <- qq[qq$corrType %in% cor, ]
          if(nrow(uu) > 1){

            # loop4: loop for position clustering
            datalist <- list()
            y <- 1
            z <- nrow(uu)
            i = 1
            for(i in y:z) {
              if(i == z){
                break
              }
              if (uu[i+1, 1] - uu[i, 1] > 1) {
                aa <- uu[y:i, ]
                aa$corrClusterName <- paste0(ge, "_", cor, "_", uu[y, 4],":",uu[y, 5],"-", uu[i, 5])
                aa$ClusterCoverage <- c(nrow(uu[y:i, ]))/aa$nCG
                aa$cgN <- nrow(uu[y:i, ])
                aa$CGsize <- uu[i, 5]-uu[y, 5]+1
                aa$CGstartcluster <- uu[y, 5]
                aa$CGstopcluster <- uu[i, 5]
                #aa$gene <- ge
                annloop <- cgAnnotation[cgAnnotation$ID %in% aa$ID, c(1, 5, 6)]
                annloop$refGeneIso <- paste0(annloop$UCSC_RefGene_Accession, "_", annloop$UCSC_RefGene_Group)
                aa <- aa %>% mutate(Corr.median.cluster = median(cor), Corr.IQR.cluster = IQR(cor),
                                    refGeneIso = paste0(unique(annloop$refGeneIso), collapse = ", "))

                if (nrow(aa) > nCGcluster) {
                  datalist[[i]]<- aa
                }
                y <- i+1
              }
            }
            #end loop clustering
            result= do.call(rbind, datalist)
            # end loop4
            resultbyCorr[[cor]] <- result
          }
        }
        resultcorrelation = do.call(rbind, resultbyCorr)
        #end  of loop3
        ww <- resultcorrelation
        resultbyGene[[ge]] <- ww
      }
      resultLoop2 = do.call(rbind, resultbyGene)
      #end of loop2:
      loop1 <- resultLoop2
      RefGeneClusterResult[[j]] <- loop1
      print(j)
    }
    # end of loop1
    resultdef <- do.call(rbind, RefGeneClusterResult)

    finalDF <- left_join(CGclusteAnn[ , c(1:4, 13:18)],corrMatrix, by = c("ID","UCSC_RefGene_Name"))
    finalDF2 <- left_join(finalDF, resultdef[ , c(2, 3, 14, 23:31)], by = c("ID","UCSC_RefGene_Name", "clusterName"))

  }
  corrCluster_Accession <- function(cgAnnotation, CGclusterAnn, corrMatrix, by = "RefGene_Accession", nCGcluster = nCGclusters){

    if(unique(CGclusterAnn$clusterType) == "RefGene_Accession"){"Valid CGclusterAnn"}else { stop("Invalid CGclusterAnn")}
    CGclusterAnn$clusterType <- NULL

    allgeneCG <- unique(cgAnnotation[cgAnnotation$UCSC_RefGene_Name != ".", c(1,4)])

    geneCG <- subset(allgeneCG, UCSC_RefGene_Name %in% unique(corrMatrix$UCSC_RefGene_Name))
    NMclustlist <- subset(CGclusterAnn, ID %in% unique(corrMatrix$ID))
    NMclustlist <- unique(NMclustlist$clusterName)

    # loop 1 by chrcluster
    chrClusterResult <- list()
    for (j in NMclustlist) {

      annmatrix <- CGclusterAnn[CGclusterAnn$clusterName == j, ]
      genelist <- unique(annmatrix$UCSC_RefGene_Name)

      # loop2: selection by gene within loop1 (chrcluster)
      resultbyGene <- list()
      for (ge in genelist){
        genematrix <- corrMatrix[corrMatrix$UCSC_RefGene_Name %in% ge, c(1,3,4) ]
        subresult <- merge(annmatrix, genematrix, by = "ID", all.x = TRUE)
        subresult <-  subresult[order(subresult$CGposition), ]
        subresult$corrType <- ifelse(subresult$cor > 0.7, "strongPosCorr",
                                     ifelse(subresult$cor >= 0.3 & subresult$cor <= 0.7, "moderatePosCorr",
                                            ifelse(subresult$cor > 0 &subresult$cor < 0.3, "weakPosCorr",
                                                   ifelse(subresult$cor < -0.7, "strongNegCorr",
                                                          ifelse(subresult$cor >= - 0.7  & subresult$cor <= -0.3, "moderateNegCorr",
                                                                 ifelse(subresult$cor > -0.3 & subresult$cor < 0, "weakNegCorr", ""))))))

        subresult<- cbind(prog = seq_len(nrow(subresult)), subresult)
        qq <- subset(subresult, !is.na(cor))
        corrType <- unique(qq$corrType)

        # loop 3: selection by corcordance corrValues within loop 1 and 2 selection
        resultbyCorr <- list()
        for (cor in corrType){
          uu <- qq[qq$corrType %in% cor, ]
          if(nrow(uu) > 1){

            # loop4: loop for position clustering
            datalist <- list()
            y <- 1
            z <- nrow(uu)
            i = 1
            for(i in y:z) {
              if (i == z) { t <- 0
              w <- 0  } else { t <- 1
              w <- 2   }
              if (uu[i+t, 1] - uu[i, 1] >= w) {
                aa <- uu[y:i, ]
                aa$corrClusterName <- paste0(ge, "_", cor, "_", uu[y, 3],":",uu[y, 4],"-", uu[i, 4])
                aa$ClusterCoverage <- c(nrow(uu[y:i, ]))/aa$nCG
                aa$cgN <- nrow(uu[y:i, ])
                aa$CGsize <- uu[i, 4]-uu[y, 4]+1
                aa$CGstartcluster <- uu[y, 4]
                aa$CGstopcluster <- uu[i, 4]
                aa <- aa %>% mutate(Corr.median.cluster = median(cor), Corr.IQR.cluster = IQR(cor))

                if (nrow(aa) > nCGcluster) {
                  datalist[[i]]<- aa
                }
                y <- i+1
              }
            }
            #end loop clustering
            result= do.call(rbind, datalist)
            # end loop4
            resultbyCorr[[cor]] <- result
          }
        }
        resultcorrelation = do.call(rbind, resultbyCorr)
        #end  of loop3
        ww <- resultcorrelation
        resultbyGene[[ge]] <- ww
      }
      resultLoop2 = do.call(rbind, resultbyGene)
      #end of loop2:
      loop1 <- resultLoop2
      chrClusterResult[[j]] <- loop1
      print(paste0("Processing: ", j))
    }
    # end of loop1

    resultdef <- do.call(rbind, chrClusterResult)

    finalDF <- left_join(CGclusterAnn[ , c(1:6, 16:21)],corrMatrix, by = c("ID","UCSC_RefGene_Name"))
    finalDF2 <- left_join(finalDF, resultdef[ , c(2, 5:7, 26:33)], by = c("ID","UCSC_RefGene_Name", "UCSC_RefGene_Accession","UCSC_RefGene_Group"))
  }

  if(by == "RefGene_Name"){
    results <- corrCluster_Name(cgAnnotation, CGclusterAnn, corrMatrix, by = "RefGene_Name", nCGcluster = nCGclusters)
  } else if(by == "RefGene_Accession"){
    results <- corrCluster_Accession(cgAnnotation, CGclusterAnn, corrMatrix, by = "RefGene_Accession", nCGcluster = nCGclusters)
  } else { stop("Invalid input for by = ")}
}
