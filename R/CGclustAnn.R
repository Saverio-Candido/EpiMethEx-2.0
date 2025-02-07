CGclustAnn <- function(cgAnnotation, by = "genomePosition",   nCGclusters = 3, CGclusterSize = 200){
  CGclustAnn_Pos <- function(cgAnnotation,  nCGcluster = nCGclusters, clusterSize = CGclusterSize){

    chrlist<- intersect(c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21","chr22", "chrX", "chrY"), unique(cgAnnotation$chr))

    matrix1 <- cgAnnotation[!is.na(cgAnnotation$CGposition), ]
    matrix2 <- unique(matrix1[ , -c(4:6)])
    resultlist <- list()
    for (j in chrlist){

      df <- matrix2[matrix2$chr == j, ]
      cc <- df[order(df$CGposition,decreasing=FALSE),]

      datalist <- list()
      y <- 1
      z <- nrow(cc)
      i = 1
      for(i in y:z) {

        if (i == z) { t <- 0
        w <- -1  } else { t <- 1
        w <- clusterSize   }

        if (cc[i+t, 3] - cc[i, 3] > w) {
          if(y !=i){
            aa <- cc[y:i, ]
            aa$clusterName <- paste0(cc[i, 2],"_",c(i),"_",cc[i, 2],":", cc[y, 3],"-", cc[i, 3])
            aa$nCG <- nrow(cc [y:i, ])
            aa$size <- paste0(cc[i, 3]-cc[y, 3]+1)
            aa$chrcluster <- paste0(cc[i, 2])
            aa$startcluster <- paste0(cc[y, 3])
            aa$stopcluster <- paste0(cc[i, 3])
            aa$clusterType <- "genomePosition"
            if (nrow(aa) > nCGcluster) {
              datalist[[i]]<- aa
            }
          }
          y <- i+1
        }
        print(paste0("Processing: ", j, " - ", y, " of ", z))
        }
      result = do.call(rbind, datalist)
      rr <- result
      resultlist[[j]] <- rr
    }
    resultdef= do.call(rbind, resultlist)
    }
  CGclustAnn_RefSeqName <- function(cgAnnotation,  nCGcluster = nCGclusters, clusterSize = CGclusterSize){

    Genematrix2 <- cgAnnotation[!is.na(cgAnnotation$CGposition), ]
    Genematrix1 <- Genematrix2[Genematrix2$UCSC_RefGene_Name != ".", ]
    Genematrix <- unique(Genematrix1[, c(1:4, 7:14)])
    Genelist <- unique(Genematrix$UCSC_RefGene_Name)
    nn <- length(Genelist)

    resultlist <- list()

    for (j in Genelist){
      df <- Genematrix[Genematrix$UCSC_RefGene_Name == j, ]
      cc <- df[order(df$CGposition,decreasing=FALSE),]


      datalist <- list()
      y <- 1
      z <- nrow(cc)
      i = 1
      for(i in y:z) {

        if (i == z) { t <- 0
        w <- -1  } else { t <- 1
        w <- clusterSize }

        if (cc[i+t, 3] - cc[i, 3] > w) {
          if(y != i){
            aa <- cc[i:y, ]
            aa$clusterName <- paste0(cc[i, 4],"_",cc[i, 2],":",cc[y, 3],"-",cc[i, 3])
            aa$nCG <- nrow(cc[y:i, ])
            aa$size <- paste0(cc[i, 3]-cc[y, 3]+1)
            aa$chrcluster <- paste0(cc[i, 2])
            aa$startcluster <- paste0(cc[y, 3])
            aa$stopcluster <- paste0(cc[i, 3])
            aa$clusterType <- "RefGene_Name"
                        if (nrow(aa) > nCGcluster) {
              datalist[[i]]<- aa
            }
          }
          y <- i+1

        }
      }
      result= do.call(rbind, datalist)
      rr <- result
      resultlist[[j]] <- rr
      print (paste0("Processing: ", which(Genelist == j), " of ", nn))
    }
    resultdef= do.call(rbind, resultlist)
  }
  CGclustAnn_RefSeqAcc <- function(cgAnnotation,  nCGcluster = nCGclusters, clusterSize = CGclusterSize){

    NMmatrix <- cgAnnotation[!is.na(cgAnnotation$CGposition), ]
    NMmatrix$merge <- paste0(NMmatrix$UCSC_RefGene_Accession,"_",NMmatrix$UCSC_RefGene_Group)
    NMmatrix2 <- NMmatrix[NMmatrix$UCSC_RefGene_Name != ".", ]
    NMlist <- unique(NMmatrix2$merge)

    nn <- length(NMlist)

    resultlist <- list()

    for (j in NMlist){
      df <- NMmatrix2[NMmatrix2$merge == j, ]
      cc <- df[order(df$CGposition,decreasing=FALSE),]


      datalist <- list()
      y <- 1
      z <- nrow(cc)
      i = 1
      for(i in y:z) {

        if (i == z) { t <- 0
        w <- -1  } else { t <- 1
        w <- clusterSize }

        if (cc[i+t, 3] - cc[i, 3] > w) {
          if(y != i){
            aa <- cc[i:y, ]
            aa$clusterName <- paste0(cc[i, 15],"_",cc[i, 2],":",cc[y, 3],"-",cc[i, 3])
            aa$nCG <- nrow(cc[y:i, ])
            aa$size <- paste0(cc[i, 3]-cc[y, 3]+1)
            aa$chrcluster <- paste0(cc[i, 2])
            aa$startcluster <- paste0(cc[y, 3])
            aa$stopcluster <- paste0(cc[i, 3])
            aa$clusterType <- "RefGene_Accession"
            if (nrow(aa) > nCGcluster) {
              datalist[[i]]<- aa
            }
          }
          y <- i+1

        }
      }
      result= do.call(rbind, datalist)
      rr <- result
      resultlist[[j]] <- rr
      print (paste0("Processing: ", which(NMlist == j), " of ", nn))
    }
    resultdef= do.call(rbind, resultlist)
    }

  if (by ==  "genomePosition"){
    CGclustAnn <- CGclustAnn_Pos(cgAnnotation,nCGcluster = nCGclusters, clusterSize = CGclusterSize)
    }else if(by == "RefGene_Name"){
    CGclustAnn <- CGclustAnn_RefSeqName(cgAnnotation,nCGcluster = nCGclusters, clusterSize = CGclusterSize)
    }else if(by == "RefGene_Accession"){
    CGclustAnn <- CGclustAnn_RefSeqAcc(cgAnnotation,  nCGcluster = nCGclusters, clusterSize = CGclusterSize)
    } else { stop("Invalid input for by = ")}
}
