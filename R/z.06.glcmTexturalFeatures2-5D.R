#' GLCM2.5D  features
#'
#' @description  Extract the GLCM2.5D features
#' @import radiomics data.table plyr reshape2
#' @export
glcmTexturalFeatures25D <- function(imgObj,n_grey){

  
  if(length(dim(imgObj))<3) { 
    imgObj <- array(imgObj, dim=c(dim(imgObj)[1],dim(imgObj)[2],1))
  } 
  
  # compute number of non-NA voxels
  
  nVoxel <- dim(imgObj)[1]*dim(imgObj)[2]*dim(imgObj)[3] - sum(is.na(imgObj))
  grayLevlesMax <- max(imgObj,na.rm = T)
  ### compute Gray Levels Cooccurrence Matrices
  
  G_list_0 <- list()
  G_list_45 <- list()
  G_list_90 <- list()
  G_list_135 <- list()
  

  
  #Compute grey level cooccurrence matrices for 4 different directions within each slice
  for (i in 1:dim(imgObj)[3]){
    if(length(imgObj[,,i])*.005 >= sum(!is.na(imgObj[,,i]))) next
    imgObj[,,i] <- round(imgObj[,,i])
    #if(length(table(unique(imgObj[,,i])))<n_grey) n_grey <- length(table(unique(imgObj[,,i])))-1
    if (min(imgObj[,,i],na.rm = T) < 0) {imgObj[,,i] <- imgObj[,,i] + abs(min(imgObj[,,i],na.rm = T))}
    G_list_0[[i]] <- as.matrix(glcm(imgObj[,,i], angle = 0, normalize = F,verbose=F,n_grey = n_grey))
    G_list_45[[i]] <- as.matrix(glcm(imgObj[,,i], angle = 45, normalize = F,verbose=F,n_grey = n_grey))
    G_list_90[[i]] <- as.matrix(glcm(imgObj[,,i], angle = 90, normalize = F,verbose=F,n_grey = n_grey))
    G_list_135[[i]] <- as.matrix(glcm(imgObj[,,i], angle = 135, normalize = F,verbose=F,n_grey = n_grey))
  }
  
  #elimino gli elementi NULL della lista
  if(length(G_list_0)!=0){
    if(length(which(sapply(G_list_0, is.null)))!=0){
      G_list_0 = G_list_0[-which(sapply(G_list_0, is.null))]
    }
  }
  
  #elimino gli elementi NULL della lista
  if(length(G_list_45)!=0){
    if(length(which(sapply(G_list_45, is.null)))!=0){
      G_list_45 = G_list_45[-which(sapply(G_list_45, is.null))]
    }
  }
  
  #elimino gli elementi NULL della lista
  if(length(G_list_90)!=0){
    if(length(which(sapply(G_list_90, is.null)))!=0){
      G_list_90 = G_list_90[-which(sapply(G_list_90, is.null))]
    }
  }
  
  #elimino gli elementi NULL della lista
  if(length(G_list_135)!=0){
    if(length(which(sapply(G_list_135, is.null)))!=0){
      G_list_135 = G_list_135[-which(sapply(G_list_135, is.null))]
    }
  }
  
  sumtot <- list()
  ### DIRECTION 0
  matrix.df <- ldply(G_list_0, data.table::melt, varnames=c("row", "col"))
  sumtot[[1]] <- acast(matrix.df, row ~ col, sum)
  ### DIRECTION 45
  matrix.df <- ldply(G_list_45, data.table::melt, varnames=c("row", "col"))
  sumtot[[2]] <- acast(matrix.df, row ~ col, sum)
  ### DIRECTION 90
  matrix.df <- ldply(G_list_90, data.table::melt, varnames=c("row", "col"))
  sumtot[[3]] <- acast(matrix.df, row ~ col, sum)
  ### DIRECTION 135
  matrix.df <- ldply(G_list_135, data.table::melt, varnames=c("row", "col"))
  sumtot[[4]] <- acast(matrix.df, row ~ col, sum)
  
  
  
  #Initialise data table for storing GLCM features; I have added a few
  featNames <- c("F_cm_2.5D.joint.max", "F_cm_2.5D.joint.avg", "F_cm_2.5D.joint.var", "F_cm_2.5D.joint.entr",
                 "F_cm_2.5D.diff.avg", "F_cm_2.5D.diff.var", "F_cm_2.5D.diff.entr",
                 "F_cm_2.5D.sum.avg", "F_cm_2.5D.sum.var", "F_cm_2.5D.sum.entr", "F_cm_2.5D.energy","F_cm_2.5D.contrast","F_cm_2.5D.dissimilarity",
                 "F_cm_2.5D.inv.diff","F_cm_2.5D.inv.diff.norm","F_cm_2.5D.inv.diff.mom","F_cm_2.5D.inv.diff.mom.norm","F_cm_2.5D.inv.var",
                 "F_cm_2.5D.corr","F_cm_2.5D.auto.corr","F_cm_2.5D.clust.tend","F_cm_2.5D.clust.shade","F_cm_2.5D.clust.prom",
                 "F_cm_2.5D.info.corr.1","F_cm_2.5D.info.corr.2")
  F_cm <- data.table(matrix(NA, nrow=length(sumtot), ncol=length(featNames)))
  colnames(F_cm) <- featNames
  
  for (iter in seq_len(length(sumtot))){
    
    Ng <- ncol(sumtot[[iter]])
    
    #Convert matrix to data table
    df.G    <- data.table(sumtot[[iter]])
    
    #Add row grey level intensity
    df.G$i <- as.numeric(row.names(sumtot[[iter]]))
    
    #Convert from wide to long table. This is the preferred format for data tables and data frames
    df.G   <- data.table::melt(df.G, id.vars="i", variable.name="j", value.name="n", variable.factor=FALSE)
    
    #Convert j from string to numeric
    df.G$j <- as.numeric(df.G$j)
    
    #Remove combinations with 0 counts
    df.G   <- df.G[n>0,]
    
    #Convert Grey level coccurrence matrix to joint probability
    df.p_ij <- df.G[,.(p_ij=n/sum(df.G$n)), by=.(i,j)] #joint probability
    df.p_i  <- df.p_ij[,.(p_i=sum(p_ij)), by=i]        #marginal probability over columns
    df.p_j  <- df.p_ij[,.(p_j=sum(p_ij)), by=j]        #marginal probability over rows
    
    #Diagonal probabilities (p(i-j))
    #First, we create a new column k which contains the absolute value of i-j.
    #Second, we sum the joint probability where k is the same.
    #This can written as one line by chaining the operations.
    df.p_imj <- copy(df.p_ij)
    df.p_imj <- df.p_imj[,"k":=abs(i-j)][,.(p_imj=sum(p_ij)), by=k]
    
    #Cross-diagonal probabilities (p(i+j))
    #Again, we first create a new column k which contains i+j
    #Second, we sum the probability where k is the same.
    #This is written in one line by chaining the operations.
    df.p_ipj <- copy(df.p_ij)
    df.p_ipj <- df.p_ipj[,"k":=i+j][,.(p_ipj=sum(p_ij)), by=k]
    
    #Merger of df.p_ij, df.p_i and df.p_j
    df.p_ij <- merge(x=df.p_ij, y=df.p_i, by="i")
    df.p_ij <- merge(x=df.p_ij, y=df.p_j, by="j")
    
    #Thus we have five probability matrices
    #Joint probability:          df.p_ij with probabilities p_ij, p_i and p_j, and indices i, j
    #Marginal probability:       df.p_i with probability p_i, and index i
    #Marginal probability:       df.p_j with probability p_j, and index j
    #Diagonal probability:       df.p_imj with probability p_imj and index k
    #Cross-diagonal probability: df.p_ipj with probability p_ipj and index k
    
    #Calculate features
    
    #Joint maximum
    F_cm$F_cm_2.5D.joint.max[iter]         <- max(df.p_ij$p_ij)
    
    #Joint average
    F_cm$F_cm_2.5D.joint.avg[iter]         <- sum(df.p_ij$i * df.p_ij$p_ij)
    
    #Joint variance
    mu <- sum(df.p_ij$i * df.p_ij$p_ij)
    F_cm$F_cm_2.5D.joint.var[iter]         <- sum((df.p_ij$i-mu)^2 * df.p_ij$p_ij)
    
    #Joint entropy
    F_cm$F_cm_2.5D.joint.entr[iter]        <- - sum(df.p_ij$p_ij * log2(df.p_ij$p_ij))
    
    #Difference average
    F_cm$F_cm_2.5D.diff.avg[iter]          <- sum((df.p_imj$k) * df.p_imj$p_imj)
    
    #Difference variance
    mu <- sum((df.p_imj$k) * df.p_imj$p_imj)
    F_cm$F_cm_2.5D.diff.var[iter]          <- sum((df.p_imj$k-mu)^2 * df.p_imj$p_imj)
    
    #Difference entropy
    F_cm$F_cm_2.5D.diff.entr[iter]         <- - sum(df.p_imj$p_imj * log2(df.p_imj$p_imj))
    
    #Sum average
    F_cm$F_cm_2.5D.sum.avg[iter]           <- sum(df.p_ipj$k * df.p_ipj$p_ipj)
    
    #Sum variance
    mu <- sum(df.p_ipj$k * df.p_ipj$p_ipj)
    F_cm$F_cm_2.5D.sum.var[iter]           <- sum((df.p_ipj$k-mu)^2 * df.p_ipj$p_ipj)
    
    #Sum entropy
    F_cm$F_cm_2.5D.sum.entr[iter]          <- - sum(df.p_ipj$p_ipj * log2(df.p_ipj$p_ipj))
    
    #Angular second moment
    F_cm$F_cm_2.5D.energy[iter]          <- sum(df.p_ij$p_ij^2)
    
    #Contrast
    F_cm$F_cm_2.5D.contrast[iter]          <- sum((df.p_ij$i - df.p_ij$j)^2 * df.p_ij$p_ij)
    
    #Dissimilarity
    F_cm$F_cm_2.5D.dissimilarity[iter]          <- sum(abs(df.p_ij$i - df.p_ij$j) * df.p_ij$p_ij)
    
    #Inverse difference
    F_cm$F_cm_2.5D.inv.diff[iter]          <- sum( df.p_ij$p_ij/(1+abs(df.p_ij$i - df.p_ij$j)) )
    
    #Inverse difference normalized
    D <- abs(df.p_ij$i - df.p_ij$j)/grayLevlesMax
    F_cm$F_cm_2.5D.inv.diff.norm[iter]          <- sum(df.p_ij$p_ij / (1+D) )
    
    #Inverse difference moment
    F_cm$F_cm_2.5D.inv.diff.mom[iter]          <- sum(df.p_ij$p_ij/(1+(df.p_ij$i - df.p_ij$j)^2))
    
    #Inverse difference moment normalized
    DD <- ((df.p_ij$i - df.p_ij$j)/grayLevlesMax)^2
    F_cm$F_cm_2.5D.inv.diff.mom.norm[iter]          <- sum(df.p_ij$p_ij / (1+DD))
    
    #Inverse variance
    df.jmorethani <- df.p_ij[which(df.p_ij$j>df.p_ij$i),]
    F_cm$F_cm_2.5D.inv.var[iter]          <- 2*sum(df.jmorethani$p_ij/(df.jmorethani$i-df.jmorethani$j)^2)
    
    #Correlation
    mu.i <- sum(df.p_i$i * df.p_i$p_i)
    sigma.i <- sqrt(sum((df.p_i$i - mu.i)^2 * df.p_i$p_i))
    F_cm$F_cm_2.5D.corr[iter] <- (1/sigma.i)^2 * sum((df.p_ij$i - mu.i) * (df.p_ij$j - mu.i) * df.p_ij$p_ij)
    
    #Autocorrelation
    F_cm$F_cm_2.5D.auto.corr[iter] <- sum(df.p_ij$i * df.p_ij$j * df.p_ij$p_ij)
    
    #Cluster tendency
    F_cm$F_cm_2.5D.clust.tend[iter] <- sum((df.p_ij$i + df.p_ij$j -2*mu.i)^2 * df.p_ij$p_ij)
    
    #Cluster shade
    F_cm$F_cm_2.5D.clust.shade[iter] <- sum((df.p_ij$i + df.p_ij$j -2*mu.i)^3 * df.p_ij$p_ij)
    
    #Cluster prominence
    F_cm$F_cm_2.5D.clust.prom[iter] <- sum((df.p_ij$i + df.p_ij$j -2*mu.i)^4 * df.p_ij$p_ij)
    
    #First measure of information correlation
    HXY <- - sum(df.p_ij$p_ij * log2(df.p_ij$p_ij))
    HX <- - sum(df.p_i$p_i * log2(df.p_i$p_i))
    HXY1 <- - sum(df.p_ij$p_ij * log2(df.p_ij$p_i * df.p_ij$p_j))
    F_cm$F_cm_2.5D.info.corr.1[iter] <- (HXY - HXY1) / HX
    
    #Second measure of information correlation
    #HXY2 <- - sum(df.p_ij$p_i * df.p_ij$p_j * log2(df.p_ij$p_i * df.p_ij$p_j))
    # if (HXY > HXY2 ) F_cm$F_cm.info.corr.2[iter] <- 0
    # else if (HXY <= HXY2 ) F_cm$F_cm.info.corr.2[iter] <- sqrt(1-exp(-2*(HXY2 - HXY)))
    
    HXY2 <- - sum(df.p_ij$p_i * df.p_ij$p_j * log2(df.p_ij$p_i * df.p_ij$p_j))
    F_cm$F_cm_2.5D.info.corr.2[iter] <- sqrt(1-exp(-2*(abs(HXY2 - HXY))))
  }
  
  return(F_cm)
}

