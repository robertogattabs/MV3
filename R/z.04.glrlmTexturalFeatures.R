#' GLRM  features
#'
#' @description  Extract the GLRM features
#' @param imgObj a 3D matrix
#' @export
#' @import radiomics data.table
glrlmTexturalFeatures <- function(imgObj, n_grey){
  
  # compute number of non-NA voxels
  
  #nVoxel <- dim(imgObj)[1]*dim(imgObj)[2]*dim(imgObj)[3] - sum(is.na(imgObj))
  Nv <- length(imgObj[which(!is.na(imgObj))])
  
  ### compute Gray Levels Cooccurrence Matrices
  
  if(length(dim(imgObj))<3) { 
    imgObj <- array(imgObj, dim=c(dim(imgObj)[1],dim(imgObj)[2],1))
  } 
  
  R_list <- list()
  Nv <- numeric()
  #Compute grey level cooccurrence matrices for 4 different directions within each slice
  for (i in 1:dim(imgObj)[3]){
    if(length(imgObj[,,i])*.005 >= sum(!is.na(imgObj[,,i]))) next
    #if(length(table(unique(imgObj[,,i])))<n_grey) n_grey <- length(table(unique(imgObj[,,i])))-1
    if (min(imgObj[,,i],na.rm = T) < 0) {imgObj[,,i] <- imgObj[,,i] + abs(min(imgObj[,,i],na.rm = T))}
    R_list[[(i-1)*4+1]] <- as.matrix(glrlm(imgObj[,,i], angle = 0, verbose=F,truncate = F, n_grey = n_grey))
    R_list[[(i-1)*4+2]] <- as.matrix(glrlm(imgObj[,,i], angle = 45, verbose=F,truncate = F, n_grey = n_grey))
    R_list[[(i-1)*4+3]] <- as.matrix(glrlm(imgObj[,,i], angle = 90, verbose=F,truncate = F, n_grey =n_grey))
    R_list[[(i-1)*4+4]] <- as.matrix(glrlm(imgObj[,,i], angle = 135, verbose=F,truncate = F, n_grey = n_grey))
    Nv[(i-1)*4+1] <- length(imgObj[,,i][which(!is.na(imgObj[,,i]))])
    Nv[(i-1)*4+2] <- length(imgObj[,,i][which(!is.na(imgObj[,,i]))])
    Nv[(i-1)*4+3] <- length(imgObj[,,i][which(!is.na(imgObj[,,i]))])
    Nv[(i-1)*4+4] <- length(imgObj[,,i][which(!is.na(imgObj[,,i]))])
  }
  
  #elimino gli elementi NULL della lista
  if(length(R_list)!=0 & length(Nv)!=0){
    if(length(which(sapply(R_list, is.null)))!=0){
      R_list = R_list[-which(sapply(R_list, is.null))]
    }
    if(length(which(sapply(Nv, is.null)))!=0){
      Nv = Nv[-which(sapply(Nv, is.null))]
    }
  }
  #Initialise data table for storing GLCM features; I have added a few
  featNames <- c("F_rlm.sre","F_rlm.lre","F_rlm.lgre","F_rlm.hgre","F_rlm.srlge",
                 "F_rlm.srhge","F_rlm.lrlge","F_rlm.lrhge","F_rlm.glnu","F_rlm.glnu.norm","F_rlm.rlnu","F_rlm.rlnu.norm", "F_rlm.r.perc",
                 "F_rlm.gl.var","F_rlm.rl.var","F_rlm.rl.entr")
  F_rlm <- data.table(matrix(NA, nrow=length(R_list), ncol=length(featNames)))
  colnames(F_rlm) <- featNames
  
  #Iterate over grey level cooccurrence matrices
  #The idea is that basically every GLCM is the same, i.e. we can just perform the same operations on every glcm.
  for (iter in seq_len(length(R_list))){
    
    # if(iter==1 | iter==5 | iter==9 | iter==13) {    Nv <- length(imgObj[,,1][which(!is.na(imgObj[,,1]))]) }
    # else if(iter==2 | iter==6 | iter==10 | iter==14) {    Nv <- length(imgObj[,,2][which(!is.na(imgObj[,,2]))]) }
    # else if(iter==3 | iter==7 | iter==11 | iter==15) {    Nv <- length(imgObj[,,3][which(!is.na(imgObj[,,3]))]) }
    # else if(iter==4 | iter==8 | iter==12 | iter==16) {    Nv <- length(imgObj[,,4][which(!is.na(imgObj[,,4]))]) }
    
    Ns <- sum(R_list[[iter]],na.rm=T)
    
    #Convert matrix to data table
    df.R    <- data.table(R_list[[iter]])
    
    #Add row grey level intensity
    df.R$i <- as.numeric(row.names(R_list[[iter]]))
    
    #Convert from wide to long table. This is the preferred format for data tables and data frames
    df.R   <- melt(df.R, id.vars="i", variable.name="j", value.name="n", variable.factor=FALSE)
    
    #Convert j from string to numeric
    df.R$j <- as.numeric(df.R$j)
    
    #Remove combinations with 0 counts
    df.R   <- df.R[n>0,]
    
    df.R <- df.R[order(rank(i))]
    
    #Convert Grey level coccurrence matrix to joint probability
    #df.r_ij <- df.R[,.(r_ij=n/sum(df.R$n)), by=.(i,j)] #joint probability
    
    df.r_i  <- df.R[,.(r_i=sum(n)), by=i]        #marginal probability over columns
    df.r_j  <- df.R[,.(r_j=sum(n)), by=j]        #marginal probability over rows
    
    #Diagonal probabilities (p(i-j))
    #First, we create a new column k which contains the absolute value of i-j.
    #Second, we sum the joint probability where k is the same.
    #This can written as one line by chaining the operations.
    # df.p_imj <- copy(df.p_ij)
    # df.p_imj <- df.p_imj[,"k":=abs(i-j)][,.(p_imj=sum(p_ij)), by=k]
    
    #Cross-diagonal probabilities (p(i+j))
    #Again, we first create a new column k which contains i+j
    #Second, we sum the probability where k is the same.
    #This is written in one line by chaining the operations.
    # df.p_ipj <- copy(df.p_ij)
    # df.p_ipj <- df.p_ipj[,"k":=i+j][,.(p_ipj=sum(p_ij)), by=k]
    
    #Merger of df.p_ij, df.p_i and df.p_j
    df.R <- merge(x=df.R, y=df.r_i, by="i")
    df.R <- merge(x=df.R, y=df.r_j, by="j")
    
    #Thus we have five probability matrices
    #Joint probability:          df.p_ij with probabilities p_ij, p_i and p_j, and indices i, j
    #Marginal probability:       df.p_i with probability p_i, and index i
    #Marginal probability:       df.p_j with probability p_j, and index j
    #Diagonal probability:       df.p_imj with probability p_imj and index k
    #Cross-diagonal probability: df.p_ipj with probability p_ipj and index k
    
    #Calculate features
    #Short runs emphasis
    F_rlm$F_rlm.sre[iter]         <- (1/Ns) * sum(df.r_j$r_j/(df.r_j$j^2))
    
    #Long runs emphasis
    F_rlm$F_rlm.lre[iter]         <- (1/Ns) * sum(df.r_j$r_j*(df.r_j$j^2))
    
    #Low grey level run emphasis
    F_rlm$F_rlm.lgre[iter]         <- (1/Ns) * sum(df.r_i$r_i/(df.r_i$i^2))
    
    #High grey level run emphasis
    F_rlm$F_rlm.hgre[iter]         <- (1/Ns) * sum(df.r_i$r_i*(df.r_i$i^2))
    
    #Short run low grey level emphasis
    F_rlm$F_rlm.srlge[iter] <- (1/Ns) * sum(df.R$n/((df.R$i^2)*(df.R$j^2)))
    
    #Short run high grey level emphasis
    F_rlm$F_rlm.srhge[iter] <- (1/Ns) * sum((df.R$n)*(df.R$i^2)/(df.R$j^2))
    
    #Long run low grey level emphasis
    F_rlm$F_rlm.lrlge[iter] <- (1/Ns) * sum((df.R$n)*(df.R$j^2)/(df.R$i^2))
    
    #Long run high grey level emphasis
    F_rlm$F_rlm.lrhge[iter] <- (1/Ns) * sum((df.R$n)*(df.R$j^2)*(df.R$i^2))
    
    #Grey level non-uniformity
    F_rlm$F_rlm.glnu[iter] <- (1/Ns) * sum(df.r_i$r_i^2)
    
    #Grey level non-uniformity normalized
    F_rlm$F_rlm.glnu.norm[iter] <- (1/Ns^2) * sum(df.r_i$r_i^2)
    
    #Run length non-uniformity
    F_rlm$F_rlm.rlnu[iter] <- (1/Ns) * sum(df.r_j$r_j^2)
    
    #Run length non-uniformity normalized
    F_rlm$F_rlm.rlnu.norm[iter] <- (1/Ns^2) * sum(df.r_j$r_j^2)
    
    #Run percentage
    
    F_rlm$F_rlm.r.perc[iter] <- Ns/Nv[iter]
    
    #Grey level variance
    p_ij <- df.R$n/sum(df.R$n)
    mu_i <- sum(df.R$i * p_ij)
    F_rlm$F_rlm.gl.var[iter] <- sum((df.R$i - mu_i)^2 * p_ij)
    
    #Run length variance
    mu_j <- sum(df.R$j * p_ij)
    F_rlm$F_rlm.rl.var[iter] <- sum((df.R$j - mu_j)^2 * p_ij)
    
    #Run entropy
    F_rlm$F_rlm.rl.entr[iter] <- - sum(p_ij * log2(p_ij))
  }
  
  return(F_rlm)
}


# old.glrlmTexturalFeatures <- function(imgObj, n_grey){
#   
#   # compute number of non-NA voxels
#   
#   #nVoxel <- dim(imgObj)[1]*dim(imgObj)[2]*dim(imgObj)[3] - sum(is.na(imgObj))
#   Nv <- length(imgObj[which(!is.na(imgObj))])
#   
#   ### compute Gray Levels Cooccurrence Matrices
#   
#   R_list <- list()
#   Nv <- numeric()
#   #Compute grey level cooccurrence matrices for 4 different directions within each slice
#   for (i in 1:dim(imgObj)[3]){
#     if(length(imgObj[,,i])*.005 >= sum(!is.na(imgObj[,,i]))) next
#     #if(length(table(unique(imgObj[,,i])))<n_grey) n_grey <- length(table(unique(imgObj[,,i])))-1
#     if (min(imgObj[,,i],na.rm = T) < 0) {imgObj[,,i] <- imgObj[,,i] + abs(min(imgObj[,,i],na.rm = T))}
#     R_list[[(i-1)*4+1]] <- as.matrix(glrlm(imgObj[,,i], angle = 0, verbose=F,truncate = F, n_grey = n_grey))
#     R_list[[(i-1)*4+2]] <- as.matrix(glrlm(imgObj[,,i], angle = 45, verbose=F,truncate = F, n_grey = n_grey))
#     R_list[[(i-1)*4+3]] <- as.matrix(glrlm(imgObj[,,i], angle = 90, verbose=F,truncate = F, n_grey =n_grey))
#     R_list[[(i-1)*4+4]] <- as.matrix(glrlm(imgObj[,,i], angle = 135, verbose=F,truncate = F, n_grey = n_grey))
#     Nv[(i-1)*4+1] <- length(imgObj[,,i][which(!is.na(imgObj[,,i]))])
#     Nv[(i-1)*4+2] <- length(imgObj[,,i][which(!is.na(imgObj[,,i]))])
#     Nv[(i-1)*4+3] <- length(imgObj[,,i][which(!is.na(imgObj[,,i]))])
#     Nv[(i-1)*4+4] <- length(imgObj[,,i][which(!is.na(imgObj[,,i]))])
#   }
#   
#   #elimino gli elementi NULL della lista
#   if(length(R_list)!=0 & length(Nv)!=0){
#     if(length(which(sapply(R_list, is.null)))!=0){
#       R_list = R_list[-which(sapply(R_list, is.null))]
#     }
#     if(length(which(sapply(Nv, is.null)))!=0){
#       Nv = Nv[-which(sapply(Nv, is.null))]
#     }
#   }
#   #Initialise data table for storing GLCM features; I have added a few
#   featNames <- c("F_rlm.sre","F_rlm.lre","F_rlm.lgre","F_rlm.hgre","F_rlm.srlge",
#                  "F_rlm.srhge","F_rlm.lrlge","F_rlm.lrhge","F_rlm.glnu","F_rlm.glnu.norm","F_rlm.rlnu","F_rlm.rlnu.norm", "F_rlm.r.perc",
#                  "F_rlm.gl.var","F_rlm.rl.var","F_rlm.rl.entr")
#   F_rlm <- data.table(matrix(NA, nrow=length(R_list), ncol=length(featNames)))
#   colnames(F_rlm) <- featNames
#   
#   #Iterate over grey level cooccurrence matrices
#   #The idea is that basically every GLCM is the same, i.e. we can just perform the same operations on every glcm.
#   for (iter in seq_len(length(R_list))){
#     
#     # if(iter==1 | iter==5 | iter==9 | iter==13) {    Nv <- length(imgObj[,,1][which(!is.na(imgObj[,,1]))]) }
#     # else if(iter==2 | iter==6 | iter==10 | iter==14) {    Nv <- length(imgObj[,,2][which(!is.na(imgObj[,,2]))]) }
#     # else if(iter==3 | iter==7 | iter==11 | iter==15) {    Nv <- length(imgObj[,,3][which(!is.na(imgObj[,,3]))]) }
#     # else if(iter==4 | iter==8 | iter==12 | iter==16) {    Nv <- length(imgObj[,,4][which(!is.na(imgObj[,,4]))]) }
#     
#     Ns <- sum(R_list[[iter]],na.rm=T)
#     
#     #Convert matrix to data table
#     df.R    <- data.table(R_list[[iter]])
#     
#     #Add row grey level intensity
#     df.R$i <- as.numeric(row.names(R_list[[iter]]))
#     
#     #Convert from wide to long table. This is the preferred format for data tables and data frames
#     df.R   <- melt(df.R, id.vars="i", variable.name="j", value.name="n", variable.factor=FALSE)
#     
#     #Convert j from string to numeric
#     df.R$j <- as.numeric(df.R$j)
#     
#     #Remove combinations with 0 counts
#     df.R   <- df.R[n>0,]
#     
#     df.R <- df.R[order(rank(i))]
#     
#     #Convert Grey level coccurrence matrix to joint probability
#     #df.r_ij <- df.R[,.(r_ij=n/sum(df.R$n)), by=.(i,j)] #joint probability
#     
#     df.r_i  <- df.R[,.(r_i=sum(n)), by=i]        #marginal probability over columns
#     df.r_j  <- df.R[,.(r_j=sum(n)), by=j]        #marginal probability over rows
#     
#     #Diagonal probabilities (p(i-j))
#     #First, we create a new column k which contains the absolute value of i-j.
#     #Second, we sum the joint probability where k is the same.
#     #This can written as one line by chaining the operations.
#     # df.p_imj <- copy(df.p_ij)
#     # df.p_imj <- df.p_imj[,"k":=abs(i-j)][,.(p_imj=sum(p_ij)), by=k]
#     
#     #Cross-diagonal probabilities (p(i+j))
#     #Again, we first create a new column k which contains i+j
#     #Second, we sum the probability where k is the same.
#     #This is written in one line by chaining the operations.
#     # df.p_ipj <- copy(df.p_ij)
#     # df.p_ipj <- df.p_ipj[,"k":=i+j][,.(p_ipj=sum(p_ij)), by=k]
#     
#     #Merger of df.p_ij, df.p_i and df.p_j
#     df.R <- merge(x=df.R, y=df.r_i, by="i")
#     df.R <- merge(x=df.R, y=df.r_j, by="j")
#     
#     #Thus we have five probability matrices
#     #Joint probability:          df.p_ij with probabilities p_ij, p_i and p_j, and indices i, j
#     #Marginal probability:       df.p_i with probability p_i, and index i
#     #Marginal probability:       df.p_j with probability p_j, and index j
#     #Diagonal probability:       df.p_imj with probability p_imj and index k
#     #Cross-diagonal probability: df.p_ipj with probability p_ipj and index k
#     
#     #Calculate features
#     #Short runs emphasis
#     F_rlm$F_rlm.sre[iter]         <- (1/Ns) * sum(df.r_j$r_j/(df.r_j$j^2))
#     
#     #Long runs emphasis
#     F_rlm$F_rlm.lre[iter]         <- (1/Ns) * sum(df.r_j$r_j*(df.r_j$j^2))
#     
#     #Low grey level run emphasis
#     F_rlm$F_rlm.lgre[iter]         <- (1/Ns) * sum(df.r_i$r_i/(df.r_i$i^2))
#     
#     #High grey level run emphasis
#     F_rlm$F_rlm.hgre[iter]         <- (1/Ns) * sum(df.r_i$r_i*(df.r_i$i^2))
#     
#     #Short run low grey level emphasis
#     F_rlm$F_rlm.srlge[iter] <- (1/Ns) * sum(df.R$n/((df.R$i^2)*(df.R$j^2)))
#     
#     #Short run high grey level emphasis
#     F_rlm$F_rlm.srhge[iter] <- (1/Ns) * sum((df.R$n)*(df.R$i^2)/(df.R$j^2))
#     
#     #Long run low grey level emphasis
#     F_rlm$F_rlm.lrlge[iter] <- (1/Ns) * sum((df.R$n)*(df.R$j^2)/(df.R$i^2))
#     
#     #Long run high grey level emphasis
#     F_rlm$F_rlm.lrhge[iter] <- (1/Ns) * sum((df.R$n)*(df.R$j^2)*(df.R$i^2))
#     
#     #Grey level non-uniformity
#     F_rlm$F_rlm.glnu[iter] <- (1/Ns) * sum(df.r_i$r_i^2)
#     
#     #Grey level non-uniformity normalized
#     F_rlm$F_rlm.glnu.norm[iter] <- (1/Ns^2) * sum(df.r_i$r_i^2)
#     
#     #Run length non-uniformity
#     F_rlm$F_rlm.rlnu[iter] <- (1/Ns) * sum(df.r_j$r_j^2)
#     
#     #Run length non-uniformity normalized
#     F_rlm$F_rlm.rlnu.norm[iter] <- (1/Ns^2) * sum(df.r_j$r_j^2)
#     
#     #Run percentage
#     
#     F_rlm$F_rlm.r.perc[iter] <- Ns/Nv[iter]
#     
#     #Grey level variance
#     p_ij <- df.R$n/sum(df.R$n)
#     mu_i <- sum(df.R$i * p_ij)
#     F_rlm$F_rlm.gl.var[iter] <- sum((df.R$i - mu_i)^2 * p_ij)
#     
#     #Run length variance
#     mu_j <- sum(df.R$j * p_ij)
#     F_rlm$F_rlm.rl.var[iter] <- sum((df.R$j - mu_j)^2 * p_ij)
#     
#     #Run entropy
#     F_rlm$F_rlm.rl.entr[iter] <- - sum(p_ij * log2(p_ij))
#   }
#   
#   return(F_rlm)
# }