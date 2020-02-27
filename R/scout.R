#' scout a directory 
#' @export
#' @import stringr progress
scout <- function( ) {
  internalAttributes <- list()
  cacheArea <- list()
  counter <- 0
  
  
  openDICOMFolder <- function( pathToOpen , recursive = TRUE ) {

    # res <- build.preMM( pathToOpen = pathToOpen , recursive = recursive )
    
    # browser()
    load("/projects/moddicom/MV3/res.RData")
    
    MM <- res$MM
    # MM <- cbind(MM , rep("" ,nrow(MM))  )
    # MM <- cbind(MM , rep("" ,nrow(MM))  )
    # MM <- cbind(MM , rep("" ,nrow(MM))  )
    # MM <- cbind(MM , rep("" ,nrow(MM))  )
    # MM <- cbind(MM , rep("" ,nrow(MM))  )
    # MM <- cbind(MM , rep("" ,nrow(MM))  )
    # colnames(MM) <- c(colnames( res$MM ),c("PatientID","StudyDescription","SeriesDescription","BodyPartExamined","ImageOrientationPatient","SeriesDate","CorrectedImage"))
    
    colonne.da.aggiungere <- c("PatientID","StudyDescription","SeriesDescription","BodyPartExamined","ImageOrientationPatient","SeriesDate","ConvolutionKernel","CorrectedImage","Radiopharmaceutical","Units")
    for( i in 1:length(colonne.da.aggiungere)) {MM <- cbind(MM , rep("" ,nrow(MM))  )}
    colnames(MM) <- c(colnames( res$MM ),colonne.da.aggiungere)

    for( riga in rownames( MM ) ) {
      cat("\n", riga )
      extractFurtherInfo <- FALSE
      isa.CT <- FALSE; isa.MRI<- FALSE; isa.PET <- FALSE;
      if( "CTImageStorage" %in% colnames(MM) ) {  if(MM[riga,"CTImageStorage"]>0) {extractFurtherInfo <- TRUE ; isa.CT <- TRUE} }
      if( "MRImageStorage" %in% colnames(MM) ) {  if(MM[riga,"MRImageStorage"]>0) {extractFurtherInfo <- TRUE ; isa.MRI <- TRUE} }
      if( "PositronEmissionTomographyImageStorage" %in% colnames(MM) ) {  if(MM[riga,"PositronEmissionTomographyImageStorage"]>0) {extractFurtherInfo <- TRUE ; isa.PET <- TRUE } }

      if(extractFurtherInfo == TRUE ) {
        PatientID <- getTagOntheFly(fileName = paste( c(riga,"//",list.files(riga)[1]),collapse='')  ,tag = "0010,0020")[1]
        StudyDescription <- getTagOntheFly(fileName = paste( c(riga,"//",list.files(riga)[1]),collapse='')  ,tag = "0008,1030")[1]
        SeriesDescription <- getTagOntheFly(fileName = paste( c(riga,"//",list.files(riga)[1]),collapse='')  ,tag = "0008,103e")[1]
        BodyPartExamined <- getTagOntheFly(fileName = paste( c(riga,"//",list.files(riga)[1]),collapse='')  ,tag = "0018,0015")[1]
        ImageOrientationPatient <- getTagOntheFly(fileName = paste( c(riga,"//",list.files(riga)[1]),collapse='')  ,tag = "0020,0037")[1]
        SeriesDate <- getTagOntheFly(fileName = paste( c(riga,"//",list.files(riga)[1]),collapse='')  ,tag = "0008,0021")[1]
        if( isa.PET ) {
          CorrectedImage <- getTagOntheFly(fileName = paste( c(riga,"//",list.files(riga)[1]),collapse='')  ,tag = "0028,0051")[1]
          Radiopharmaceutical  <- getTagOntheFly(fileName = paste( c(riga,"//",list.files(riga)[1]),collapse='')  ,tag = "0018,0031")[1]
          Units  <- getTagOntheFly(fileName = paste( c(riga,"//",list.files(riga)[1]),collapse='')  ,tag = "0054,1001")[1]
        } else {
          CorrectedImage <- "";  Radiopharmaceutical  <- "";  Units  <- ""
        }
        if( isa.CT ) { 
          ConvolutionKernel  <- getTagOntheFly(fileName = paste( c(riga,"//",list.files(riga)[1]),collapse='')  ,tag = "0018,1210")[1] 
        } else {
          ConvolutionKernel <- ""
        }
        
        
        if( !(length(PatientID)>0)) PatientID<-"NA"
        if( !(length(StudyDescription)>0)) StudyDescription<-""
        if( !(length(SeriesDescription)>0)) SeriesDescription<-""
        if( !(length(BodyPartExamined)>0)) BodyPartExamined<-""
        if( !(length(ImageOrientationPatient)>0)) ImageOrientationPatient<-""
        if( !(length(SeriesDate)>0)) SeriesDate<-""
        if( !(length(CorrectedImage)>0)) CorrectedImage<-""
        if( !(length(Radiopharmaceutical)>0)) Radiopharmaceutical<-""
        if( !(length(Units)>0)) Units<-""
        if( !(length(ConvolutionKernel)>0)) Units<-""
        
        MM[riga,"PatientID"] <- PatientID
        MM[riga,"SeriesDate"] <- SeriesDate
        MM[riga,"StudyDescription"] <- StudyDescription
        MM[riga,"SeriesDescription"] <- SeriesDescription
        MM[riga,"BodyPartExamined"] <- BodyPartExamined
        MM[riga,"ImageOrientationPatient"] <- ImageOrientationPatient
        MM[riga,"CorrectedImage"] <- CorrectedImage
        MM[riga,"Radiopharmaceutical"] <- Radiopharmaceutical
        MM[riga,"Units"] <- Units
        MM[riga,"ConvolutionKernel"] <- ConvolutionKernel
      }
    }
    return(MM)
  }
    
  build.preMM <- function( pathToOpen , recursive = TRUE ) {
    lst.res <- list()
    arr.dirs <- list.dirs(pathToOpen , recursive = T)
    
    if( internalAttributes$verbose == TRUE  )    pb <- progress_bar$new(total = length(arr.dirs))
    
    # arr.dirs <- "/media/localadmin/DATA/images/melanoma//sub-970277/ses-20161122131432/01209-ResultsMMOncologyReading"
    
    for( folder in arr.dirs ) {
      arr.folder <- list.dirs(path = folder)
      
      cat("\n",folder)
      
      # Se sono in una foglia
      if( length(arr.folder) == 1 ) {
        lst.res[[ folder ]] <- list()
        arr.files <- list.files(path = arr.folder,full.names = TRUE)
        comando <- "dcmdump"
        
        for( fileName in arr.files) {
          is.error <- FALSE
          cat("\n\t ",fileName)
          
          # if( fileName=="/projects/moddicom/MV3/img.testing.unit//test.01/axial.01.MRI/GUESS.nii.gz") browser()
          
          options(warn=-1)
          stringa <- fileName
          resul <- suppressWarnings(system2(comando,stringa,stdout = TRUE,stderr = TRUE))
          if(substr(resul[1],1, str_length("E: ")) == "E: " | substr(resul[1],1, str_length("W: ")) == "W: " ) {
            is.error <- TRUE
          }
          options(warn=-0)
          
          if( is.error == FALSE ){
            posizione <- which(unlist(lapply(str_locate_all(string = resul,pattern = "(0008,0016)"),length))>0)
            kindOfValue <- str_trim(str_replace_all(resul[posizione],"\\(0008,0016\\) UI =",""))
            if( length(kindOfValue) > 1 ) kindOfValue <- kindOfValue[1]
            kindOfValue <- str_trim(substr(x = kindOfValue, start = 1, stop = str_locate(string = kindOfValue,pattern = "#")-1))
            
            
            if( !( kindOfValue %in% names(lst.res[[ folder ]])) ) lst.res[[ folder ]][[ kindOfValue ]] <- 0
            lst.res[[ folder ]][[ kindOfValue ]] <- lst.res[[ folder ]][[ kindOfValue ]] + 1          
            
            counter <- counter + 1
          }
        }
        if( internalAttributes$verbose == TRUE ) pb$tick()
      }
    }
    # prendi la lista delle SOPClassUIDs
    arr.SOPClassUIDS <- unique(unlist(lapply( 1:length(lst.res), function(x){  unlist(lapply(  names(lst.res[x]), function( y ) { return(   paste(c("SOPUID.",names(lst.res[x][[y]])),collapse = '')  )} ))       }  )))
    MM <- matrix(0,ncol=length(arr.SOPClassUIDS)+1, nrow=length(names(lst.res)))
    colnames(MM) <- c("folder",arr.SOPClassUIDS)
    rownames(MM) <- names(lst.res)
    MM[,"folder"]<- names(lst.res)
    # Ora costruisci la tabella
    for( riga in names(lst.res)) {
      lapply(colnames(MM)[-c(1)] , function(SOP) {  if(!is.null(lst.res[[riga]][[SOP]])) {MM[riga,SOP] <<-lst.res[[riga]][[SOP]]}   }   )
    }
    return( list("MM"=MM, "arr.SOPClassUID"=arr.SOPClassUIDS ) )
  }
  getTagOntheFly<-function( tag, fileName ) {
    comando <- "dcmdump"

    stringa <- fileName
    options(warn=-1)
    resul <- suppressWarnings(system2(comando,stringa,stdout = TRUE,stderr = TRUE))
    if(substr(resul[1],1, str_length("E: dcmdump")) == "E: dcmdump" ) {
      browser()
    }
    options(warn=0)

    posizione <- which(unlist(lapply(str_locate_all(string = resul,pattern = paste(c("(",tag,")"),collapse = '')   ),length))>0)
    kindOfValue <- str_trim(str_replace_all(resul[posizione],paste( c("\\(",tag,"\\)") ,collapse = ''),""))
    kindOfValue <- str_trim(substr(x = kindOfValue, start = 1, stop = str_locate(string = kindOfValue,pattern = "#")-1))
    return(kindOfValue);
  }

  #=================================================================================
  # Constructor
  #=================================================================================
  constructor<-function( ) {
    
    # Attributes - set by user
    internalAttributes$rootDir <<- ""
    internalAttributes$verbose <<- FALSE

    # Internal Structures and objs
    logObj <<- logHandler()                                   # log/error handler Object
    cacheArea <<- list(  )
  }
  constructor( )
  return( list(
    "openDICOMFolder"=openDICOMFolder
  ))
  
}