#' service class
#'
#' @description  It handles messages from script to a chosen output (screen, file, etc.)
#' @useDynLib moddicomV3
#' @export
#' @import stringr XML
services<-function() {
  # ------------------------------------------------
  # rotateMatrix
  # ------------------------------------------------
  rotateMatrix<-function( m , rotations = 1) {
    if(rotations == 1 ) m<-t(m[nrow(m):1,])
    if(rotations == 2 ) m<-m[nrow(m):1,ncol(m):1]
    if(rotations == 3 ) m<-t(m)[ncol(m):1,]
    return( m )
  }
  # ------------------------------------------------
  # getPointPlaneDistance
  # ------------------------------------------------
  getPointPlaneDistance<-function(Punto,Piano) {
    return(abs(Piano[1]*Punto[1]+Piano[2]*Punto[2]+Piano[3]*Punto[3]+Piano[4])/sqrt(Piano[1]^2+Piano[2]^2+Piano[3]^2))
  }
  # ------------------------------------------------
  # get3DPosFromNxNy
  # ------------------------------------------------
  get3DPosFromNxNy<-function(Nx,Ny,oM) {
    return(xy<-t(oM%*%c(Nx,Ny,0,1)))
  }
  # ------------------------------------------------
  # getPlaneEquationBetween3Points
  # ------------------------------------------------
  getPlaneEquationBetween3Points<-function(Pa,Pb,Pc) {
    ac<-(Pb[2]-Pa[2])*(Pc[3]-Pa[3])-(Pc[2]-Pa[2])*(Pb[3]-Pa[3])
    bc<-(Pb[3]-Pa[3])*(Pc[1]-Pa[1])-(Pc[3]-Pa[3])*(Pb[1]-Pa[1])
    cc<-(Pb[1]-Pa[1])*(Pc[2]-Pa[2])-(Pc[1]-Pa[1])*(Pb[2]-Pa[2])
    dc<--(ac*Pa[1]+bc*Pa[2]+cc*Pa[3])
    return(c(ac,bc,cc,dc))
  }

  # ===============================================================
  # getXMLStructureFromDICOMFile
  # ===============================================================
  getXMLStructureFromDICOMFile<-function(fileName, folderCleanUp = TRUE) {
    # aggiungi estensione xml
    fileNameXML<-paste(fileName,".xml")
    fileNameXML<-str_replace_all(string = fileNameXML , pattern = " .xml",replacement = ".xml")
    # salva solo path senza nome oggetto DICOM
    pathToStore<-substr(fileName,1,tail(which(strsplit(fileName, '')[[1]]=='/'),1)-1)
    # se file con estensione xml gia' esiste nella cartella non fare nulla altrimenti lo aggiunge
    if(!file.exists( fileNameXML ) | folderCleanUp==TRUE) {
      stringa1<-"dcm2xml";
      p.fileNameXML <- paste( c("'",fileNameXML,"'") ,collapse = '')
      p.fileName <- paste( c("'",fileName,"'") ,collapse = '')
      stringa2<-paste(" +M  ",p.fileName,p.fileNameXML,collapse='')
      options(warn=-1)
      system2(stringa1,stringa2,stdout=NULL)
      options(warn=0)
    }
    # Load the XML file: restituisce il file xml nella variabile doc
    doc = xmlInternalTreeParse(fileNameXML)
    return(doc);
  }

  # ===============================================================
  # getDICOMTag
  # ===============================================================
  getDICOMTag<-function(tag=tag, fileName="", folderCleanUp = TRUE) {
    obj.S<-services();

    # exemption: you want an Image!
    # if(tag == "7fe0,0010") return( getImageFromRAW(fileName)  );
    if(tag == "7fe0,0010") stop("Not available Yet! (for tag 7fe0,0010)")

    # build the XML file and get the XML structure
    doc<-getXMLStructureFromDICOMFile(fileName = fileName, folderCleanUp = folderCleanUp)

    # build the QUERY
    stringaQuery<-''
    if(tag=="0018,0031" | tag=="0018,1072" | tag=="0018,1074" | tag=="0018,1075") {
      stringaQuery<-paste(c('/file-format/data-set/sequence[@tag="0054,0016"]/item/element[@tag="',tag,'"]'),collapse='');
    }
    if(stringaQuery=='') stringaQuery<-paste(c('/file-format/data-set/element[@tag="',tag,'"]'),collapse='');

    # execute the QUERY
    valore<-xpathApply(doc,stringaQuery,xmlValue);
    if(length(valore)==2) logObj$handle( "error" , "a tag in DICOM file seems to be duplicated"  );
    if(length(valore)==0) return(NA);

    valore<-valore[[1]]
    return(valore);
  }

  # ===============================================================
  # splittaTAG
  # ===============================================================
  splittaTAG<-function(stringa) {
    return( as.numeric(strsplit(stringa,split = "\\\\")[[1]])   )
  }

  # ===============================================================
  # new.trilinearInterpolator
  # ===============================================================
  new.trilinearInterpolator<-function( voxelCube , pixelSpacing.new  ,pixelSpacing.old  ) {
    Nx.old<-dim(voxelCube)[1];	Ny.old<-dim(voxelCube)[2];	Nz.old<-dim(voxelCube)[3]
    xDim.old<-pixelSpacing.old[1];	yDim.old<-pixelSpacing.old[2];	zDim.old<-pixelSpacing.old[3]
    xDim.new<-pixelSpacing.new[1];	yDim.new<-pixelSpacing.new[2];	zDim.new<-pixelSpacing.new[3]

    fattoreDiScalaX<-pixelSpacing.old[1]/pixelSpacing.new[1];
    fattoreDiScalaY<-pixelSpacing.old[2]/pixelSpacing.new[2];
    fattoreDiScalaZ<-pixelSpacing.old[3]/pixelSpacing.new[3];

    Nx.new<-ceiling(Nx.old * fattoreDiScalaX)
    Ny.new<-ceiling(Ny.old * fattoreDiScalaY)
    Nz.new<-ceiling(Nz.old * fattoreDiScalaZ)

    result<-array(rep( 0 , Nx.new * Ny.new * Nz.new ))

    res<-.C("newnewtrilinearInterpolator",
            as.integer(Nx.old),as.integer(Ny.old),as.integer(Nz.old),
            as.integer(Nx.new),as.integer(Ny.new),as.integer(Nz.new),
            as.double(pixelSpacing.old[1]),as.double(pixelSpacing.old[2]),as.double(pixelSpacing.old[3]),
            as.double(pixelSpacing.new[1]),as.double(pixelSpacing.new[2]),as.double(pixelSpacing.new[3]),
            as.double(voxelCube),as.double(result) );
    result<-array( res[[14]] , dim=c(Nx.new,Ny.new,Nz.new) )
    return( result )
  }
  # ========================================================================================
  # cropCube: crop a voxel cube in order to limit its dimension to the needs
  # ========================================================================================
  cropCube<-function( bigCube ) {
    matPos<-which(bigCube!=0,arr.ind = T)
    min.x<-min(matPos[,1]);     max.x<-max(matPos[,1])
    min.y<-min(matPos[,2]);     max.y<-max(matPos[,2])
    min.z<-min(matPos[,3]);     max.z<-max(matPos[,3])
    newCube<-bigCube[ min.x:max.x, min.y:max.y , min.z:max.z]
    location<-list( "min.x"=min.x, "max.x"=max.x, "min.y"=min.y, "max.y"=max.y, "min.z"=min.z, "max.z"=max.z  )
    return( list ( "voxelCube"=newCube, "location"=location) )
  }

  return( list(
    "get3DPosFromNxNy"=get3DPosFromNxNy,
    "getPlaneEquationBetween3Points"=getPlaneEquationBetween3Points,
    "getPointPlaneDistance"=getPointPlaneDistance,
    "getXMLStructureFromDICOMFile" = getXMLStructureFromDICOMFile,
    "getDICOMTag"=getDICOMTag,
    "splittaTAG"=splittaTAG,
    "new.trilinearInterpolator"=new.trilinearInterpolator,
    "rotateMatrix"=rotateMatrix,
    "cropCube"=cropCube
  ))
}
