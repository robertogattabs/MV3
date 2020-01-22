#' class for loading and presenting DICOM data
#'
#' @description  Instantiate an object of the class \code{geoLet}.This represents just the classname,
#'               methods are exposed with the technique of 'closure'.
#'               In order to see manuals for the single mathods, consider the vignette or use the
#'               available for the following wrapping functions:
#'               \itemize{
#'               \item \code{GLT.openDICOMFolder( );} : to load a DICOM series into an geoLet object
#'               \item \code{GLT.getImageVoxelCube( );} : to get the ImageVoxelCube stored into a geoLet object
#'               \item \code{GLT.getPixelSpacing( );} : to get the pixelSpacing (x,y,z) of the main ImageVoxelCube stored into a geoLet object
#'               \item \code{GLT.getROIList( );} : to get the list of the ROI defined in a geoLet object
#'               \item \code{GLT.getTag( );} : to get a single DICOM-tag of a DICOM file loaded into a geoLet object
#'               \item \code{GLT.getROIVoxels( );} : to get the IMAGE Voxels geometrically located into a ROI, for a given geoLet object
#'               \item \code{GLT.extractDoseVoxels( );} : to get the DOSE Voxels geometrically located into a ROI, for a given geoLet object
#'               \item \code{GLT.calculateDVH( );} : to get the DVH calculated from a geoLet object
#'               }
#'               The original methods for the class geoLet can also be invocked using the same name without the previx 'GTL.', i.e.:
#'               misc3d rgl Rvcg oce rmarkdown moments
#' @export
#' @useDynLib moddicomV3
#' @import stringr XML oro.nifti
geoLet<-function() {

  # global variables
  internalAttributes<-list()                             # Attributes
  logObj<-logHandler()                                   # log/error handler Object
  dataStorage<-list()                                    # memory data structure
  SOPClassUIDList<-c()
  global_tableROIPointList<-c()

  #=================================================================================
  # openDICOMFolder
  # Loads a Folder containing one or more DICOM Studies
  #=================================================================================
  # Open a folder and load the content
  openDICOMFolder<-function(pathToOpen) {

    if(!dir.exists(pathToOpen)) logObj$handle( "error" , "The indicate Path does not exist"  );

    # ----------------------------------------------
    # get the dcm file type
    # ----------------------------------------------
    if( internalAttributes$verbose == TRUE ) cat("\n Dir scouting: ")
    SOPClassUIDList<<-getFolderContent( pathToOpen );

    # ----------------------------------------------
    # Load CT/RMN Scans
    # ----------------------------------------------
    if( internalAttributes$verbose == TRUE ) cat("\n Image Loading: ")
    loadCTRMNRDScans( );

    # ----------------------------------------------
    # Carica l'RTStruct, se presente
    # ----------------------------------------------
    if( internalAttributes$verbose == TRUE ) cat("\n RTStruct Loading: ")
    loadRTStructFiles()
    if( internalAttributes$verbose == TRUE ) cat("Done")

    # ----------------------------------------------
    # Carica i nifti, se presenti
    # ----------------------------------------------
    if( internalAttributes$verbose == TRUE ) cat("\n nifti files Loading: ")
    loadNIFTIFileDescription()
    if( internalAttributes$verbose == TRUE ) cat("Done")
  }

  #=================================================================================
  # loadNIFTIFiles
  # Loads the nifti files in the folder
  #=================================================================================
  loadNIFTIFileDescription<-function( ) {
    for( riga in 1:nrow(SOPClassUIDList) ) {
      if( SOPClassUIDList[ riga , "kind"] == "nifti" ) {
        fileNameWithPath <- SOPClassUIDList[ riga , "fileName"]
        lastBS <- rev(unlist(str_locate_all(fileNameWithPath,"/")))[1]
        ROIName <- str_trim(str_replace_all(str_trim(str_sub(fileNameWithPath, lastBS+1)),".nii.gz",""))
        ROIName <- paste(c(ROIName,".nii"), collapse = '')

        # fileNameWithPath<-SOPClassUIDList[ riga , "fileName"]
        # aaa <- readNIfTI(fname = fileNameWithPath)
        # tmpVC <- slot(aaa,".Data")
        # if( slot(aaa,"reoriented") == FALSE ) logObj$sendLog(  "In the NIFTI file, 'reoriented' is set to FALSE" ,"ERR" );
        # if( slot(aaa,"scl_slope") != 1 ) logObj$sendLog(  "In the NIFTI file, 'slope' is no 1" ,"ERR" );
        # if( slot(aaa,"scl_inter") != 0 ) logObj$sendLog(  "In the NIFTI file, 'intercept' is no 1" ,"ERR" );
        # dim.x <- slot(aaa,"dim_")[2]
        # dim.y <- slot(aaa,"dim_")[3]
        # dim.z <- slot(aaa,"dim_")[4]
        dataStorage$structures[[ROIName]] <<- list()
        if( !("structures" %in% names(dataStorage$info)) ) dataStorage$info$structures <<- list()
        dataStorage$info$structures[[ROIName]]$type <<- "NIFTI"
        dataStorage$info$structures[[ROIName]]$fileName <<- fileNameWithPath
        dataStorage$info$structures[[ROIName]]$loaded <<- FALSE
      }
    }
  }
  #=================================================================================
  # loadRTStructFiles
  # Loads a DICOM RT Struct (one x folder)
  #=================================================================================
  loadRTStructFiles<-function( ) {

    imageSerie<-list();   listaPuntiROI<-list()
    explicitRTStructFileName <- internalAttributes$explicitRTStructFileName

    TMP<-list()
    if( is.na(explicitRTStructFileName)) {
      righe.RTStruct <- which(SOPClassUIDList[,"kind"] == "RTStructureSetStorage")
      for( riga in righe.RTStruct ) {
        fileName <- SOPClassUIDList[riga ,"fileName"]
        SOPInstanceUID <- SOPClassUIDList[riga ,"SOPInstanceUID"]
        TMP[[ SOPInstanceUID ]]<-getStructuresFromXML( fileName );
      }
    } else {
      TMP[[ explicitRTStructFileName ]]<-getStructuresFromXML( explicitRTStructFileName );
    }

    # now let me use some more easy to handle variable names
    matrice2<-c(); matrice3<-c(); FORUID.m<-NA;
    for(i in names(TMP)) {
      matrice2<-cbind(matrice2,TMP[[i]]$IDROINameAssociation)
      matrice3<-rbind(matrice3,TMP[[i]]$tableROIPointList)
      # Aggiungi le informazioni relative al FrameOfReferenceUID delle ROI caricate
      # ed il ReferencedROINumber!!!!! (xè non è il numero della ROI di moddicom, possono essere diversi)
      if( !("structures" %in% names(dataStorage$info)) ) dataStorage$info$structures<<-list();
      for( nomeROI in TMP[[i]]$IDROINameAssociation[2,] ) {
        dataStorage$info$structures[[nomeROI]]<<-list();
        dataStorage$info$structures[[nomeROI]]$FrameOfReferenceUID<<-TMP[[i]]$FORUID.m
        dataStorage$info$structures[[nomeROI]]$SeriesInstanceUID<<-TMP[[i]]$RTStructSeriesInstanceUID
        dataStorage$info$structures[[nomeROI]]$type<<-"DICOMRTStruct"
        dataStorage$info$structures[[nomeROI]]$associatedSlices<<-c()
      }
    }
    listaROI<-list()

    # for each ROI
    for(i in matrice2[2,]) {

      # get the points
      subMatrix<-matrice3[which(matrice3[,2]==i,arr.ind = TRUE),]
      # if some points exist
      quantiElementiTrovati <- -1

      if(is.list(subMatrix) & !is.array(subMatrix)) quantiElementiTrovati<-1
      if(length(subMatrix)==4 & !is.array(subMatrix) & is.matrix(subMatrix)==FALSE) subMatrix <- t(subMatrix)
      if(is.matrix(subMatrix) & is.array(subMatrix)) quantiElementiTrovati<-dim(subMatrix)[1]

      if(quantiElementiTrovati==-1) {
        logObj$sendLog( "Unexpected error in loading slices. No slices found.", "ERR"  );
      }

      if( quantiElementiTrovati >0 ) {
        listaROI[[i]]<-list()
        # add properly the points to the 'listaROI' structure
        for(contatore in seq(1,quantiElementiTrovati) ) {
          if( quantiElementiTrovati == 1) {
            ROIPointStringList<-subMatrix[[3]]
            SOPInstance<-subMatrix[[4]]
          }
          else {
            ROIPointStringList<-subMatrix[contatore,3][[1]]
            SOPInstance<-subMatrix[contatore,4][[1]]
          }
          listaCoords<-strsplit(ROIPointStringList,"\\\\");
          listaCoords<-as.numeric(listaCoords[[1]])
          # if a ROI already exists for the slice, well append it to the list
          if( !( SOPInstance  %in% names(listaROI[[i]])  ) )  listaROI[[i]][[   SOPInstance  ]]<-list()
          listaROI[[i]][[   SOPInstance  ]][[ length(listaROI[[i]][[   SOPInstance  ]])+1  ]]<-matrix(listaCoords,ncol=3,byrow=T)
          # Add the first one as last (close the loop)
          listaROI[[i]][[   SOPInstance  ]][[ length(listaROI[[i]][[   SOPInstance  ]]) ]]<-rbind(listaROI[[i]][[   SOPInstance  ]][[length(listaROI[[i]][[   SOPInstance  ]])]],listaROI[[i]][[   SOPInstance  ]][[length(listaROI[[i]][[   SOPInstance  ]])]][1,])
        }
      } else {
        listaROI[[i]]<-NA
      }
    }

    for( tmpSOPIUID in names(TMP)) {
      tableROIPointList <- TMP[[tmpSOPIUID]]$tableROIPointList
      ROINamesDaSistemare <- unique(tableROIPointList[,"ROIName"])
      for(ROIName in ROINamesDaSistemare ) {
        SOTTOMATRICE <- tableROIPointList[ which(tableROIPointList[,"ROIName"]==ROIName), ]
        dataStorage$info$structures[[ROIName]]$associatedSlices <<- SOTTOMATRICE
      }
    }

    dataStorage[["structures"]]<<-listaROI
  }
  #####################################################################################
  # getStructuresFromXML: carica il file xml del RT struct
  #
  # INPUT:
  #   - fileName: il nome del file del RT struct
  # OUTPUT:
  #   -
  #################################################################################
  getStructuresFromXML<-function( fileName ) {

    obj.S<-services();
    massimo<-0
    folderCleanUp <- internalAttributes$attr_folderCleanUp
    arr.SeriesInstanceUID <- giveBackImageSeriesInstanceUID()

    # Load the XML file if not in cache
    doc <- obj.S$getXMLStructureFromDICOMFile( fileName = fileName, folderCleanUp = folderCleanUp )

    # prima di tutto controlla che la FrameOfReferenceUID sia la stessa OVUNQUE e che punti
    # ad una serie di immagini ESISTENTE!
    # E' un chiodo ma .... ragionevole, almeno per ora
    # Estrae la seriesIstanceUID, la frame of reference UID e il Referenced Frame of Reference UID
    RTStructSeriesInstanceUID <- xpathApply(doc,'//element[@tag="0020,000e" and @name="SeriesInstanceUID"]',xmlValue)[[1]]
    FORUID.m <- xpathApply(doc,'//element[@tag="0020,0052" and @name="FrameOfReferenceUID"]',xmlValue)[[1]]
    FORUID.d <- xpathApply(doc,'//element[@tag="3006,0024" and @name="ReferencedFrameOfReferenceUID"]',xmlValue)
    # Analizza tutti i valori del Referenced frame of reference UID e controlla se tutti i valori coincidono con il
    # frame of reference UID
    for( FORUID.d_index in seq( 1, length(FORUID.d) ) ) {
      if( FORUID.d[[ FORUID.d_index ]] !=  FORUID.m ) {
        logObj$sendLog(  "FrameOfReferenceUID not aligned in RTStruct file" , "ERR" );
      }
    }

    if( length(unique(unlist(FORUID.d))) > 1 ) ogObj$sendLog(  "more than 1 FrameOfReferenceUID in RTStruct file" , "ERR" );
    if( (unique(unlist(FORUID.d)) == FORUID.m ) == FALSE) ogObj$sendLog(  "FrameOfReferenceUID not aligned (?) in RTStruct file" , "ERR" );
    referencedFORUID <-  unique(unlist(FORUID.d))[1]

    # Guarda se ci sono delle immagini con quel FrameOfReferenceUID
    immagini.associabili <- which(SOPClassUIDList[ ,"FrameOfReferenceUID"] == FORUID.m & SOPClassUIDList[ ,"type"] =="IMG")
    if( length(immagini.associabili) == 0 ) {
      logObj$sendLog(  "the FrameOfReferenceUID of the RTStruct is not associated to an image" , "ERR"  );
    }

    # SEQUENCES: the one with the attribute  tag="3006,0020"  and name="StructureSetROISequence"
    # is the one with association NAME<->ID
    # Estrazione della parte di xml contenente il nome,numero delle varie ROI
    n2XML<-getNodeSet(doc,'/file-format/data-set/sequence[@tag="3006,0020" and @name="StructureSetROISequence"]/item')

    # SEQUENCES: now get the true coords
    # Estrazione delle coordinate di ciascuna ROI
    n3XML<-getNodeSet(doc,'/file-format/data-set/sequence[@tag="3006,0039" and @name="ROIContourSequence"]/item')

    # ROI Names
    matrice2<-c()
    # Per ciascuna ROI viene estratto il numero, il nome e organizzati in una matrice
    for(i in n2XML) {
      ROINumber<-xpathApply(xmlDoc(i),'/item/element[@tag="3006,0022"]',xmlValue)[[1]]
      ROIName<-xpathApply(xmlDoc(i),'/item/element[@tag="3006,0026"]',xmlValue)[[1]]
      matrice2<-rbind(matrice2,c(ROINumber,ROIName))
    }

    matrice2<-t(matrice2)
    # ROI Point list
    massimo<-0
    matrice3<-c()

    # esegue una interazione per ogni ROI
    # Non fa altro che estrarre da 'i' la parte che contiene i vertici della ROI.
    for(i in n3XML) {
      ROINumber<-xpathApply(xmlDoc(i),'/item/element[@tag="3006,0084"]',xmlValue)[[1]]
      ROIName<-matrice2[2,which(matrice2[1,]==ROINumber)]

      # la funzione getNodeSet() restituisce l'insieme delle coordinate dei vertici di una sola ROI
      # (sto intendendo per ROI l'insieme di tutti i contorni con lo stesso ROIName)
      listaPuntiDaRavanare<-getNodeSet(xmlDoc(i),'/item/sequence/item')
      numero.Punti.semiperimetro.massimo<-0

      # logObj$sendLog( paste( c( "\n Loading the ROI '".ROIName."' : ",length(listaPuntiDaRavanare)," polylines ") , collapse = '') )

      # Estrae solo un punto (fPoint.x, fPoint.y, fPoint.z) da listaPuntiDaRavanare in quanto viene assunto
      # che la ROI rispetto alla immagine sono tra loro paralleli (non sempre necessariamente vero)
      for(i2 in listaPuntiDaRavanare)   {

        # ReferencedSOPInstanceUID<-xpathApply(xmlDoc(i2),'//element[@tag="0008,1155"]',xmlValue)[[1]]
        # salva le coordinate delle ROI in una lista
        ROIPointList<-xpathApply(xmlDoc(i2),'/item/element[@tag="3006,0050"]',xmlValue)
        # organizza le coordinate in un vettore
        splittedROIPointList<-as.numeric(strsplit(ROIPointList[[1]],split = "\\\\")[[1]])
        # Separa in vettori diversi le coordinate dei punti x, y, z
        fPoint.x<-splittedROIPointList[1]
        fPoint.y<-splittedROIPointList[2]
        fPoint.z<-splittedROIPointList[3]

        # Calcola di quanti punti consta il "semiperimetro" della ROI
        # (questo serve per avere un'idea di dove prendere 3 punti abbstanza distanti per calcolare
        # il piano su cui giace)
        numero.Punti.semiperimetro <- as.integer((length(splittedROIPointList)/3)/2)

        # Se sei in presenza del numero di punti massimo, visto finora, prendi 3 punti
        if(numero.Punti.semiperimetro>numero.Punti.semiperimetro.massimo &
           numero.Punti.semiperimetro>10) {

          # prendi i tre punti sperabilmente più "distanti"
          # E' una stima, lo sa il cielo quali siano in realtà: dovrei
          # calcolare tutte le distanze reciproche! (ma anche no...)
          # Considero il primo punto della sequenza, quello ad un quarto e quello a metà
          p1 <- 1;    p2 <- numero.Punti.semiperimetro/2;      p3 <- numero.Punti.semiperimetro;
          # Costruisci la matrice di tutti i punti (così è più facile estrarre i 3 punti)
          matrice.punti <- matrix(splittedROIPointList,ncol=3,byrow = TRUE)
          p1 <- matrice.punti[p1,];
          p2 <- matrice.punti[p2,];
          p3 <- matrice.punti[p3,];
          # memorizza l'equazione del piano ed il numero di punti massimo della ROI
          numero.Punti.semiperimetro.massimo<-numero.Punti.semiperimetro
        }

        arr.ReferencedSOPInstanceUID<-c()
        if( length(arr.SeriesInstanceUID) == 0 )  {
          logObj$sendLog( "giveBackImageSeriesInstanceUID(); gave back nothing" , "ERR" );
        }

        # Cerca di assegnarlo ad una slice di immagine
        # Considerando un punto nella matrice della ROI calcola per ogni slice la distanza punto piano
        # in modo che se la distanza risulta minore di 0.2 assegna la ROI a quella determinata slice
        # cicla su tutte le immagini con quel FrameOfReferenceUID e cerca di capire a quale potresti associarla
        righe <- which(SOPClassUIDList[,"FrameOfReferenceUID"] == referencedFORUID & SOPClassUIDList[,"type"] == "IMG")
        for( riga in righe ) {
          tmp.SerInstUID <- SOPClassUIDList[riga,"SeriesInstanceUID"]
          slice.index <- as.character(SOPClassUIDList[riga,"SOPInstanceUID"])
          distanza <- obj.S$getPointPlaneDistance( c( fPoint.x, fPoint.y, fPoint.z ), dataStorage$info[[tmp.SerInstUID]][[slice.index]]$planeEquation)
          if( abs( distanza ) < internalAttributes$maxDistanceForImageROICoupling ) {
            arr.ReferencedSOPInstanceUID <- c( arr.ReferencedSOPInstanceUID, slice.index )
          }
        }

        # Assumento le slice di immagine fra loro parallele, vedi se è su un piano parallelo ad esse
        # se 'numero.Punti.semiperimetro.massimo'>0 significa che i tre punti della ROI sono stati estratti
        RSOPIUID.da.rimuovere <- c()
        for( RSOPIUID in arr.ReferencedSOPInstanceUID) {
          tmp.SerInstUID <- SOPClassUIDList[which(SOPClassUIDList[,"SOPInstanceUID"] == RSOPIUID),"SeriesInstanceUID"]

          if(numero.Punti.semiperimetro.massimo>0) {
            distanza1<-obj.S$getPointPlaneDistance(p1,dataStorage$info[[tmp.SerInstUID]][[RSOPIUID]]$planeEquation)
            distanza2<-obj.S$getPointPlaneDistance(p2,dataStorage$info[[tmp.SerInstUID]][[RSOPIUID]]$planeEquation)
            distanza3<-obj.S$getPointPlaneDistance(p3,dataStorage$info[[tmp.SerInstUID]][[RSOPIUID]]$planeEquation)
            # calcola il massimo gap fra i 3 punti.
            tot <- c(distanza1,distanza2,distanza3)
            tot <- max(tot)-min(tot)
            # Se è maggiore di un errore indicato ( .1 mm ) dichiara non paralleli i piani
            # (in realtà punti troppo vicini potrebbero fregarmi. Speriamo di no!)
            if(tot > .1 ) {
              # annovera fra le SOPInstanceUID da togliere dalla lista delle associabili
              logObj$sendLog ( paste( c("\n ROI ",ROIName," is not coplanar with images)" ),collapse = '')  );
              RSOPIUID.da.rimuovere <- c(RSOPIUID.da.rimuovere,RSOPIUID)
            }
          }
        }
        # Trattieni solo le SOPClassUID delle immagini che sono risultate paralele
        # (togli quelle non paralle dalla lista costruita precedentemente)
        arr.ReferencedSOPInstanceUID <- arr.ReferencedSOPInstanceUID[  which(!(arr.ReferencedSOPInstanceUID %in% RSOPIUID.da.rimuovere)) ]

        # Aggiorna la tabella delle associazioni ROI, immagini
        if( length(arr.ReferencedSOPInstanceUID) > 0 ) {
          for( tmp.RSOPUID in arr.ReferencedSOPInstanceUID ) {
            matrice3 <- rbind( matrice3, unlist(c( ROINumber, ROIName, ROIPointList, tmp.RSOPUID )) )
          }
        } else {
          logObj$sendLog (   paste(c("ROI ",ROIName,": the point (",fPoint.x,",",fPoint.y,",",fPoint.z,") has no image slice! "),collapse='')  );
        }
        # matrice3<-rbind(matrice3,c(ROINumber,ROIName,ROIPointList,ReferencedSOPInstanceUID))
      }
    }

    if( length(matrice3)>0) {
      if(is.matrix(matrice3)) {
        colnames(matrice3)<-c( "ROINumber", "ROIName", "ROIPointList", "ReferencedSOPInstanceUID"  )
      } else {
        names(matrice3)<-c( "ROINumber", "ROIName", "ROIPointList", "ReferencedSOPInstanceUID"  )
      }
    }

    colnames(matrice3) <- c("ROINumber", "ROIName", "ROIPointList", "ReferencedSOPInstanceUID")
    global_tableROIPointList <<- matrice3
    return(list("IDROINameAssociation"=matrice2,"tableROIPointList"=matrice3,"FORUID.m"=FORUID.m,"RTStructSeriesInstanceUID"=RTStructSeriesInstanceUID))
  }
  #=================================================================================
  # getFolderContent
  # Entra in una cartella ed estrai la lista dei file DICOM al suo interno, ricavandone
  # il tipo. Inoltre filtra gli oggetti DICOM solo alle SOPClassUID note.
  #=================================================================================
  getFolderContent <- function( pathToOpen ) {
    objS <- services()
    # if no path is given, use the set one
    if(!dir.exists(pathToOpen)) logObj$sendLog(  "The indicate Path does not exist" ,"ERR" );

    # salva in un array tutti i DICOM presenti nella cartella
    DCMFilenameArray<-list.files(pathToOpen,internalAttributes$defaultExtension.dicom)
    NIFTIFilenameArray<-list.files(pathToOpen,internalAttributes$defaultExtension.nifti)

    # lista con la SOP Class UID di ciascun DICOM
    SOPClassUIDList<-list()
    nomiColonne <- c("fileName","tag","kind","type","IPP.x","IPP.y","IPP.z","FrameOfReferenceUID","ImageOrder","field2Order","p.x","p.y","p.z","SOPInstanceUID","ImageOrientationPatient","SeriesInstanceUID")
    MMatrix <- matrix(ncol=length(nomiColonne),nrow=0)
    colnames(MMatrix)<-nomiColonne

    ImagingPositionArray <- c()
    Iteration <- 0
    ImageOrder <- 1
    for(i in 1:length(DCMFilenameArray) ) {
      fileNameWithPath<-paste(pathToOpen,"/",DCMFilenameArray[i] , sep="")
      # Nel caso in cui sia DICOM
      # if( substr(fileNameWithPath,nchar(fileNameWithPath)-3,nchar(fileNameWithPath))=='.dcm' ) {
      if( substr( fileNameWithPath, nchar( fileNameWithPath ) - 3,nchar(fileNameWithPath))!='.xml' &
          substr( fileNameWithPath, nchar( fileNameWithPath ) - 3,nchar(fileNameWithPath))!='.raw' &
          substr( fileNameWithPath, nchar( fileNameWithPath ) - (str_length(internalAttributes$defaultExtension.nifti)-1),nchar(fileNameWithPath))!=internalAttributes$defaultExtension.nifti
          ) {
# browser()
        if( internalAttributes$verbose == TRUE ) cat(".")

        riga <- nrow(MMatrix) + 1
        MMatrix <- rbind(MMatrix, rep("",ncol(MMatrix))  )
        valore<-getDICOMTag( fileName = fileNameWithPath, tag = "0008,0016")
        MMatrix[riga, "fileName"] <- fileNameWithPath
        MMatrix[riga, "tag"] <- valore
        FrameOfReferenceUID<-getDICOMTag( fileName = fileNameWithPath, tag = "0020,0052")

        MMatrix[riga, "FrameOfReferenceUID"]<-FrameOfReferenceUID
        MMatrix[riga, "kind"]<-"Unknown"
        MMatrix[riga, "type"]<-"Unknown"
        if( valore == "1.2.840.10008.5.1.4.1.1.2" ) MMatrix[riga, "kind"]<-"CTImageStorage"
        if( valore == "1.2.840.10008.5.1.4.1.1.481.2" ) MMatrix[riga, "kind"]<-"RTDoseStorage"
        if( valore == "1.2.840.10008.5.1.4.1.1.481.3" ) MMatrix[riga, "kind"]<-"RTStructureSetStorage"
        if( valore == "1.2.840.10008.5.1.4.1.1.481.3" ) MMatrix[riga, "type"]<-"RTStruct"
        if( valore == "1.2.840.10008.5.1.4.1.1.481.5" ) MMatrix[riga, "kind"]<-"RTPlanStorage"
        if( valore == "1.2.840.10008.5.1.4.1.1.4" ) MMatrix[riga, "kind"]<-"MRImageStorage"
        if( valore == "1.2.840.10008.5.1.4.1.1.128" ) MMatrix[riga, "kind"]<-"PositronEmissionTomographyImageStorage"
        if( valore == "1.2.840.10008.5.1.4.1.1.2.1" ) MMatrix[riga, "kind"]<-"CTImageStorage"

        if( MMatrix[riga, "kind"] == "CTImageStorage" |
            MMatrix[riga, "kind"] == "MRImageStorage" |
            MMatrix[riga, "kind"] == "PositronEmissionTomographyImageStorage" )  {
          MMatrix[riga, "type"]<-"IMG"
          MMatrix[riga, "ImageOrder"]<-ImageOrder
          ImageOrder <- ImageOrder + 1

          # Prendi il Pixel Spacing e lo slice thickness
          pixelSpacing<-objS$splittaTAG(getDICOMTag(fileName = fileNameWithPath,tag = "0028,0030"))
          sliceThickness<-objS$splittaTAG(getDICOMTag(fileName = fileNameWithPath,tag = "0018,0050"))
          MMatrix[riga, "p.x"]<-pixelSpacing[1]
          MMatrix[riga, "p.y"]<-pixelSpacing[2]
          MMatrix[riga, "p.z"]<-sliceThickness

          # ImageOrientation
          ImageOrientationPatient <- getDICOMTag(fileName = fileNameWithPath , tag = "0020,0037")
          MMatrix[riga, "ImageOrientationPatient"] <- ImageOrientationPatient
        }

        SOPInstanceUID <- getDICOMTag( tag = "0008,0018", fileName = fileNameWithPath)
        MMatrix[riga, "SOPInstanceUID"] <- SOPInstanceUID

        ImagePositionPatient <- getDICOMTag( fileName = fileNameWithPath, tag = "0020,0032")
        if( !is.na(ImagePositionPatient) ) {
          ImagingPosition <- objS$splittaTAG(ImagePositionPatient)
          MMatrix[riga, "IPP.x"] <- ImagingPosition[1]
          MMatrix[riga, "IPP.y"] <- ImagingPosition[2]
          MMatrix[riga, "IPP.z"] <- ImagingPosition[3]
        }
        SeriesInstanceUID<-getDICOMTag(tag = "0020,000e", fileName = fileNameWithPath )
        MMatrix[riga, "SeriesInstanceUID"] <- SeriesInstanceUID
      }
    }
    # Ora gestisti eventuali NIFTI
    for(i in 1:length(NIFTIFilenameArray) ) {
      fileNameWithPath<-paste(pathToOpen,"/",NIFTIFilenameArray[i] , sep="")
      # if( substr(fileNameWithPath,nchar(fileNameWithPath)-str_length("nii.gz"),nchar(fileNameWithPath))=='.nii.gz' ) {
      if(substr( fileNameWithPath, nchar( fileNameWithPath ) - (str_length(internalAttributes$defaultExtension.nifti)-1),nchar(fileNameWithPath))==internalAttributes$defaultExtension.nifti) {
        riga <- nrow(MMatrix)+1
        MMatrix <- rbind(MMatrix, rep("",ncol(MMatrix))  )
        MMatrix[riga, "fileName"] <- fileNameWithPath
        MMatrix[riga, "tag"] <- ""
        MMatrix[riga, "kind"] <- "nifti"
        MMatrix[riga, "type"] <- "nifti"
      }
    }

    # Ora devi ordinare le immagini in funzione delle loro coordinate!
    arr.SeriesInstanceUID <- unique(MMatrix[  which( MMatrix[,"type"]=="IMG" ), "SeriesInstanceUID" ])
    for( SIUID in arr.SeriesInstanceUID ) {
      righe <- which( MMatrix[,"type"]=="IMG" & MMatrix[,"SeriesInstanceUID"]==SIUID )
      sottomatrice <- MMatrix[  righe,  c("SOPInstanceUID","IPP.x","IPP.y","IPP.z") ]
      FOV.z <- diff(range(as.numeric(sottomatrice[,"IPP.z"])))
      FOV.y <- diff(range(as.numeric(sottomatrice[,"IPP.y"])))
      FOV.x <- diff(range(as.numeric(sottomatrice[,"IPP.x"])))

      # Ordinale per il gap maggiore (para-assiale/coronale/saggittale)
      winner <- order(c(FOV.x,FOV.y,FOV.z),decreasing = T)[1]
      if( winner == 1 ) field2Order <- "IPP.x"
      if( winner == 2 ) field2Order <- "IPP.y"
      if( winner == 3 ) field2Order <- "IPP.z"
      nuovo.ordine <- order(as.numeric(sottomatrice[,field2Order]),decreasing = F)
      tmp.SOPIUID <- sottomatrice[nuovo.ordine,"SOPInstanceUID"]
      for( i in 1:length(tmp.SOPIUID )) {
        MMatrix[ which( MMatrix[,"SOPInstanceUID"]==tmp.SOPIUID[i]  ) , "ImageOrder"] <- i
        MMatrix[ which( MMatrix[,"SOPInstanceUID"]==tmp.SOPIUID[i]  ) , "field2Order"] <- field2Order
      }
    }

    SOPClassUIDList <- MMatrix
    return(SOPClassUIDList);
  }
  #=================================================================================
  # loadCTRMRDNScans
  # Loads a DICOM CT/MR Scans
  #=================================================================================
  loadCTRMNRDScans<-function( ) {
    objS <- services()
    # Cicla sulle sole righe corrispondenti a delle immagini
    righe.con.immagini <- which(SOPClassUIDList[,"type"]=="IMG")

    for( riga in righe.con.immagini ) {

      fileName <- SOPClassUIDList[riga,"fileName"]
      SeriesInstanceUID<-SOPClassUIDList[riga, "SeriesInstanceUID"]

      if( internalAttributes$verbose == TRUE ) cat("\n\t ",fileName)

      FrameOfReferenceUID<-SOPClassUIDList[riga,"FrameOfReferenceUID"]
      ImageOrder <- SOPClassUIDList[riga,"ImageOrder"]
      SOPInstanceUID <- SOPClassUIDList[riga,"SOPInstanceUID"]
      kind.of.SOPClassUID <- SOPClassUIDList[riga, "kind"]
      IPP.x <- SOPClassUIDList[riga,"IPP.x"]
      IPP.y <- SOPClassUIDList[riga,"IPP.y"]
      IPP.z <- SOPClassUIDList[riga,"IPP.z"]
      ImageOrientationPatient <- SOPClassUIDList[riga,"ImageOrientationPatient"]
      ImagePositionPatient <- c(SOPClassUIDList[riga,"IPP.x"],SOPClassUIDList[riga,"IPP.y"],SOPClassUIDList[riga,"IPP.z"])
      pixelSpacing <- as.numeric(c(SOPClassUIDList[riga,"p.x"],SOPClassUIDList[riga,"p.y"]))

      # Se non era ancora stato settato, setta il FrameOfReferenceUID di riferimento
      if(is.na(internalAttributes$attr_mainFrameOfReferenceUID)) internalAttributes$attr_mainFrameOfReferenceUID <<- FrameOfReferenceUID

      # Verifica che il FORUID sia compatibile con quello caricato (se no, dai errore)
      if( FrameOfReferenceUID != internalAttributes$attr_mainFrameOfReferenceUID ) stop("FORUID differente!")

      # if( SeriesInstanceUID %in% names(dataStorage[["info"]]) )  dataStorage[["info"]][[SeriesInstanceUID]] <- list()
      # if( SeriesInstanceUID %in% names(dataStorage[["info"]]) )
      if( !(SeriesInstanceUID %in% names(dataStorage[["info"]])) )  dataStorage[["info"]][[SeriesInstanceUID]] <<- list()

      dataStorage[["info"]][[SeriesInstanceUID]][[SOPInstanceUID]]<<-list()
      dataStorage[["info"]][[SeriesInstanceUID]][[SOPInstanceUID]][["PatientPosition"]]<<-c(SOPClassUIDList[riga,"IPP.x"],SOPClassUIDList[riga,"IPP.y"],SOPClassUIDList[riga,"IPP.z"])

      # Carica l'immagine
      SingleSliceLoader <- loadCTRMNRDScans.SingleSlice( fileName = fileName )
      immagine <- SingleSliceLoader$immagine

      # now update the structure in memory
      if( length( dataStorage$img ) == 0 ) dataStorage$img <<- list()
      if( length( dataStorage$img[[SeriesInstanceUID]] ) == 0 ) dataStorage$img[[ SeriesInstanceUID ]] <<- list()
      dataStorage$img[[ SeriesInstanceUID ]][[ SOPInstanceUID ]] <<- immagine

      # Costruisci la matrice per le trasformazioni affini
      iPP<-as.numeric(c(IPP.x,IPP.y,IPP.z))
      iOP<-objS$splittaTAG( ImageOrientationPatient )
      oM<-matrix(c(iOP[1],iOP[2],iOP[3],0,iOP[4],iOP[5],iOP[6],0,0,0,0,0,iPP[1],iPP[2],iPP[3],1),ncol=4);

      dataStorage[[ "info" ]][[ SeriesInstanceUID ]][[ SOPInstanceUID ]][[ "ImagePositionPatient" ]] <<- ImagePositionPatient
      dataStorage[[ "info" ]][[ SeriesInstanceUID ]][[ SOPInstanceUID ]][[ "ImageOrientationPatient" ]] <<- ImageOrientationPatient
      oM[1,1]<-oM[1,1]*pixelSpacing[1]
      oM[2,1]<-oM[2,1]*pixelSpacing[1]
      oM[3,1]<-oM[3,1]*pixelSpacing[1]
      oM[1,2]<-oM[1,2]*pixelSpacing[2]
      oM[2,2]<-oM[2,2]*pixelSpacing[2]
      oM[3,2]<-oM[3,2]*pixelSpacing[2]
      dataStorage[["info"]][[SeriesInstanceUID]][[SOPInstanceUID]][["orientationMatrix"]]<<-oM
      dataStorage[["info"]][[SeriesInstanceUID]][[SOPInstanceUID]][["pixelSpacing"]]<<-pixelSpacing
      dataStorage[["info"]][[SeriesInstanceUID]][[SOPInstanceUID]][["Rows"]]<<-getDICOMTag(tag = "0028,0010", fileName = fileName)
      dataStorage[["info"]][[SeriesInstanceUID]][[SOPInstanceUID]][["Columns"]]<<-getDICOMTag(tag = "0028,0011", fileName = fileName)
      if( kind.of.SOPClassUID == "MRImageStorage" ) {
        dataStorage[["info"]][[SeriesInstanceUID]][[SOPInstanceUID]][["RepetitionTime"]]<<-getDICOMTag( tag = "0018,0080", fileName = fileName)
        dataStorage[["info"]][[SeriesInstanceUID]][[SOPInstanceUID]][["EchoTime"]]<<-getDICOMTag(tag = "0018,0081", fileName = fileName)
      }
      # Aggiungi eventuali campi specifici per quel tipo di imaging (vero sopratutto per le PET/CT)
      for( campo in names(SingleSliceLoader$fields)) {
        dataStorage[[ "info" ]][[ SeriesInstanceUID ]][[ SOPInstanceUID ]][[ campo ]] <<- SingleSliceLoader$fields[[ campo ]]
      }

      if( kind.of.SOPClassUID == "PositronEmissionTomographyImageStorage" ) {
        # Fai qualche controllo di qualita', sulla reale possibilita' di importare le immagini PET
        rescale.type <-  SingleSliceLoader$fields$rescale.type
        UM <-  SingleSliceLoader$fields$UM
        CountsSource <- SingleSliceLoader$fields$CountsSource
        DecayCorrection <- SingleSliceLoader$fields$DecayCorrection
        deltaT <- SingleSliceLoader$fields$deltaT

        if(is.na(rescale.type)) rescale.type <- UM
        if( UM != rescale.type) {  logObj$sendLog(  "in PET image the rescale slope/intercept have different UM than the one used in the image (0054,1001) vs (0028,1054)" , "ERR"  )  }
        if(CountsSource!='EMISSION' | DecayCorrection!='START') { logObj$sendLog(  c("\n ERROR: CountsSource!='EMISSION' or DecayCorrection!='START' ! This modality is not yet supported") , "ERR")  }
        if(is.na(deltaT) | is.null(deltaT) | deltaT==0) {  logObj$sendLog( "\n Error: deltaT between RadiopharmaceuticalStartTime and AcquisitionTime seems to be invalid" , "ERR" )  }
      }

      # three points to find out plane equation
      Pa<-c(oM[1,4],oM[2,4],oM[3,4])
      Pb<-objS$get3DPosFromNxNy(1000,0,oM)
      Pc<-objS$get3DPosFromNxNy(0,1000,oM)

      abcd<-objS$getPlaneEquationBetween3Points(Pa,Pb,Pc)
      piano<-matrix(abcd,nrow=1)
      colnames(piano)<-c("a","b","c","d")
      dataStorage[["info"]][[SeriesInstanceUID]][[SOPInstanceUID]][["planeEquation"]]<<-piano

      if(  kind.of.SOPClassUID == "RTDoseStorage" ) {
        stop("not yet implemented")
      }
    }
  }

  #=================================================================================
  # loadCTRMNRDScans.SingleSlice
  # Si occupa del caricamento di una singola slice. E' stato scorporato in quanto deve
  # essere potenzialmente evocato da più punti del programma ( gestione cache )
  #=================================================================================
  loadCTRMNRDScans.SingleSlice<-function( fileName ) {

    objServ<-services()
    fields <- list()

    # get the image data
    immagine<-getDICOMTag(tag = "7fe0,0010", fileName = fileName);
    SOPClassUID <- SOPClassUIDList[ SOPClassUIDList[,"fileName"]==fileName, "kind"]

    # apply rescaleSlope and rescaleIntercept, if needed
    rescale.intercept<-as.numeric(getDICOMTag(fileName = fileName,tag ="0028,1052" )); # rescale Intercept
    rescale.slope<-as.numeric(getDICOMTag(fileName = fileName,tag ="0028,1053" )); # rescale Slope
    rescale.type<-getDICOMTag(fileName = fileName,tag ="0028,1054" ); # rescale Type

    if(is.na(rescale.intercept)) rescale.intercept = 0;
    if(is.na(rescale.slope)) rescale.slope = 1;
    immagine <- immagine * rescale.slope + rescale.intercept

    immagine <- objServ$rotateMatrix( immagine, rotations = 1 )

    if( SOPClassUID == "PositronEmissionTomographyImageStorage" ) {
      res <- calculate.SUVCoefficient.BW( fileName = fileName)
      SUVCoefficient.BW <- res$SUVCoefficient.BW
      immagine <- SUVCoefficient.BW * immagine

      # copia i campi 'fields' acquisisti dalla funzione per il calcolo del SUV e copiali nei campi
      # fields della lista che verra' restituita
      for( campo in names(res$fields)) {
        fields[[ campo ]] <- res$fields[[ campo ]]
      }
    }

    fields$rescale.intercept <- rescale.intercept
    fields$rescale.slope <- rescale.slope
    fields$rescale.type <- rescale.type

    return( list( "immagine" = immagine,
                  "fields" = fields ) )
  }
  #=================================================================================
  # calculate.SUVCoefficient.BW
  # Calcola il coefficiente moltiplicativo per il SUV
  #=================================================================================
  calculate.SUVCoefficient.BW<-function( fileName ) {
    fields <- list()
    AcquisitionTime<-getDICOMTag(fileName = fileName,tag ="0008,0032" );
    RadiopharmaceuticalStartTime<-getDICOMTag(fileName = fileName,tag ="0018,1072" );
    PatientWeight<-as.numeric(getDICOMTag(fileName = fileName,tag ="0010,1030" ));
    RadionuclideTotalDose<-as.numeric(getDICOMTag(fileName = fileName,tag ="0018,1074" )); # RadionuclideTotalDose
    RadionuclideHalfLife<-as.numeric(getDICOMTag(fileName = fileName,tag ="0018,1075" )); # RadionuclideHalfLife
    UM<-getDICOMTag(fileName = fileName, tag ="0054,1001" ); # UM of voxel Cube
    CountsSource<-getDICOMTag(fileName = fileName, tag ="0054,1002" ); # CountsSource
    DecayCorrection<-getDICOMTag(fileName = fileName, tag ="0054,1102" ); # DecayCorrection
    rescale.type<-getDICOMTag(fileName = fileName,tag ="0028,1054" ); # rescale Type)

    # deltaT<-as.numeric(difftime(as.POSIXct(AcquisitionTime, format = "%H:%M:%S"),as.POSIXct(RadiopharmaceuticalStartTime, format = "%H:%M:%S"),units = 'secs'))
    deltaT<-as.numeric(difftime(as.POSIXct(AcquisitionTime, format = "%H%M%S"),as.POSIXct(RadiopharmaceuticalStartTime, format = "%H%M%S"),units = 'secs'))

    # Una bella considerazione sul rescale.type ci potrebbe anche stare. Nel frattempo pongo il rescale a 1
    rescaleDueToUM<-1

    #SUVCoefficient.BW<-PatientWeight/( RadionuclideTotalDose * exp( -deltaT *log(2)/(RadionuclideHalfLife) ) )
    SUVCoefficient.BW<-PatientWeight/( RadionuclideTotalDose * 2^( -deltaT / RadionuclideHalfLife ) ) * 1000
    SUVCoefficient.BW<-SUVCoefficient.BW*rescaleDueToUM

    fields$AcquisitionTime <- AcquisitionTime
    fields$RadiopharmaceuticalStartTime <- RadiopharmaceuticalStartTime
    fields$PatientWeight <- PatientWeight
    fields$RadionuclideTotalDose <- RadionuclideTotalDose
    fields$RadionuclideHalfLife <- RadionuclideHalfLife
    fields$UM <- UM
    fields$CountsSource <- CountsSource
    fields$DecayCorrection <- DecayCorrection
    fields$rescale.type <- rescale.type
    fields$deltaT <- deltaT
    fields$SUVCoefficient.BW <- SUVCoefficient.BW

    return( list( "SUVCoefficient.BW" = SUVCoefficient.BW,
                  "fields" = fields ) );
  }
  #=================================================================================
  # getImageFromRAW
  # build a row data from a DICOM file stored on filesystem and load it
  # into memory (using DCMTK)
  #=================================================================================
  getImageFromRAW<-function(fileName) {

    objSV<-services()
    fileNameRAW<-paste(fileName,".0.raw")
    fileNameRAW<-str_replace_all(string = fileNameRAW , pattern = " .0.raw",replacement = ".0.raw")

    if(!file.exists(fileName)) logObj$sendLog( " the fileName is missing in geoLet::getImageFromRAW()", "ERR"  );

    pathToStore<-substr(fileName,1,tail(which(strsplit(fileName, '')[[1]]=='/'),1)-1)
    if(!file.exists( fileNameRAW )  | internalAttributes$attr_folderCleanUp==TRUE) {
      stringa1<-"dcmdump";
      if ( Sys.info()["sysname"] == "Windows") {
        fileNameFS<-chartr("\\","/",fileName);
        stringa2<-chartr("/","\\\\",stringa1)
      }
      else fileNameFS<-fileName;
      fileNameFS <- paste(c("'",fileNameFS,"'"),collapse = '')
      stringa2<-paste(" +W  ",pathToStore,fileNameFS,collapse='')
      options(warn=-1)
      stringone<-as.character(paste( c(stringa1," ",stringa2),collapse=''))

      # gestisci le system call in maniera diversa in funzione che sia WINDOWS o LINUX
      if ( Sys.info()["sysname"] == "Windows") {
        res<-.C("executeCMDLine",  as.character(stringone), as.integer(str_length(stringone))  )
      }
      else {
        system2(stringa1,stringa2,stdout=NULL)
      }
      options(warn=0)
    }
    rowsDICOM<-as.numeric(getDICOMTag(fileName = fileName,tag = '0028,0010'))
    columnsDICOM<-as.numeric(getDICOMTag(fileName = fileName,tag = '0028,0011'))
    bitsAllocated<-as.numeric(getDICOMTag(fileName = fileName,tag = '0028,0100'))
    SOPClassUID <- SOPClassUIDList[ SOPClassUIDList[,"fileName"]==fileName, "kind"]

    if( SOPClassUID != "RTDoseStorage" ) {
      if(bitsAllocated!=16)
        logObj$sendLog( "16bit pixel are allowed only for non-RTDoseStorage", "ERR"  );
      if ( Sys.info()["sysname"] == "Windows") {
        fileNameRAWFS<-chartr("\\","/",fileNameRAW);
        fileNameRAWFS<-chartr("/","\\\\",fileNameRAWFS);
      }
      else fileNameRAWFS<-fileNameRAW;

      if(!file.exists(fileNameRAWFS)) logObj$sendLog( "problem in creating image binary file in geoLet::getImageFromRAW()", "ERR"  );

      rn<-readBin(con = fileNameRAWFS, what="integer", size=2, endian="little",n=rowsDICOM*columnsDICOM)
      rn<-matrix(rn,ncol=columnsDICOM, byrow = TRUE)
    }
    if( SOPClassUID == "RTDoseStorage" ) {
      if(bitsAllocated==32) {
        if( SOPClassUID != "RTDoseStorage" )
          logObj$sendLog(  "32bit pixel are allowed only for RTDoseStorage" ,"ERR" );
        numberOfFrames<-as.numeric(getDICOMTag(fileName = fileName, tag = '0028,0008'))
        if ( Sys.info()["sysname"] == "Windows") {
          fileNameRAWFS<-chartr("\\","/",fileNameRAW);
          fileNameRAWFS<-chartr("/","\\\\",fileNameRAWFS);
        }
        else fileNameRAWFS<-fileNameRAW;

        if(!file.exists(fileNameRAWFS)) logObj$sendLog( "problem in creating image binary file in geoLet::getImageFromRAW()" ,"ERR" );

        rn<-readBin(con = fileNameRAWFS, what="integer", size=4, endian="little",n=rowsDICOM*columnsDICOM*numberOfFrames)
        # per ora va via come ciclo FOR, poi ci ragioniamo.... (per le performances)
        matRN<-array(0,c(rowsDICOM,columnsDICOM,numberOfFrames))
        ct<-1
        for( z in seq(1,numberOfFrames)) {
          for(x in seq(1,rowsDICOM)) {
            for(y in seq(1,columnsDICOM)) {
              matRN[x,columnsDICOM-y,z]<-rn[ct]
              ct<-ct+1
            }
          }
        }
        new_atRN<-array(0,c(columnsDICOM,rowsDICOM,numberOfFrames))
        for(ct in seq(1:dim(matRN)[3]  )) {
          new_atRN[,,ct]<-t(objSV$rotateMatrix( matRN[,,ct], rotations=2 ))
        }
        rn<-new_atRN
      } else  {
        numberOfFrames<-as.numeric(getDICOMTag(fileName = fileName, tag = '0028,0008'))

        if ( Sys.info()["sysname"] == "Windows") {
          fileNameRAWFS<-chartr("\\","/",fileNameRAW);
          fileNameRAWFS<-chartr("/","\\\\",fileNameRAWFS)
        }
        else fileNameRAWFS<-fileNameRAW;

        if(!file.exists(fileNameRAWFS)) logObj$sendLog( "problem in creating image binary file in geoLet::getImageFromRAW()", "ERR"  );

        rn<-readBin(con = fileNameRAWFS, what="integer", size=2, endian="little",n=rowsDICOM*columnsDICOM*numberOfFrames)
        matRN<-array(0,c(rowsDICOM,columnsDICOM,numberOfFrames))
        ct<-1
        for( z in seq(1,numberOfFrames)) {
          for(x in seq(1,rowsDICOM)) {
            for(y in seq(1,columnsDICOM)) {
              matRN[x,columnsDICOM-y,z]<-rn[ct]
              ct<-ct+1
            }
          }
        }
        new_atRN<-array(0,c(columnsDICOM,rowsDICOM,numberOfFrames))
        for(ct in seq(1:dim(matRN)[3]  )) {
          new_atRN[,,ct]<-t(objSV$rotateMatrix( matRN[,,ct], rotations=2 ))
        }
        rn<-new_atRN
      }
    }

    return(rn)
  }
  #=================================================================================
  # getROIList
  # restituisce la lista delle ROI
  #=================================================================================
  getROIList<-function() {
    mat2Ret<-matrix( c(seq(1,length(names(dataStorage$structures))),names(dataStorage$structures)),nrow=2 ,byrow=T )
    return(mat2Ret[2,])
  }
  #=================================================================================
  # getImageVoxelCube
  # give back the greyLevel voxel cube. If no ps.x/y/z are specified it gives back
  # the voxelCube of the original dimensions, otherwise it gives back the interpolated
  # voxelCube according to the wished pixelSpacing along x,y or z
  #=================================================================================
  getImageVoxelCube<-function( ps.x=NA, ps.y=NA, ps.z=NA , SeriesInstanceUID = NA) {
    objS<-services();

    if( length(giveBackImageSeriesInstanceUID()) > 1 &
        is.na(SeriesInstanceUID) ) {
      logObj$sendLog(  "There are more than one Series, please specify which SeriesInstanceUID you want" ,"ERR" );
    }

    # prendi il cubone
    voxelCube<-createImageVoxelCube( SeriesInstanceUID = SeriesInstanceUID)

    # se non  server interpolare
    if(is.na(ps.x) && is.na(ps.y) && is.na(ps.z) ) return(voxelCube)

    # se invece serve interpolare: prendi i pixelSpacing lungo la X, la Y e la Z (slice thickness)
    oldPixelSpacing<-getPixelSpacing( SeriesInstanceUID = SeriesInstanceUID);

    if(is.na(ps.x))  ps.x <- oldPixelSpacing[1];
    if(is.na(ps.y))  ps.y <- oldPixelSpacing[2];
    if(is.na(ps.z))  ps.z <- oldPixelSpacing[3];

    voxelCube<-objS$new.trilinearInterpolator(
      voxelCube = voxelCube,
      pixelSpacing.new = c(ps.x,ps.y,ps.z),
      pixelSpacing.old = oldPixelSpacing )

    invisible(gc())

    return( voxelCube )
  }
  #=================================================================================
  # NAME: getROIVoxels
  # restituisce i voxel interni ad una data ROI
  #=================================================================================
  getROIVoxels<-function( Structure  , new.pixelSpacing=c(), SeriesInstanceUID = NA, croppedCube  = TRUE, onlyVoxelCube = FALSE, voxel.inclusion.threshold = 0.5) {

    if( dataStorage$info$structures[[Structure]]$type == "DICOMRTStruct" ) {
      res <- getROIVoxels.DICOM(Structure = Structure , new.pixelSpacing=new.pixelSpacing,
                                SeriesInstanceUID = SeriesInstanceUID, croppedCube  = croppedCube,
                                onlyVoxelCube = onlyVoxelCube)
    }
    if( dataStorage$info$structures[[Structure]]$type == "NIFTI" ) {
      res <- getROIVoxels.NIFTI(Structure = Structure , new.pixelSpacing=new.pixelSpacing,
                                SeriesInstanceUID = SeriesInstanceUID, croppedCube  = croppedCube,
                                onlyVoxelCube = onlyVoxelCube, voxel.inclusion.threshold = voxel.inclusion.threshold)
    }

    return( res )
  }
  #=================================================================================
  # NAME: getROIVoxels.NIFTI
  # restituisce i voxel interni ad una data ROI - DICOM
  #=================================================================================
  getROIVoxels.NIFTI<-function( Structure  , new.pixelSpacing=c(), SeriesInstanceUID = NA, croppedCube  = TRUE, onlyVoxelCube = FALSE, voxel.inclusion.threshold ) {
    objS<-services();
    # browser()

    if( length(new.pixelSpacing) > 1) stop("\n Interpolation not yet supported")

    fileNameWithPath <- dataStorage$info$structures[[ Structure ]]$fileName
    aaa <- readNIfTI(fname = fileNameWithPath)
    nifti.VC <- slot(aaa,".Data")
    if( slot(aaa,"reoriented") == FALSE ) logObj$sendLog(  "In the NIFTI file, 'reoriented' is set to FALSE" ,"ERR" );
    if( slot(aaa,"scl_slope") != 1 ) logObj$sendLog(  "In the NIFTI file, 'slope' is no 1" ,"ERR" );
    if( slot(aaa,"scl_inter") != 0 ) logObj$sendLog(  "In the NIFTI file, 'intercept' is no 1" ,"ERR" );
    dim.x <- slot(aaa,"dim_")[2];    dim.y <- slot(aaa,"dim_")[3];    dim.z <- slot(aaa,"dim_")[4]

    # slice <- 10;  nifti.VC[ which(nifti.VC==0) ] <- NA;    image(VC[,,slice]);    image(nifti.VC[,,35* slice/(dim(VC)[3]) ], add=T,col='green')
    VC <- getImageVoxelCube(SeriesInstanceUID = SeriesInstanceUID)
    masked.array <- array(0,dim = dim(VC) )

    if(  dim.x < dim(VC)[1]  | dim.y < dim(VC)[2] | dim.z < dim(VC)[3] ) stop("one of the dimensions of the nifti cube is less than the voxel cube, please check")

    matrice.Punti <- which( nifti.VC!=0,arr.ind = T )

    for( riga in 1:nrow(matrice.Punti )) {
      pos.old.VC <- c(  ( matrice.Punti[riga,1] / dim(nifti.VC)[1]) * dim( VC )[1],
                        ( matrice.Punti[riga,2] / dim(nifti.VC)[2]) * dim( VC )[2],
                        ( matrice.Punti[riga,3] )  )
      pos.old.VC <- round( pos.old.VC )
      masked.array[ pos.old.VC[1], pos.old.VC[2], pos.old.VC[3]  ] <- masked.array[ pos.old.VC[1], pos.old.VC[2], pos.old.VC[3]  ] + 1
    }

    # calcola il rapporto in volume dei voxels e guarda quali sono sopra o sotto soglia
    # rapporto <- round( (dim(VC)[1] * dim(VC)[2] * dim(VC)[3]) / ( dim.x * dim.y * dim.z ) * internalAttributes$threshold.4.niftiROI )
    rapporto <- round( ( dim.x * dim.y * dim.z )  / (dim(VC)[1] * dim(VC)[2] * dim(VC)[3])  * voxel.inclusion.threshold )
    # browser()

    masked.array[ which( masked.array < rapporto ) ] <- 0
    masked.array[ which( masked.array != 0 ) ] <- 1

    return(masked.array)
  }
  #=================================================================================
  # NAME: getROIVoxels.dicom
  # restituisce i voxel interni ad una data ROI - DICOM
  #=================================================================================
  getROIVoxels.DICOM<-function( Structure  , new.pixelSpacing=c(), SeriesInstanceUID = NA, croppedCube  = TRUE, onlyVoxelCube = FALSE) {
    objS<-services();
# browser()
    if(!(Structure %in% getROIList())) logObj$sendLog(  paste(c( Structure," not present."  ),collapse = ''), "ERR"  );
    if( length(SeriesInstanceUID) > 1 ) logObj$sendLog(  paste( "Error, too many SeriesInstanceUIDs. No more than one is admitted."  ), "ERR"  );

    # try to find out which Series is the CT/MR serie
    if(is.na(SeriesInstanceUID)) {
      arr.SeriesInstanceUID <- giveBackImageSeriesInstanceUID()
      if( length(arr.SeriesInstanceUID) == 0 ) {
        logObj$sendLog( "No image series seem to be available.. do you have any CT/MRI/PET in the indicated folder?" , "ERR" );
      }
      if( length(arr.SeriesInstanceUID) >1 ) {
        logObj$sendLog( "too many series are available: please pick up one" , "ERR" );
      }
      SeriesInstanceUID <- arr.SeriesInstanceUID[1]
    }

    # Questo va fatto solo se e' chiaro che non si tratta di un nifti'
    res <- getROIVoxelsFromCTRMN( Structure = Structure, SeriesInstanceUID = SeriesInstanceUID,
                                  new.pixelSpacing = new.pixelSpacing)

    if(sum(!is.na(res$masked.images)) == 0 ) return(list())

    croppedRes<-list()
    if( croppedCube == TRUE ) {
      croppedRes$DOM<-res$DOM
      croppedRes$geometricalInformationOfImages<-res$geometricalInformationOfImages
      croppedRes$masked.images<-objS$cropCube( bigCube = res$masked.images)
      croppedRes$masked.images$location$fe<-dim(res$masked.images)[1]
      croppedRes$masked.images$location$se<-dim(res$masked.images)[2]
      croppedRes$masked.images$location$te<-dim(res$masked.images)[3]
      croppedRes$geometricalInformationOfImages$koc<-"croppedCube"
      croppedRes$resamplingInformation<-res$resamplingInformation
    } else {
      croppedRes$DOM<-res$DOM
      croppedRes$geometricalInformationOfImages<-res$geometricalInformationOfImages
      croppedRes$masked.images$location$fe<-dim(res$masked.images)[1]
      croppedRes$masked.images$location$se<-dim(res$masked.images)[2]
      croppedRes$masked.images$location$te<-dim(res$masked.images)[3]
      croppedRes$geometricalInformationOfImages$koc<-"uncroppedCube"
      croppedRes$resamplingInformation<-res$resamplingInformation
      croppedRes$masked.images$location$min.x <- 1;      croppedRes$masked.images$location$min.y <- 1;      croppedRes$masked.images$location$min.z <- 1
      croppedRes$masked.images$location$max.x <- dim(res$masked.images)[1];
      croppedRes$masked.images$location$max.y <- dim(res$masked.images)[2];
      croppedRes$masked.images$location$max.z <- dim(res$masked.images)[3];
      croppedRes$masked.images$voxelCube <- res$masked.images
    }
    class(croppedRes)<-"geoLetStructureVoxelList"

    # Se si vuole indietro SOLO il VOXELCube, prevediamo un ritorno semplificato
    if( onlyVoxelCube == TRUE ) croppedRes <- croppedRes$masked.images$voxelCube

    invisible(gc())
    return( croppedRes )
  }
  #=================================================================================
  # NAME: getROIVoxelsFromCTRMN
  # Estrae i voxel da scansioni CT,MR
  #=================================================================================
  getROIVoxelsFromCTRMN<-function( Structure = Structure, SeriesInstanceUID = SeriesInstanceUID,
                                   new.pixelSpacing=c() ) {
    objService<-services()

    # define some variables to make more clear the code
    numberOfRows<-as.numeric(dataStorage$info[[SeriesInstanceUID]][[1]]$Rows);
    numberOfColumns<-as.numeric(dataStorage$info[[SeriesInstanceUID]][[1]]$Columns);
    numberOfSlices<-length(dataStorage$img[[SeriesInstanceUID]]);

    old.ps <- getPixelSpacing( SeriesInstanceUID =  SeriesInstanceUID)
    if( length(new.pixelSpacing)>0 ) {
      if(new.pixelSpacing[1] == old.ps[1] & new.pixelSpacing[2]==old.ps[2]) new.pixelSpacing<-c()
    }
    if(  length(new.pixelSpacing)>0 ){
      new.pixelSpacing <- c(new.pixelSpacing,old.ps[3])
      ratio.ps.x <- old.ps[1] / new.pixelSpacing[1];   ratio.ps.y <- old.ps[2] / new.pixelSpacing[2]
    } else {
      new.numberOfRows <- numberOfRows;   new.numberOfColumns <- numberOfColumns
    }

    # initialize the image array with the right dimension
    image.arr<-array( data = -1, dim = c(numberOfColumns, numberOfRows, numberOfSlices ) )

    # index values listed as characters: creates empty array of DICOM orientation matrices
    index <- as.character(sort(as.numeric(SOPClassUIDList[  which(SOPClassUIDList[,"SeriesInstanceUID"]==SeriesInstanceUID  ),  "ImageOrder" ])))

    # create and fill the vector DOM, of DICOM orientation matrices
    DOM<-c();  nn<-0
    for (n in index) {
      nn<-nn+1
      tmpSOPIUID <- SOPClassUIDList[  which(SOPClassUIDList[,"ImageOrder"]==n  & SOPClassUIDList[,"SeriesInstanceUID"]==SeriesInstanceUID ),  "SOPInstanceUID" ]
      image.arr[,,nn] <- dataStorage$img[[SeriesInstanceUID]][[tmpSOPIUID]]
    }

    resampled <- FALSE
    if( length(new.pixelSpacing) > 0) {
      if( sum(new.pixelSpacing == old.ps) != 3 ) {

        image.arr<-objService$new.trilinearInterpolator(voxelCube = image.arr,
                                                        pixelSpacing.new = new.pixelSpacing,
                                                        pixelSpacing.old = old.ps)

        new.pixelSpacing <- c( numberOfColumns/dim(image.arr)[1] * old.ps[1],
                               numberOfRows/dim(image.arr)[2] * old.ps[2],
                               new.pixelSpacing[3])

        new.numberOfColumns <- dim(image.arr)[1]
        new.numberOfRows <- dim(image.arr)[2]
        resampled <- TRUE
      }
    }

    for (n in index) {
      tmpSOPIUID <- SOPClassUIDList[  which(SOPClassUIDList[,"ImageOrder"]==n  & SOPClassUIDList[,"SeriesInstanceUID"]==SeriesInstanceUID ),  "SOPInstanceUID" ]
      Dic.Or.Mat <- dataStorage$info[[SeriesInstanceUID]][[tmpSOPIUID]]$orientationMatrix[c(1:3,5:7,13:15)]
      if(all(new.pixelSpacing==old.ps)==FALSE) {
        Dic.Or.Mat[1:3] <- Dic.Or.Mat[1:3]* (new.pixelSpacing[1]/old.ps[1])
        Dic.Or.Mat[4:6] <- Dic.Or.Mat[4:6]* (new.pixelSpacing[2]/old.ps[2])
      }
      DOM<-c(DOM, Dic.Or.Mat)
    }

    # fills the vectors of X and Y coordinates
    # and other Vectors 'associatedInstanceNumberVect' and 'arrayInstanceNumberWithROI'
    TotalX<- -10000;  TotalY<- -10000;  arrayAssociationROIandSlice<- -10000;
    OriginX<- -10000;   OriginY<- -10000
    associatedInstanceNumberVect<- -10000
    contatoreROI<-1; indiceDOM<-1;
    referencedFORUID <- dataStorage$info$structures[[Structure]]$FrameOfReferenceUID

    # for each instance number
    for (n in index) {

      tmpSOPIUID <- SOPClassUIDList[  which(SOPClassUIDList[,"ImageOrder"]==n & SOPClassUIDList[,"SeriesInstanceUID"]==SeriesInstanceUID  ),  "SOPInstanceUID" ]
      numeroROIComplanari <- which(dataStorage$info$structures[[Structure]]$associatedSlices[,"ReferencedSOPInstanceUID"] == tmpSOPIUID)
      if( length( numeroROIComplanari ) > 0 ) {

        elementoListROI <- dataStorage$structures[[Structure]][[tmpSOPIUID]]

        for( xyz in 1:length(elementoListROI)) {
          field2Order <- SOPClassUIDList[  which(SOPClassUIDList[,"ImageOrder"]==n  & SOPClassUIDList[,"SeriesInstanceUID"]==SeriesInstanceUID ),  "field2Order" ]

          if(field2Order == "IPP.y") {
            TotalX<-c(TotalX, elementoListROI[[xyz]][,3])
            TotalY<-c(TotalY, elementoListROI[[xyz]][,1])
          }
          if(field2Order == "IPP.z") {
            TotalX<-c(TotalX, elementoListROI[[xyz]][,1])
            TotalY<-c(TotalY, elementoListROI[[xyz]][,2])
          }
# browser()
          numeroPunti <- length( elementoListROI[[xyz]][,1])
          # for each point write which is the related Slice in the cube-matrix
          arrayAssociationROIandSlice<-c(arrayAssociationROIandSlice,rep(   which( index == n ) -1  , numeroPunti ))

          # Usa OriginX and OriginY as terminator
          TotalX <- c( TotalX, OriginX )
          TotalY <- c( TotalY, OriginY )
          arrayAssociationROIandSlice <- c( arrayAssociationROIandSlice, OriginX )

          contatoreROI <- contatoreROI + 1
        }
      }
    }
    # browser()
    quanti.attesi <- length(unique(global_tableROIPointList[ global_tableROIPointList[,"ROIName"]==Structure , "ROIPointList"] ))
    quanti.assegnati <- contatoreROI - 1
    if( quanti.attesi !=  quanti.assegnati)  {
      logObj$sendLog( paste(  c( "only ",quanti.assegnati," polyline have been associated (expected: ",quanti.attesi,")" )   , collapse = '')  );
    }

    # ok, call the Wrapper!
    final.array<-NewMultiPointInPolyObl(
      # array of DICOM Orientation Matrices
      DICOMOrientationVector = DOM,
      # X and Y vector Points
      totalX = TotalX, totalY = TotalY,
      # association between ROIs and Slices in the 3D Matrix
      arrayAssociationROIandSlice = arrayAssociationROIandSlice,
      # matrices dimensions (rows and columns)
      nX = new.numberOfColumns,
      nY = new.numberOfRows,
      nZ = numberOfSlices
    )

    final.array<-array(data = final.array, dim = c(   new.numberOfColumns, new.numberOfRows, numberOfSlices )   )

    # ROTATE THE MATRIX
    for ( i in seq(1,dim(image.arr)[3] )) {
      final.array[,,i]<-t(objService$rotateMatrix(final.array[,,i],rotations=3))
    }

    immagineMascherata <- array(NA, dim=c(  dim(image.arr)[1],dim(image.arr)[2],dim(image.arr)[3]  ))
    immagineMascherata[ which( final.array == 1, arr.ind = TRUE ) ] <- image.arr[ which( final.array == 1, arr.ind = TRUE) ]
    return(list(
      "DOM"=array(DOM, dim = c(3,3,length(index))),
      "final.array"=final.array,
      "masked.images"=immagineMascherata,
      "geometricalInformationOfImages"=list(),
      "resamplingInformation"=list(
        "px" = new.pixelSpacing[1], "py" = new.pixelSpacing[1], "pz" = new.pixelSpacing[3], "resampled" = resampled,
        "numberOfRows" = new.numberOfRows, "numberOfColumns" = new.numberOfColumns)
    )
    )
  }
  old.getROIVoxelsFromCTRMN<-function( Structure = Structure, SeriesInstanceUID = SeriesInstanceUID,
                                         new.pixelSpacing=c() ) {
    objService<-services()

    # define some variables to make more clear the code
    numberOfRows<-as.numeric(dataStorage$info[[SeriesInstanceUID]][[1]]$Rows);
    numberOfColumns<-as.numeric(dataStorage$info[[SeriesInstanceUID]][[1]]$Columns);
    numberOfSlices<-length(dataStorage$img[[SeriesInstanceUID]]);

    browser()

    old.ps <- getPixelSpacing( SeriesInstanceUID =  SeriesInstanceUID)
    if( length(new.pixelSpacing)>0 ) {
      if(new.pixelSpacing[1] == old.ps[1] & new.pixelSpacing[2]==old.ps[2]) new.pixelSpacing<-c()
    }
    if(  length(new.pixelSpacing)>0 ){
      new.pixelSpacing <- c(new.pixelSpacing,old.ps[3])
      ratio.ps.x <- old.ps[1] / new.pixelSpacing[1];   ratio.ps.y <- old.ps[2] / new.pixelSpacing[2]
    } else {
      new.numberOfRows <- numberOfRows;   new.numberOfColumns <- numberOfColumns
    }

    # initialize the image array with the right dimension
    image.arr<-array( data = -1, dim = c(numberOfColumns, numberOfRows, numberOfSlices ) )

    # index values listed as characters: creates empty array of DICOM orientation matrices
    index <- as.character(sort(as.numeric(SOPClassUIDList[  which(SOPClassUIDList[,"SeriesInstanceUID"]==SeriesInstanceUID  ),  "ImageOrder" ])))

    # create and fill the vector DOM, of DICOM orientation matrices
    DOM<-c();  nn<-0
    for (n in index) {
      nn<-nn+1
      tmpSOPIUID <- SOPClassUIDList[  which(SOPClassUIDList[,"ImageOrder"]==n  & SOPClassUIDList[,"SeriesInstanceUID"]==SeriesInstanceUID ),  "SOPInstanceUID" ]
      image.arr[,,nn] <- dataStorage$img[[SeriesInstanceUID]][[tmpSOPIUID]]
    }

    resampled <- FALSE
    if( length(new.pixelSpacing) > 0) {
      if( sum(new.pixelSpacing == old.ps) != 3 ) {

        image.arr<-objService$new.trilinearInterpolator(voxelCube = image.arr,
                                                        pixelSpacing.new = new.pixelSpacing,
                                                        pixelSpacing.old = old.ps)

        new.pixelSpacing <- c( numberOfColumns/dim(image.arr)[1] * old.ps[1],
                               numberOfRows/dim(image.arr)[2] * old.ps[2],
                               new.pixelSpacing[3])

        new.numberOfColumns <- dim(image.arr)[1]
        new.numberOfRows <- dim(image.arr)[2]
        resampled <- TRUE
      }
    }

    for (n in index) {
      tmpSOPIUID <- SOPClassUIDList[  which(SOPClassUIDList[,"ImageOrder"]==n  & SOPClassUIDList[,"SeriesInstanceUID"]==SeriesInstanceUID ),  "SOPInstanceUID" ]
      Dic.Or.Mat <- dataStorage$info[[SeriesInstanceUID]][[tmpSOPIUID]]$orientationMatrix[c(1:3,5:7,13:15)]
      if(all(new.pixelSpacing==old.ps)==FALSE) {
        Dic.Or.Mat[1:3] <- Dic.Or.Mat[1:3]* (new.pixelSpacing[1]/old.ps[1])
        Dic.Or.Mat[4:6] <- Dic.Or.Mat[4:6]* (new.pixelSpacing[2]/old.ps[2])
      }
      DOM<-c(DOM, Dic.Or.Mat)
    }

    # fills the vectors of X and Y coordinates
    # and other Vectors 'associatedInstanceNumberVect' and 'arrayInstanceNumberWithROI'
    TotalX<- -10000;  TotalY<- -10000;  arrayAssociationROIandSlice<- -10000;
    OriginX<- -10000;   OriginY<- -10000
    associatedInstanceNumberVect<- -10000
    contatoreROI<-1; indiceDOM<-1;
    referencedFORUID <- dataStorage$info$structures[[Structure]]$FrameOfReferenceUID

    # for each instance number
    for (n in index) {

      tmpSOPIUID <- SOPClassUIDList[  which(SOPClassUIDList[,"ImageOrder"]==n & SOPClassUIDList[,"SeriesInstanceUID"]==SeriesInstanceUID  ),  "SOPInstanceUID" ]
      numeroROIComplanari <- which(dataStorage$info$structures[[Structure]]$associatedSlices[,"ReferencedSOPInstanceUID"] == tmpSOPIUID)
      if( length( numeroROIComplanari ) > 0 ) {

        elementoListROI <- dataStorage$structures[[Structure]][[tmpSOPIUID]]

        for( xyz in 1:length(elementoListROI)) {
          TotalX<-c(TotalX, elementoListROI[[xyz]][,1])
          TotalY<-c(TotalY, elementoListROI[[xyz]][,2])

          numeroPunti <- length( elementoListROI[[xyz]][,1])
          # for each point write which is the related Slice in the cube-matrix
          arrayAssociationROIandSlice<-c(arrayAssociationROIandSlice,rep(   which( index == n ) -1  , numeroPunti ))

          # Usa OriginX and OriginY as terminator
          TotalX <- c( TotalX, OriginX )
          TotalY <- c( TotalY, OriginY )
          arrayAssociationROIandSlice <- c( arrayAssociationROIandSlice, OriginX )

          contatoreROI <- contatoreROI + 1
        }
      }
    }

    quanti.attesi <- length(unique(global_tableROIPointList[ global_tableROIPointList[,"ROIName"]==Structure , "ROIPointList"] ))
    quanti.assegnati <- contatoreROI - 1
    if( quanti.attesi !=  quanti.assegnati)  {
      logObj$sendLog( paste(  c( "only ",quanti.assegnati," polyline have been associated (expected: ",quanti.attesi,")" )   , collapse = '')  );
    }

    # ok, call the Wrapper!
    final.array<-NewMultiPointInPolyObl(
      # array of DICOM Orientation Matrices
      DICOMOrientationVector = DOM,
      # X and Y vector Points
      totalX = TotalX, totalY = TotalY,
      # association between ROIs and Slices in the 3D Matrix
      arrayAssociationROIandSlice = arrayAssociationROIandSlice,
      # matrices dimensions (rows and columns)
      nX = new.numberOfColumns,
      nY = new.numberOfRows,
      nZ = numberOfSlices
    )

    final.array<-array(data = final.array, dim = c(   new.numberOfColumns, new.numberOfRows, numberOfSlices )   )

    # ROTATE THE MATRIX
    for ( i in seq(1,dim(image.arr)[3] )) {
      final.array[,,i]<-t(objService$rotateMatrix(final.array[,,i],rotations=3))
    }

    immagineMascherata <- array(NA, dim=c(  dim(image.arr)[1],dim(image.arr)[2],dim(image.arr)[3]  ))
    immagineMascherata[ which( final.array == 1, arr.ind = TRUE ) ] <- image.arr[ which( final.array == 1, arr.ind = TRUE) ]
    return(list(
      "DOM"=array(DOM, dim = c(3,3,length(index))),
      "final.array"=final.array,
      "masked.images"=immagineMascherata,
      "geometricalInformationOfImages"=list(),
      "resamplingInformation"=list(
        "px" = new.pixelSpacing[1], "py" = new.pixelSpacing[1], "pz" = new.pixelSpacing[3], "resampled" = resampled,
        "numberOfRows" = new.numberOfRows, "numberOfColumns" = new.numberOfColumns)
    )
    )
  }
  #=================================================================================
  # NewMultiPointInPolyObl
  #=================================================================================
  NewMultiPointInPolyObl<-function(DICOMOrientationVector,totalX,totalY,arrayAssociationROIandSlice,nX,nY,nZ ) {

    maxX<-max(totalX)
    minX<-min(totalX[which(totalX>-10000)])
    maxY<-max(totalY)
    minY<-min(totalY[which(totalY>-10000)])

    # creates the PIPvector
    PIPvector<-rep.int(x = 0, times = nX * nY * nZ)
    numberOfPoints<-length(totalX);
    result<-.C("NewMultiPIPObl",
               as.integer(PIPvector), as.double(totalX), as.double(totalY), as.integer(numberOfPoints),
               as.integer(nX), as.integer(nY), as.integer(nZ),
               as.integer(arrayAssociationROIandSlice),
               as.double(DICOMOrientationVector),as.double(minX),as.double(maxX),as.double(minY),as.double(maxY))
    return(result[[1]])
  }
  #=================================================================================
  # getPixelSpacing
  #=================================================================================
  getPixelSpacing <- function( SeriesInstanceUID  = NA) {
    if(is.na(SeriesInstanceUID)) {
      SeriesInstanceUID <- giveBackImageSeriesInstanceUID()
      if( length(SeriesInstanceUID) > 1) stop("Too many SeriesIntanceUID have been found: which one?")
    }
    riga <- which(SOPClassUIDList[ ,"SeriesInstanceUID"]==SeriesInstanceUID & SOPClassUIDList[ ,"type"]=="IMG")[1]
    p.x <- SOPClassUIDList[ riga, "p.x"]
    p.y <- SOPClassUIDList[ riga, "p.y"]
    p.z <- SOPClassUIDList[ riga, "p.z"]
    return( as.numeric(c( p.x , p.y, p.z ) ))
  }
  #=================================================================================
  # giveBackImageSeriesInstanceUID
  #=================================================================================
  giveBackImageSeriesInstanceUID<-function() {
    arr.SeriesInstanceUID <- unique(SOPClassUIDList[  which(SOPClassUIDList[ ,"type"]=="IMG")    ,"SeriesInstanceUID"])
    # if( length(arr.SeriesInstanceUID) > 1 ) stop("Error, two SeriesInstanceUID seems to be present")
    if( length(arr.SeriesInstanceUID) ==0 ) stop("Error, no images seems to be loaded")
    # return(arr.SeriesInstanceUID[1])
    return( arr.SeriesInstanceUID )
  }
  #=================================================================================
  # createImageVoxelCube
  # create the imageVoxelCube for the current obj and for the image stored
  #=================================================================================
  createImageVoxelCube<-function( SeriesInstanceUID = NA ) {
    # if not passed, get the series Instance UID of the images
    if(is.na(SeriesInstanceUID)) {
      SeriesInstanceUID <- giveBackImageSeriesInstanceUID()
      if( length(SeriesInstanceUID) > 1) stop("Too many SeriesIntanceUID have been found: which one?")
    }

    Rows <- dim(dataStorage$img[[SeriesInstanceUID]][[1]])[1]
    Columns <- dim(dataStorage$img[[SeriesInstanceUID]][[1]])[2]
    Slices <- length(which(SOPClassUIDList[ ,"SeriesInstanceUID"]==SeriesInstanceUID))

    cubone<-array(data = 0,dim = c(Columns,Rows,Slices))

    for( i in seq( 1 , Slices ) )  {
      SOPInstanceUID <- SOPClassUIDList[ which(SOPClassUIDList[ ,"SeriesInstanceUID"]==SeriesInstanceUID & SOPClassUIDList[ ,"ImageOrder"]==i), "SOPInstanceUID"]
      cubone[,,i] <- dataStorage$img[[SeriesInstanceUID]][[SOPInstanceUID]]
    }
    return(cubone)
  }
  #=================================================================================
  # getDICOMTag
  # tag = which tag
  # fileName = nome del file (se presente)
  #=================================================================================
  getDICOMTag<-function( tag = tag, fileName ) {
    obj.S<-services();
    if(tag == "7fe0,0010") return( getImageFromRAW(fileName)  );
    return(obj.S$getDICOMTag(tag = tag,fileName = fileName, folderCleanUp = internalAttributes$attr_folderCleanUp ))
  }
  #=================================================================================
  # getROIImageImageAssociations
  #=================================================================================
  getROIImageImageAssociations<-function( ROIName ) {
    if( length( ROIName ) != 1 ) stop("Please, specify the interested ROI")
    kkk <- global_tableROIPointList[  which( global_tableROIPointList[,"ROIName"] %in% ROIName ) ,"ReferencedSOPInstanceUID"]
    aaa <- table(SOPClassUIDList[  which( SOPClassUIDList[,"SOPInstanceUID"] %in% kkk ) , "kind"])
    return(aaa)
  }
  #=================================================================================
  # getFilesInfo
  #=================================================================================
  getFilesInfo<-function( ) {
    return( SOPClassUIDList )
  }
  #=================================================================================
  # get.CT.SeriesInstanceUID
  #=================================================================================
  get.CT.SeriesInstanceUID<-function( ) {
    SIUID <- unique(SOPClassUIDList[  which(SOPClassUIDList[,"kind"] == "CTImageStorage"), "SeriesInstanceUID" ])
    return( SIUID )
  }
  #=================================================================================
  # get.MR.SeriesInstanceUID
  #=================================================================================
  get.MRI.SeriesInstanceUID<-function( ) {
    SIUID <- unique(SOPClassUIDList[  which(SOPClassUIDList[,"kind"] == "MRImageStorage"), "SeriesInstanceUID" ])
    return( SIUID )
  }
  #=================================================================================
  # get.MR.SeriesInstanceUID
  #=================================================================================
  get.PET.SeriesInstanceUID<-function( ) {
    SIUID <- unique(SOPClassUIDList[  which(SOPClassUIDList[,"kind"] == "PositronEmissionTomographyImageStorage"), "SeriesInstanceUID" ])
    return( SIUID )
  }

  #=================================================================================
  # class
  #=================================================================================
  class<-function( what="class") {
    if(what=="class") return("geoLet");
    if(what=="version") return("3");
    return("")
  }
  #=================================================================================
  # setAttribute
  #=================================================================================
  setAttribute<-function( folderCleanUp = NA, verbose = NA, RTStructFileName = NA, maxROIPlaneDistance = NA ) {
    if(!is.na(folderCleanUp)) internalAttributes$attr_folderCleanUp <<- folderCleanUp
    if(!is.na(verbose)) internalAttributes$verbose <<- verbose
    if(!is.na(RTStructFileName)) internalAttributes$explicitRTStructFileName <<- RTStructFileName
    if(!is.na(maxROIPlaneDistance))  internalAttributes$maxDistanceForImageROICoupling <<- maxROIPlaneDistance
  }
  #=================================================================================
  # Constructor
  #=================================================================================
  constructor<-function( ) {

    # Attributes - set by user
    internalAttributes$maxDistanceForImageROICoupling<<-0.2
    internalAttributes$explicitRTStructFileName<<-NA
    internalAttributes$verbose<<-TRUE
    internalAttributes$attr_folderCleanUp<<-FALSE                      # force to re-dump DICOM files
    internalAttributes$attr_ROIVoxelMemoryCache<<-TRUE                 # force to cache ROI Voxel
    internalAttributes$attr_ROIVoxelMemoryCacheArray<<-list();
    internalAttributes$attr_arrayXMLCache<<-list();                    # array containint XML files
    internalAttributes$attr_dataChache<<-list();
    internalAttributes$attr_attributeList<<-''
    internalAttributes$attr_loadXMLInCache<<-FALSE
    internalAttributes$attr_loadRAWInCache<<-TRUE
    internalAttributes$attr_arrayXMLCache<<-list();
    internalAttributes$attr_arrayRAWCache<<-list();
    internalAttributes$attr_ROI.non.compl<<-c();
    internalAttributes$attr_neglectedROIs<<-c()
    internalAttributes$attr_withRTStruct<<-TRUE
    internalAttributes$attr_mainFrameOfReferenceUID<<-NA              # frameOfReference (geometry)
    internalAttributes$defaultExtension.dicom <<- ""
    internalAttributes$defaultExtension.nifti<<- ".nii.gz"
    internalAttributes$threshold.4.niftiROI<<- 0.4

    # Internal Structures and objs
    logObj <<- logHandler()                                   # log/error handler Object
    dataStorage <<- list()                                    # memory data structure
    dataStorage$info <<- list()
    SOPClassUIDList <<- list()
    global_tableROIPointList<<-c()

  }
  constructor( )
  return( list(
    "openDICOMFolder"=openDICOMFolder,
    "setAttribute"=setAttribute,
    "getImageVoxelCube"=getImageVoxelCube,
    "getPixelSpacing"=getPixelSpacing,
    "getROIList"=getROIList,
    "getFilesInfo"=getFilesInfo,
    "getROIVoxels"=getROIVoxels,
    "getROIImageImageAssociations"=getROIImageImageAssociations,
    "get.CT.SeriesInstanceUID"=get.CT.SeriesInstanceUID,
    "get.MRI.SeriesInstanceUID"=get.MRI.SeriesInstanceUID,
    "get.PET.SeriesInstanceUID"=get.PET.SeriesInstanceUID
    ))
}
