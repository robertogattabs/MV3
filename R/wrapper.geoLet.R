#' open a DICOM Folder
#'
#' @description  Open a folder and load the DICOM and NIFTI content. The folder can contain:
#' \itemize{
#' \item one or more image series (CT/MRI/PET)
#' \item one or more DICOM RTStruct file
#' \item or more NIFTI file 
#' }
#' Different series must have different SeriesInstanceUIDs. Geometries are aligned on the base of the FrameOfReferenceUID. 
#' 
#' @param obj.geoLet a \code{geoLet} obj
#' @param pathToOpen the path of the folder containig the objects to load. This folder have to not to contain subfolders
#' @export
#' @examples \dontrun{
#' 
#' # instantiate two geoLet() objects
#' library(moddicomV3)
#' ooo <- geoLet()
#' aaa <- geoLet()
#' 
#' # load the series
#' geo.openDICOMFolder(obj.geoLet = aaa,pathToOpen =  "/home/localadmin/sharedFolder/GEDiscovery690/")
#' geo.openDICOMFolder(obj.geoLet = ooo,pathToOpen = "/projects/moddicom/MV3/img.testing.unit/test.01/axial.01.MRI/")
#' 
#' }
geo.openDICOMFolder <- function( obj.geoLet, pathToOpen ) {
  obj.geoLet$openDICOMFolder( pathToOpen = pathToOpen)
}


#' get the Image Voxe lCube
#'
#' @description  Return the Image Voxel Cube (IVC) of a loaded DICOM Series. If the object contains only one serie, the parameter \code{SeriesInstanceUID} is
#' optional, otherwise is necessary to specify which Series should be retrieved.
#' @param obj.geoLet a \code{geoLet} obj
#' @param ps.x (optional) a new pixel spacing in case you want back an interpolated Voxel Cube
#' @param ps.y (optional) idem
#' @param ps.z (optional) idem
#' @param SeriesInstanceUID the SeriesIntanceUID of the Images you want back. Optional, if you have loaded only one serie.
#' @export
#' @examples \dontrun{
#' 
#' # instantiate two geoLet() objects
#' 
#' library(moddicomV3)
#' ooo <- geoLet()
#' 
#' # load the series
#' 
#' geo.openDICOMFolder(obj.geoLet = ooo,pathToOpen =  "/home/localadmin/sharedFolder/GEDiscovery690/")
#' 
#' # get back the voxelCube of the CT
#' 
#' VC <- geo.getImageVoxelCube(ooo, SeriesInstanceUID = "1.2.840.113619.2.290.3.2831166708.4.1454914198.540.3" )
#' 
#' # If I don't know the SeriesInstanceUID (it's quite cryptic, actually) but I know I want back the CT, I can use :
#' 
#' SIUID.ct <- geo.get.CT.SeriesInstanceUID( ooo )
#' VC <- geo.getImageVoxelCube(ooo, SeriesInstanceUID = SIUID.ct )
#' }
geo.getImageVoxelCube <- function( obj.geoLet, ps.x=NA, ps.y=NA, ps.z=NA , SeriesInstanceUID = NA  ) {
  return( obj.geoLet$getImageVoxelCube( ps.x = ps.x, ps.y = ps.y, ps.z = ps.z , SeriesInstanceUID = SeriesInstanceUID ) )
}


#' return the stored SeriesIntanceUIDs for CTs
#'
#' @description  Return SeriesInstanceUIDs for the CT images loaded by a previous \code{openDICOMFolder}. cause
#' 
#' @param obj.geoLet a \code{geoLet} obj
#' @export
#' @examples \dontrun{
#' 
#' # instantiate two geoLet() objects
#' 
#' library(moddicomV3)
#' ooo <- geoLet()
#' 
#' # load the series
#' 
#' geo.openDICOMFolder(obj.geoLet = ooo,pathToOpen =  "/home/localadmin/sharedFolder/GEDiscovery690/")
#' 
#' # show me the SeriesIntanceUIDs of the stored CTs. (I could have more than one CT series)
#' 
#' SIUID.ct <- geo.get.CT.SeriesInstanceUID( ooo )
#' }
geo.get.CT.SeriesInstanceUID <- function( obj.geoLet ) {
  return( obj.geoLet$geo.get.CT.SeriesInstanceUID() )
}

#' return the stored SeriesIntanceUIDs for MRIs
#'
#' @description  Return SeriesInstanceUIDs for the MRI images loaded by a previous \code{openDICOMFolder}. For examples, see \code{geo.get.CT.SeriesInstanceUID}
#' 
#' @param obj.geoLet a \code{geoLet} obj
#' @export
geo.get.MRI.SeriesInstanceUID <- function( obj.geoLet ) {
  return( obj.geoLet$geo.get.MRI.SeriesInstanceUID() )
}

#' return the stored SeriesIntanceUIDs for PET
#'
#' @description  Return SeriesInstanceUIDs for the MRI images loaded by a previous \code{openDICOMFolder}. For examples, see \code{geo.get.CT.SeriesInstanceUID}
#' 
#' @param obj.geoLet a \code{geoLet} obj
#' @export
geo.get.PET.SeriesInstanceUID <- function( obj.geoLet ) {
  return( obj.geoLet$geo.get.PET.SeriesInstanceUID() )
}

#' get the Pixel Spacing
#'
#' @description  Retrieve and return the Pixel Spacing from on of the images (random) loaded by by a previous \code{openDICOMFolder}. If more series have been loaded, the SeriesInstanceUIDs  must be specified. The pixel spacing is retrieved from a random image (I actually hope you have all the images with the same pixel spacing).
#' @param obj.geoLet a \code{geoLet} obj
#' @param SeriesInstanceUID the SeriesIntanceUID of the Images you want back. Optional, if you have loaded only one serie.
#' @export
#' @examples \dontrun{
#' 
#' # instantiate two geoLet() objects
#' 
#' library(moddicomV3)
#' ooo <- geoLet()
#' 
#' # load the series
#' 
#' geo.openDICOMFolder(obj.geoLet = ooo,pathToOpen =  "/home/localadmin/sharedFolder/GEDiscovery690/")
#' 
#' # get the pixelspacingback the voxelCube of the first of the two loaded MRI series
#' 
#' geo.getPixelSpacing(ooo, SeriesInstanceUID = geo.get.MRI.SeriesInstanceUID( ooo )[1] )
#' }
geo.getPixelSpacing <- function( obj.geoLet, SeriesInstanceUID = NA ) {
  return( obj.geoLet$getPixelSpacing( SeriesInstanceUID = SeriesInstanceUID ) )
}


#' list the loaded ROIs
#'
#' @description  return a list of the loades ROIs. ROIs can be loaded from DICOM-RT struct files or NIFTI file via a \code{openDICOMFolder} method of function. The ROIs, coming from DICOM or NIFTI are presented in the same list and the NIFTI have the extension \code{.nii}
#' @param obj.geoLet a \code{geoLet} obj
#' @export
#' @examples \dontrun{
#' 
#' # instantiate two geoLet() objects
#' 
#' library(moddicomV3)
#' ooo <- geoLet()
#' 
#' # load the series
#' 
#' geo.openDICOMFolder(obj.geoLet = ooo,pathToOpen =  "/home/localadmin/sharedFolder/GEDiscovery690/")
#' 
#' # get the pixelspacingback the voxelCube of the first of the two loaded MRI series
#' 
#' geo.getROIList(ooo)
#' }
geo.getROIList <- function( obj.geoLet  ) {
  return( obj.geoLet$getROIList(  ) )
}

#' get information on the loaded files
#'
#' @description  return a table containing information on the loaded files (images, structures, etc..). Not all the columns can be filled for all the rows, it depend on the kind of object the row refers.
#' The meaning of the columns is the following:
#' \itemize{
#' \item \code{fileName} the fullpath name of the loaded file
#' \item \code{tag} the SOPClassUID
#' \item \code{kind} the SOPClassUID in a human readeable format (if recognized)
#' \item \code{type} set to \code{IMG} if the file contains an image, \code{RTStruct} for a DICOMRT struct, \code{NIFTI} for a NIFTI file.
#' \item \code{IPP.x} the value in the first position in the Image Position Patient (0020,0032)  array tag 
#' \item \code{IPP.y} the value in the second position in the Image Position Patient (0020,0032)  array tag 
#' \item \code{IPP.z} the value in the third position in the Image Position Patient (0020,0032)  array tag
#' \item \code{FrameOfReferenceUID} the FrameOfReferenceUID indicated in the file.
#' \item \code{ImageOrder} the position of the image in the series. This is calculated based on the IPP.x, IPP.y or IPP.z depending on the Image Orientation Patient
#' \item \code{field2Order} can be \code{IPP.z} for axial or para-axial series, \code{IPP.y} for sagittal or para-sagittal series or \code{IPP.x} for coronal or para-coronal images
#' \item \code{p.x} the first value of the Pixel Spacing (0028,0030) tag. This not necessarly means the x direction in the 3D geometry: it means the 'x' direction (column) in the voxel cube
#' \item \code{p.y} the second value of the Pixel Spacing (0028,0030) tag. This not necessarly means the y direction in the 3D geometry: it means the 'y' direction (rows) in the voxel cube
#' \item \code{p.z} the value of the slice tickness. This not necessarly means the z direction in the 3D geometry: it means the 'z' direction (slice number) in the voxel cube
#' \item \code{SOPInstanceUID} the SOPInstanceUID of the DICOM object
#' \item \code{ImageOrientationPatient} the ImageOrientationPatient as presented in the tag (0020,0037)
#' \item \code{SeriesInstanceUID} the SeriesInstanceUID of the DICOM object
#' }
#' @param obj.geoLet a \code{geoLet} obj
#' 
#' @export
#' @examples \dontrun{
#' 
#' # instantiate two geoLet() objects
#' 
#' library(moddicomV3)
#' ooo <- geoLet()
#' 
#' # load the series
#' 
#' geo.openDICOMFolder(obj.geoLet = ooo,pathToOpen =  "/home/localadmin/sharedFolder/GEDiscovery690/")
#' 
#' # retrieve image infos about the loaded objects
#' 
#' geo.getFilesInfo(ooo)
#' }
geo.getFilesInfo <- function( obj.geoLet  ) {
  return( obj.geoLet$getFilesInfo() ) 
}

#' get image/ROIs association infos
#'
#' @description  return a table containing the number of ROIs associated for each loaded Series. The \code{associaation} between ROIs and images can happen when a polyline of the ROI is coplanar with the image. For point is very important becase if not considered can let you miss a significant number of polyline of the ROI
#' @param obj.geoLet a \code{geoLet} obj
#' @param ROIName the ROIName you want to check
#' @param details (boolean, default to \code{FALSE} ). If set to \code{TRUE} it takes more time but also give back a set of 0 and 1 to say which slices in the imagexovelcube are associated to the ROI
#' @export
#' @examples \dontrun{
#' 
#' # instantiate two geoLet() objects
#' 
#' library(moddicomV3)
#' ooo <- geoLet()
#' 
#' # load the series
#' 
#' geo.openDICOMFolder(obj.geoLet = ooo,pathToOpen =  "/home/localadmin/sharedFolder/GEDiscovery690/")
#' 
#' # get the pixelspacingback the voxelCube of the first of the two loaded MRI series
#' 
#' geo.getROIImageImageAssociations(ooo)
#' }
geo.getROIImageImageAssociations <- function( obj.geoLet, ROIName , details = FALSE ) {
  return( obj.geoLet$getROIImageImageAssociations( ROIName = ROIName, details = details) ) 
}

#' graphs images and ROIs
#'
#' @description  return a table containing the number of ROIs associated for each loaded Series. The \code{associaation} between ROIs and images can happen when a polyline of the ROI is coplanar with the image. For point is very important becase if not considered can let you miss a significant number of polyline of the ROI
#' @param obj.geoLet a \code{geoLet} obj
#' @param arr.ROInames an arrax containing the ROIs you wanto to plot
#' @param rows the number of rows of the final grid (default = 3)
#' @param cols the number of cols of the final grid (default = 3)
#' @param grey.scale set to TRUE for a grey.level palette. FALSE for a red.level palette.
#' @param alpha The alpha value for the ROIs
#' @param SeriesInstanceUID The SeriesIntanceUID of the image series. This is not needed if only one serie is loaded and is mandatory in case of multiple series
#' @export
#' @examples \dontrun{
#' 
#' # instantiate two geoLet() objects
#' 
#' library(moddicomV3)
#' ooo <- geoLet()
#' 
#' # load the series
#' 
#' geo.openDICOMFolder(obj.geoLet = ooo,pathToOpen =  "/home/localadmin/sharedFolder/GEDiscovery690/")
#' 
#' # plot two ROIs (the former from a DICOM, the latter from a NIFTI)
#' 
#' geo.showROIs(objGeoLet = ooo,arr.ROInames = c("GTV","GUESS.nii"))
#' }
geo.showROIs <-  function( objGeoLet, arr.ROInames, rows = 3, cols = 3 , grey.scale = TRUE, alpha = 0.4, SeriesInstanceUID=NA ) {
  
  # carica il voxel Cube e le ROI
  lst.ROI <- list()
  VC <- objGeoLet$getImageVoxelCube(SeriesInstanceUID = SeriesInstanceUID)
  for( ROIName in arr.ROInames ) { lst.ROI[[ ROIName ]] <- objGeoLet$getROIVoxels(Structure = ROIName,croppedCube = F,onlyVoxelCube = T, SeriesInstanceUID = SeriesInstanceUID) }
  
  # prepara la finestra per l'output
  par(mfrow=c(rows,cols))
  par(mar=c(0,0,0,0))
  
  # prendi le z in cui compare almeno una ROI
  slices <- sort(unique(unlist(lapply( lst.ROI, function( cubo ) { unique(which(!is.na(cubo),arr.ind = T)[,3])  } ))))
  
  # plotta ogni slice
  for( z in 1:length( slices ) ) {
    if( grey.scale == TRUE) image( VC[,,slices[z]],col = grey.colors(255,start = 0))
    else image(VC[,,slices[z]])
    
    # e per ogni slice plotta le ROI che vi sono associate
    for( ROIName in arr.ROInames ) {
      # se per quella slice e quella ROI c'e' almeno un pixel != 0, allora plotta
      if( sum(!is.na(lst.ROI[[ ROIName ]][,,slices[z]])) > 0 ) {
        image( lst.ROI[[ ROIName ]][,,slices[z]],add=T, col = rgb(0.8,0,0,alpha = alpha))  
      }
    }
  } 
}

#' retrurn a DICOM tag
#'
#' @description  give back the value of the specified DICOM tag. The returned tag depends on the nature of the tag, so pay attention to the many options (e.g. a PatientId can be retrieved from any images of a series, are all equals, on the other hand, slice tichkness can change withing the serie)
#' @param obj.geoLet a \code{geoLet} obj
#' @param arr.ROInames an arrax containing the ROIs you wanto to plot
#' @param tag hte DICOM in the DICOM format. Most common tags are for example '0010,0010','0010,0020' (PatientId and Patient Name respectively).
#' @param fileName the name of the DICOM file you want to get the value of the tag. If it is not specified, the first file of the serie is considered. In more than one serie is available, the SeriesInstanceUID needs to be specified, in order to retrieve the name of the first file to use.
#' @param SeriesInstanceUID in case the fileName is missing and there are more than one serie, the SeriesInstanceUID is mandator< 
#' @export
#' @examples \dontrun{
#' 
#' # instantiate two geoLet() objects
#' 
#' library(moddicomV3)
#' ooo <- geoLet()
#' 
#' # load the series
#' 
#' geo.openDICOMFolder(obj.geoLet = ooo,pathToOpen =  "/home/localadmin/sharedFolder/GEDiscovery690/")
#' 
#' # plot two ROIs (the former from a DICOM, the latter from a NIFTI)
#' 
#' getTag(objGeoLet = ooo,tag = "0010,0010", fileName="./test.dcm")
#' }
geo.getTag<-function(objGeoLet , tag , fileName=NA, SeriesInstanceUID = NA) {   
  tagValue <- objGeoLet$getTag(tag = tag, fileName = fileName, SeriesInstanceUID = SeriesInstanceUID)
  return( tagValue )
}

