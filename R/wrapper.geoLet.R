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