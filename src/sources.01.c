#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include "R.h"
#include <Rinternals.h>
#include <complex.h>
#include <limits.h>
#define min(a, b) (((a) < (b)) ? (a) : (b));

struct pointInSpace{
  double x;
  double y;
  double z;
};
struct DICOM_OrientationMatrix{
  double a11;		// 1st element of ImageOrientationPatient
  double a21;		// 2nd element of ImageOrientationPatient
  double a31;		// 3rd element of ImageOrientationPatient
  double a12;		// 4th element of ImageOrientationPatient
  double a22;		// 5th element of ImageOrientationPatient
  double a32;		// 6th element of ImageOrientationPatient
  double Sx;		// 1st element of ImagePositionPatient
  double Sy;		// 2nd element of ImagePositionPatient
  double Sz;		// 3nd element of ImagePositionPatient
  double XpixelSpacing;	// 1st element of Pixel Spacing
  double YpixelSpacing;	// 2nd element of Pixel Spacing
  struct DICOM_OrientationMatrix *next;
  struct DICOM_OrientationMatrix *prev;
};

struct _c_data{
  // data passed as parameters
  int xNVoxel,yNVoxel,zNVoxel;
  double xDim,yDim,zDim;
  int newXNVoxel,newYNVoxel,newZNVoxel;
  double newXDim,newYDim,newZDim;
  // Voxel number for each vertex of the recognized Cube
  int xInf,xSup,yInf,ySup,zInf,zSup;
};

/*
 * Function to calculate the index position of an array which refers to a 3D matrix
 * x,y,z are the coors of the matrix you are interested in 
 * nx,ny,nz  are the dimensions of the matrix along the 3 axes
 *
 */
int posDecod(int x, int y, int z, int nx, int ny, int nz) {
  return( z*ny*nx+ y*nx + x   );
}

//
// _c_getCubeVertex
//
// Prende le coordinate in floating point del punto da mappare e restituisce l'intero
// corrispondente ai centroidi relativi alle coordinate dei vertici del cubo da interpolare
// La formula per calcolare il Kinf. è   Kinf = int( (2*xPos - dx) / (2 * dx) )
//
void _c_getCubeVertex(struct _c_data * punt, double xPos, double yPos, double zPos, int ct) {

  punt->xInf = (int)((2 * xPos - punt->xDim) / ( 2 * punt->xDim));
  punt->yInf = (int)((2 * yPos - punt->yDim) / ( 2 * punt->yDim));
  punt->zInf = (int)((2 * zPos - punt->zDim) / ( 2 * punt->zDim));
  punt->xSup = punt->xInf + 1;
  punt->ySup = punt->yInf + 1;
  punt->zSup = punt->zInf + 1;

  // se xPos < punt->xDim significa che sono nel bordo. Considera allora il valore del primo valore di X per x0 e x1
  // (in sostanza proietta all'esterno il valore dei pixel). Analogo se maggiore.
  // STO ANDANDO IN OVERRIDE RISPETTO A QUANDO EVENTUALMENTE CALCOLATO PRIMA
  if( (2*xPos-punt->xDim) < punt->xDim/2 ) {punt->xInf = 0; punt->xSup = 0;}
  if( (2*yPos-punt->yDim) < punt->yDim/2 ) {punt->yInf = 0; punt->ySup = 0;}
  if( (2*zPos-punt->zDim) < punt->zDim/2  ) {punt->zInf = 0; punt->zSup = 0;}

  // Faccio il reciproco di xPos per vedere se non sto' toccando l'altro estremo
  xPos = punt->xDim * punt->xNVoxel - xPos;
  yPos = punt->yDim * punt->yNVoxel - yPos;
  zPos = punt->zDim * punt->zNVoxel - zPos;
  if( (2*xPos-punt->xDim) < punt->xDim/2 ) {punt->xInf = punt->xNVoxel - 1; punt->xSup = punt->xNVoxel - 1;}
  if( (2*yPos-punt->yDim) < punt->yDim/2 ) {punt->yInf = punt->yNVoxel - 1; punt->ySup = punt->yNVoxel - 1;}
  if( (2*zPos-punt->zDim) < punt->zDim/2  ) {punt->zInf = punt->zNVoxel - 1; punt->zSup = punt->zNVoxel - 1;}
}

//
// _c_TrilinearInterpolation
//
// Effettua l'interpolazione considerando i valori ai vertici e le dimensioni fisiche del cubo
//
double _c_TrilinearInterpolation(double x0y0z0, double x0y0z1, double x0y1z0, double x0y1z1, double x1y0z0, double x1y0z1, double x1y1z0, double x1y1z1, double x0,double y0, double z0, double dx1x0, double dy1y0, double dz1z0, double x, double y,double z) {
  double xd,yd,zd,c00,c01,c10,c11,c0,c1,c;
  xd = (x-x0)/dx1x0;
  yd = (y-y0)/dy1y0;
  zd = (z-z0)/dz1z0;
  c00 = x0y0z0*(1-xd)+x1y0z0*xd;
  c10 = x0y1z0*(1-xd)+x1y1z0*xd;
  c01 = x0y0z1*(1-xd)+x1y0z1*xd;
  c11 = x0y1z1*(1-xd)+x1y1z1*xd;

  c0 = c00*(1-yd)+c10*yd;
  c1 = c01*(1-yd)+c11*yd;
  c = c0*(1-zd)+c1*zd;
  return c;
}

void newnewtrilinearInterpolator(
    int *NXold, int *NYold,int *NZold,
    int *NXnew, int *NYnew,int *NZnew,
    double *oldXps, double *oldYps, double *oldZps,
    double *newXps, double *newYps, double *newZps,
    double *values,double *returnMatrix) {
  int ct;
  int zVoxelProgressivo,yVoxelProgressivo,xVoxelProgressivo;
  double zPos,yPos,xPos,valoreCalcolato;
  struct _c_data *punt;

  // Alloca un puntatore alla struttura _c_data
  punt = (struct _c_data *)calloc(1,sizeof(struct _c_data));
  if( punt == NULL ) return;

  // calcola il nuovo passo e copia i valori in _c_data: l'idea è quella di ridurre il passaggio parametri
  // fra funzioni ed il clone delle variabili per preservare memoria
  punt->newXDim = *newXps;  punt->newYDim = *newYps;  punt->newZDim = *newZps;
  punt->xNVoxel = *NXold; punt->yNVoxel = *NYold; punt->zNVoxel = *NZold;
  punt->xDim = *oldXps; punt->yDim = *oldYps; punt->zDim = *oldZps;
  int maxNewYVoxel,maxNewZVoxel;

  // Pulisci la matrice di destinazione
  for(int ct=0;ct< ((*NXnew)*(*NYnew)*(*NZnew));ct++ ) returnMatrix[ct]=0;
  ct=0;
  for( zPos = punt->newZDim/2, zVoxelProgressivo=0; zVoxelProgressivo < *NZnew; zPos+=punt->newZDim, zVoxelProgressivo++ ) {
    for( yPos = punt->newYDim/2, yVoxelProgressivo=0; yVoxelProgressivo < *NYnew; yPos+=punt->newYDim, yVoxelProgressivo++ ) {
      for( xPos = punt->newXDim/2, xVoxelProgressivo=0; xVoxelProgressivo < *NXnew; xPos+=punt->newXDim, xVoxelProgressivo++ ) {

        // acquisisci i dati relativi ai voxel della vecchia matrice i cui centroidi sono
        // ai vertici del cubo da interpolare
        _c_getCubeVertex( punt , xPos , yPos , zPos ,ct );

        // reperisci i valori di tali voxel dalla vecchia matrice
        // ed effettua l'interpolazione rispetto a tali centroidi
        valoreCalcolato = _c_TrilinearInterpolation(
          values[punt->zInf*(*NXold)*(*NYold)+punt->yInf*(*NXold)+punt->xInf],  //x0y0z0 (sample value)
          values[punt->zSup*(*NXold)*(*NYold)+punt->yInf*(*NXold)+punt->xInf],  //x0y0z1 (sample value)
          values[punt->zInf*(*NXold)*(*NYold)+punt->ySup*(*NXold)+punt->xInf],  //x0y1z0 (sample value)
          values[punt->zSup*(*NXold)*(*NYold)+punt->ySup*(*NXold)+punt->xInf],  //x0y1z1 (sample value)
          values[punt->zInf*(*NXold)*(*NYold)+punt->yInf*(*NXold)+punt->xSup],  //x1y0z0 (sample value)
          values[punt->zSup*(*NXold)*(*NYold)+punt->yInf*(*NXold)+punt->xSup],  //x1y0z1 (sample value)
          values[punt->zInf*(*NXold)*(*NYold)+punt->ySup*(*NXold)+punt->xSup],  //x1y1z0 (sample value)
          values[punt->zSup*(*NXold)*(*NYold)+punt->ySup*(*NXold)+punt->xSup],  //x1y1z1 (sample value)
          punt->xInf*(*oldXps)+(*oldXps)/2, //x0
          punt->yInf*(*oldYps)+(*oldYps)/2, //y0,
          punt->zInf*(*oldZps)+(*oldZps)/2, //z0,
          *oldXps, //dx1x0,
          *oldYps, //dy1y0,
          *oldZps, //dz1z0,
          xPos, yPos, zPos);

        // memorizza il risultato nella nuova struttura
        returnMatrix[ xVoxelProgressivo + yVoxelProgressivo * (*NXnew) + zVoxelProgressivo * ((*NYnew) * (*NXnew))] = valoreCalcolato;
        // if( valoreCalcolato!=0 ) printf("\nx = %d, y= %d, z= %d",xVoxelProgressivo,yVoxelProgressivo,zVoxelProgressivo);
      }
      if(maxNewYVoxel<yVoxelProgressivo) maxNewYVoxel = yVoxelProgressivo;
    }
    if(maxNewZVoxel<zVoxelProgressivo) maxNewZVoxel = zVoxelProgressivo;
  }
  free(punt);
}


/*
 nvert :   number of vertex
 vertx :   x coords of vertex points
 verty :   y coords of vertex points
 */
int isThePointInsideThePoly(int nvert, double *vertx, double *verty, double testx, double testy, int fromPosition, int toPosition)
{
  int i, j, c = 0;
  for (i = fromPosition, j = fromPosition+nvert-1; i < fromPosition+nvert; j = i++) {
    if ( ((verty[i]>testy) != (verty[j]>testy)) &&
         (testx < (vertx[j]-vertx[i]) * (testy-verty[i]) / (verty[j]-verty[i]) + vertx[i]) )
      c = !c;
  }
  return c;
}

/*
 * Function for calculating the <x,y,z> coords of a voxel into the space
 * starting from rows (Ny) and columns (Nx) of the 2D Matrix of doses
 * Nx	:		columns
 * Ny 	:		rows
 * DOM  :		Dicom Orientation Matrix
 */
struct pointInSpace get3DPosFromNxNy(int Nx, int Ny, struct DICOM_OrientationMatrix DOM) {
  struct pointInSpace pis;
  //printf("\n a11=%lf, Nx=%d Ny=%d, a12=%lf, Sx=%lf  \n",DOM.a11,Nx, Ny,  DOM.a12, DOM.Sx);
  pis.x = DOM.a11 * Nx * DOM.XpixelSpacing + DOM.a12 * Ny * DOM.YpixelSpacing + DOM.Sx;
  pis.y = DOM.a21 * Nx * DOM.XpixelSpacing + DOM.a22 * Ny * DOM.YpixelSpacing + DOM.Sy;
  pis.z = DOM.a31 * Nx * DOM.XpixelSpacing + DOM.a32 * Ny * DOM.YpixelSpacing + DOM.Sz;
  return pis;
}

/*
 * Function for calculating the <Nx,Ny> coords of a voxel into the space
 */
/*
struct pointInSpace getNxNyPos3D( double Px, double Py, double Pz, struct DICOM_OrientationMatrix DOM) {
  struct pointInSpace pis;
  double Xx,  Xy, Xz,  Yx, Yy, Yz, Sx, Sy, Sz;
  Xx = DOM.a11; Xy = DOM.a21; Xz = DOM.a31;
  Yx = DOM.a12; Yy = DOM.a22; Yz = DOM.a32;
  
  pis.x = Xz * (Px-Sx) - Xx*(Pz-Sz) /  (Xz*Yx - Yz*Xx );
  pis.y = Yx * (Py-Sy) - Yy*(Px-Sx) /  (Yx*Xy - Yy*Xx );
  return pis;
}*/


/*
 *  Algorithm for Point In Polygon on Oblique planes
 *  calculated over a multiple series of contours each one on a different slice.
 *	Each slice is ordered according NumSlices value.
 *	Variables to be addressed to the C code are:
 *      PIPvector       vector of Points in Polygon to be arranged in an array in R
 *	totalY: 	vector of all Y coordinates of contours (one appended to the other)
 *	totalX: 	vector of all X coordinates of contours (one appended to the other)
 *      numberOfPoints: number of points declared in totalX/Y
 *      nX, nY, nZ:	x,y,z dimensions of the voxel matrix
 *      arrayAssociationROIandSlice: associates each point in totalX/Y with a slice. terminator (-100000) means the ROI is closed
 *	arrayDOM:		DICOM orientation vector, vectorized DICOM orientation matrix multiplied by pixel spacing
 *
 *      minX,maxX,minY,maxY: min, max values of coords in order to preserve time. In future they will be calculated directly into this C code...
 */
void NewMultiPIPObl (int *PIPvector, double *totalX, double *totalY, int *numberOfPoints,
                     int *nX, int *nY, int *nZ,
                     int *arrayAssociationROIandSlice,
                     double *arrayDOM,
                     double *minX, double *maxX, double *minY, double *maxY) {
  int x, y, z, c, fromPosition, toPosition;
  struct pointInSpace point;
  struct DICOM_OrientationMatrix DOM;

  int baseCursor = 0;
  for (baseCursor = 0; baseCursor < *numberOfPoints; baseCursor++ ) {

    // check if a new ROI is beginning (AND you are not looking the last point!) )
    if ( totalX[ baseCursor ] == -10000 && baseCursor < (*numberOfPoints-1) ) {

      // get the corresponding z position
      z = arrayAssociationROIandSlice[ baseCursor + 1 ];

      // get the DOM
      DOM.a11= arrayDOM[ z * 9 + 0 ]; DOM.a21= arrayDOM[ z * 9 + 1 ];
      DOM.a31= arrayDOM[ z * 9 + 2 ]; DOM.a12= arrayDOM[ z * 9 + 3 ];
      DOM.a22= arrayDOM[ z * 9 + 4 ]; DOM.a32= arrayDOM[ z * 9 + 5 ];
      DOM.Sx=  arrayDOM[ z * 9 + 6 ]; DOM.Sy=  arrayDOM[ z * 9 + 7 ];
      DOM.Sz=  arrayDOM[ z * 9 + 8 ];
      DOM.XpixelSpacing = 1;  DOM.YpixelSpacing = 1;

      // now Run to find the beginning (clear) and the end of the
      // points of the ROI, in the given array
      fromPosition = baseCursor;
      for ( toPosition = fromPosition+1; totalX[toPosition]!=-10000; toPosition++) ;;
      toPosition--;

      // check all the voxel on such z-slice
      // loop through Y axis
      for ( y = 0; y < *nY; y++ ) {
        // loop through X axis
        for ( x = 0; x < *nX; x++ ) {
          point = get3DPosFromNxNy(x, y, DOM);
//          printf("\n (%lf , %lf , %lf ) in [%lf-%lf,%lf-%lf]",point.x,point.y,point.z,*minX,*maxX,*minY,*maxY);
          // check the box, just to avoid redundant computation
          if(point.x < *minX || point.x > *maxX || point.y < *minY || point.y > *maxY) continue;

          // now check if the given point is in the poly
          c = isThePointInsideThePoly(toPosition-fromPosition-1, totalX, totalY,
                                      point.x, point.y,  fromPosition+1,  toPosition+1);

          if( c == 1 )   {
            PIPvector[z * (*nX) * (*nY) + x + y * (*nX)] = !PIPvector[z * (*nX) * (*nY) + x + y * (*nX)];
          }
        }
      }
    }
  }
}

/*
*  Algorithm for Point In Polygon on Oblique planes
*  calculated over a multiple series of contours each one on a different slice.
*	Each slice is ordered according NumSlices value.
*	Variables to be addressed to the C code are:
*      PIPvector       vector of Points in Polygon to be arranged in an array in R
*	     totalY 	       vector of all Y coordinates of contours (one appended to the other)
*	     totalX          vector of all X coordinates of contours (one appended to the other)
*      numberOfPoints  number of points declared in totalX/Y
*      nX, nY, nZ:	   x,y,z dimensions of the voxel matrix
*      arrayAssociationROIandSlice: associates each point in totalX/Y with a slice. terminator (-100000) means the ROI is closed
*	     arrayDOM     		DICOM orientation vector, vectorized DICOM orientation matrix multiplied by pixel spacing
*	     mainOrientation  1 = Z (assiale), 2 = Y (coronale), 3 = X (saggittale)
*
*      minX,maxX,minY,maxY: min, max values of coords in order to preserve time. In future they will be calculated directly into this C code...
*/
void NewMultiOrientedPIPObl (int *PIPvector, double *totalX, double *totalY, double *totalZ,
                             int *numberOfPoints,
                             int *nX, int *nY, int *nZ,
                             int *arrayAssociationROIandSlice,
                             double *arrayDOM,
                             double *minX, double *maxX, double *minY, double *maxY, double *minZ, double *maxZ,
                             int *mainOrientation ) {
  int x, y, z, c, fromPosition, toPosition, sliceNumber;
  struct pointInSpace point;
  struct DICOM_OrientationMatrix DOM;
  int riga, colonna;

  int baseCursor = 0;
  for (baseCursor = 0; baseCursor < *numberOfPoints; baseCursor++ ) {

    // check if a new ROI is beginning (AND you are not looking the last point!) )
    if ( totalX[ baseCursor ] == -10000 && baseCursor < (*numberOfPoints-1) ) {

      // get the number of the slice corresponding to the polyline we are going to analyze
      sliceNumber = arrayAssociationROIandSlice[ baseCursor + 1 ];
      // z = arrayAssociationROIandSlice[ baseCursor + 1 ];

      // get the DOM
      DOM.a11= arrayDOM[ sliceNumber * 9 + 0 ]; DOM.a21= arrayDOM[ sliceNumber * 9 + 1 ];
      DOM.a31= arrayDOM[ sliceNumber * 9 + 2 ]; DOM.a12= arrayDOM[ sliceNumber * 9 + 3 ];
      DOM.a22= arrayDOM[ sliceNumber * 9 + 4 ]; DOM.a32= arrayDOM[ sliceNumber * 9 + 5 ];
      DOM.Sx=  arrayDOM[ sliceNumber * 9 + 6 ]; DOM.Sy=  arrayDOM[ sliceNumber * 9 + 7 ];
      DOM.Sz=  arrayDOM[ sliceNumber * 9 + 8 ];
      // Set the pixelspacing to 1 because they should be already set from the calling code
      DOM.XpixelSpacing = 1;  DOM.YpixelSpacing = 1;

      // now Run to find the beginning (clear) and the end of the
      // points of the ROI, in the given array
      fromPosition = baseCursor;
      for ( toPosition = fromPosition+1; totalX[toPosition]!=-10000; toPosition++) ;;
      toPosition--;

      // check all the voxel on such slice
      for ( riga = 0; riga < *nY; riga++ ) {
        for ( colonna = 0; colonna < *nX; colonna++ ) {
          if( *mainOrientation == 1 ) {  x = colonna; y = riga; }
          if( *mainOrientation == 2 ) {  x = colonna; z = riga; }
          if( *mainOrientation == 3 ) {  y = colonna; z = riga; }
          point = get3DPosFromNxNy(colonna, riga, DOM);
          // check the box, just to avoid redundant computation

          if(point.x < *minX || point.x > *maxX ||
             point.y < *minY || point.y > *maxY ||
             point.z < *minZ || point.z > *maxZ) continue;

          c = isThePointInsideThePoly(toPosition-fromPosition-1, totalX, totalY,
                                      point.x, point.y,  fromPosition+1,  toPosition+1);
          if( c == 1 )   {
            PIPvector[z * (*nX) * (*nY) + x + y * (*nX)] = !PIPvector[z * (*nX) * (*nY) + x + y * (*nX)];
          }
        }
      }
    }
  }
}

/*imported from Nic */

/*
 * Function for calculating the DOUBLE of area of a facet in a triangular mesh
 */
double FacetSurface(double p1X, double p1Y, double p1Z, 
                    double p2X, double p2Y, double p2Z, double p3X, double p3Y, double p3Z) {
  double ax = p2X - p1X;
  double ay = p2Y - p1Y;
  double az = p2Z - p1Z;
  double bx = p3X - p1X;
  double by = p3Y - p1Y;
  double bz = p3Z - p1Z;
  double cx = ay*bz - az*by;
  double cy = az*bx - ax*bz;
  double cz = ax*by - ay*bx;
  //printf("\np2X*p3Y-p3X*p2Y=%lf",p2X * p3Y - p3X * p2Y );
  return sqrt(cx*cx + cy*cy + cz*cz);
}

/*
 * Function for calculating SIX TIMES the signed volume
 * of a tetrahedron in a triangular mesh
 */
double SignedVolumeOfTriangle(double p1X, double p1Y, double p1Z, 
                              double p2X, double p2Y, double p2Z, double p3X, double p3Y, double p3Z) {
  double v321 = p3X*p2Y*p1Z;
  double v231 = p2X*p3Y*p1Z;
  double v312 = p3X*p1Y*p2Z;
  double v132 = p1X*p3Y*p2Z;
  double v213 = p2X*p1Y*p3Z;
  double v123 = p1X*p2Y*p3Z;
  //printf("\ndVolume=%lf", (double)(1.0/6.0));
  return (-v321 + v231 + v312 - v132 - v213 + v123);	
}

/* 
 * Function for calculating the mesh surface
 */
void MeshSurface(double *X, double *Y, double *Z, int *numT, int *V1, int *V2, int *V3, double *Surface) {
  int n;			// counter
  *Surface = 0;	// initial surface
  for (n=0; n<*numT; n++) {
    *Surface = *Surface + FacetSurface(X[V1[n]], Y[V1[n]], Z[V1[n]], X[V2[n]], Y[V2[n]], Z[V2[n]], X[V3[n]], Y[V3[n]], Z[V3[n]]);
  }
  *Surface=0.5 * *Surface;
}
/* 
 * Function for calulating the volume of a mesh
 */
void MeshVolume(double *X, double *Y, double *Z, int *numT, int *V1, int *V2, int *V3, double *Volume) {
  int n;			// counter
  *Volume=0; 		// initial volume
  for (n=0; n<*numT; n++) {
    *Volume = *Volume + SignedVolumeOfTriangle(X[V1[n]], Y[V1[n]], Z[V1[n]], X[V2[n]], Y[V2[n]], Z[V2[n]], X[V3[n]], Y[V3[n]], Z[V3[n]]);
    //printf("\nX = %lf Y = %lf Z = %lf n = %d", X[V1[n]], Y[V1[n]], Z[V1[n]], n);
  }
  *Volume = fabs(*Volume * (double)(1.0/6.0));  // absolute value of a double
}


void c_getInterpolatedSlice2D(
    int *Nx, int *Ny, double *deltaZ, double *zRelativePosition,
    int *NxOut, int *NyOut,
    double *x, double *y,
    double *xOut, double *yOut,
    double *ff,
    double *ft,
    double *res
) {
  
  double deltaX, deltaY;
  double f000,f010,f100,f110,f001,f011,f101,f111;
  double xPos, yPos;
  int Xcurs, Ycurs;
  int intXcurs, intYcurs;
  int lny,lnx;
  double valoreCalcolato;

  deltaX = x[1]-x[0];  deltaY = y[1]-y[0];
  
  for( lny = 0; lny < *NyOut; lny++ ) {    
    for( lnx = 0; lnx < *NxOut; lnx++ ) {      
      xPos = xOut[lnx]; yPos = yOut[lny];
      Xcurs = (xPos - x[0]) / deltaX;  Ycurs = (yPos - y[0]) / deltaY;
      intXcurs = floor(Xcurs); intYcurs = floor(Ycurs);
      // se cade in un range ragionevole
      if( (intXcurs > 0) & ((intXcurs+1) < *Nx) & (intYcurs > 0) & ((intYcurs+1) < *Ny) ) {
        f000 = ff[  intYcurs * (*Nx) + intXcurs ];
        f100 = ff[  intYcurs * (*Nx) + (intXcurs+1) ];
        f010 = ff[  (intYcurs+1) * (*Nx) + intXcurs ];
        f110 = ff[  (intYcurs+1) * (*Nx) + (intXcurs+1) ];
        f001 = ft[  intYcurs * (*Nx) + intXcurs ];
        f101 = ft[  intYcurs * (*Nx) + (intXcurs+1) ];
        f011 = ft[  (intYcurs+1) * (*Nx) + intXcurs ];
        f111 = ft[  (intYcurs+1) * (*Nx) + (intXcurs+1) ];
        
        
        valoreCalcolato = _c_TrilinearInterpolation(
          f000,  //x0y0z0 (sample value)
          f001,  //x0y0z1 (sample value)
          f010,  //x0y1z0 (sample value)
          f011,  //x0y1z1 (sample value)
          f100,  //x1y0z0 (sample value)
          f101,  //x1y0z1 (sample value)
          f110,  //x1y1z0 (sample value)
          f111,  //x1y1z1 (sample value)
          0, //x0
          0, //y0,
          0, //z0,
          deltaX, //dx1x0,
          deltaY, //dy1y0,
          *deltaZ, //dz1z0,
          xPos-x[intXcurs], yPos-y[intYcurs], *zRelativePosition);        

          res[ lny * (*NxOut) + lnx ] = valoreCalcolato;
      }
    }
  }
}



/*
 *  Algorithm for Point In Polygon on a flat 2D plane
 */
void bidimentionalFlatPointInPolygon (
                     int *nX, int *nY,
                     int *numVertex,
                     double *xVertex, 
                     double *yVertex,
                     int *PIPvector 
) {
  int x,y,c;
  for ( y = 0; y < *nY; y++ ) {    // loop through X axis
    for ( x = 0; x < *nX; x++ ) {
      c = isThePointInsideThePoly(*numVertex , xVertex , yVertex, x, y , 0 , 0);
      if( c == 1 )   {
        PIPvector[ x + y * (*nX)] = !PIPvector[ x + y * (*nX)];
      }
    }
  }
}
   /*     
void toDelete_c_getInterpolatedSlice2D(
    int *Nx, int *Ny,
    int *NxOut, int *NyOut,
    double *x, double *y,
    double *xOut, double *yOut,
    double *f,
    double *res
) {
  
  double deltaX, deltaY;
  double f00,f01,f10,f11,festx_0,festx_1;
  double xPos, yPos;
  int Xcurs, Ycurs;
  int intXcurs, intYcurs;
  int lny,lnx;
  
  deltaX = x[1]-x[0];  deltaY = y[1]-y[0];
  
  for( lny = 0; lny < *NyOut; lny++ ) {    
    printf("\n y = %d",lny);
    for( lnx = 0; lnx < *NxOut; lnx++ ) {      
      xPos = xOut[lnx]; yPos = yOut[lny];
      Xcurs = (xPos - x[0]) / deltaX;  Ycurs = (yPos - y[0]) / deltaY;
      intXcurs = ceil(Xcurs); intYcurs = ceil(Ycurs);
      printf("\n x = %d",lnx);
      // se cade in un range ragionevole
      if( (intXcurs > 0) & ((intXcurs+1) < *Nx) & (intYcurs > 0) & ((intYcurs+1) < *Ny) ) {
        f00 = f[  intYcurs * (*Nx) + intXcurs ];
        f01 = f[  intYcurs * (*Nx) + (intXcurs+1) ];
        f10 = f[  (intYcurs+1) * (*Nx) + intXcurs ];
        f11 = f[  (intYcurs+1) * (*Nx) + (intXcurs+1) ];
        festx_0 = ((f01-f00)/deltaX) * (xPos-x[intXcurs]) + f00;
        festx_1 = ((f11-f10)/deltaX) * (xPos-x[intXcurs]) + f10;
        res[ lny * (*NxOut) + lnx ] = ((festx_0-festx_1)/deltaY) * (yPos-y[intYcurs]) + festx_1;
      }
    }
  }
}
*/

/*
 *  Algorithm for eroding margins in voxelcube structures
 *  cube : a pointer to the array of the cube structure
 *  nX,nY,nZ : dimensions of the cube
 *  mx,my,mz : erosion along the three axis
 *  iterator : must be set to '0', it is used to check the deep of the iteration tree.
 *  minValue : which is the minimum value that should be considered 'water' *  
 */
void erosion( double *cube, 
              int *nX, int *nY, int *nZ, 
              int *mx, int *my, int *mz, 
              int *iterator, int *minValue) {
  int x,y,z,center,ct;
  if( *iterator >= 10) return;  // just to avoid infinite loops
  if(*mx == 0 && *my ==0 && *mz ==0 ) return;
  
  
  // loop per ogni elemento del cubo
  for( z=0; z<*nZ; z++ ) {
    for( y=0; y<*nY; y++ ) {
      for( x=0; x<*nX; x++) {
        // prendi l'offset relativo al punto in esame
        center = posDecod(x,y,z,*nX,*nY,*nZ);
        // se != 0 vediamo l'intorno
        if(cube[center]>(*minValue+1)) {
          if( *mx>0 ){
            //if(x==0 || x==*nX) cube[center]=-1;
            if(x==0 || x>=(*nX-1)) cube[center]=-1;
            else if( cube[posDecod(x-1,y,z,*nX,*nY,*nZ)]<(*minValue+1) || cube[posDecod(x+1,y,z,*nX,*nY,*nZ)]<(*minValue+1)  ) cube[center]=-1;
          }
          if( *my>0 ){
            // if(y==0 || y==*nY) cube[center]=-1;
            if(y==0 || y==(*nY-1)) cube[center]=-1;
            else if( cube[posDecod(x,y-1,z,*nX,*nY,*nZ)]<(*minValue+1) || cube[posDecod(x,y+1,z,*nX,*nY,*nZ)]<(*minValue+1)  ) cube[center]=-1;
          }        
          if( *mz>0 ){
            // if(z==0 || z==*nZ) cube[center]=-1;
            if(z==0 || z==(*nZ-1)) cube[center]=-1;
            else if( cube[posDecod(x,y,z-1,*nX,*nY,*nZ)]<(*minValue+1) || cube[posDecod(x,y,z+1,*nX,*nY,*nZ)]<(*minValue+1)  )  cube[center]=-1;
          }
        }
      }
    }
  }
  // decrementa i vincoli sui margini
  if( *mx>0 ) *mx = *mx - 1;
  if( *my>0 ) *my = *my - 1;  
  if( *mz>0 ) *mz = *mz - 1;
  // trasforma tutti i '-1' in '0' per l'iterazione successiva
  //for(ct=0; ct<= posDecod((*nX-1),(*nY-1),(*nZ-1),*nX,*nY,*nZ); ct++) {
  for(ct=0; ct<= ((*nX)*(*nY)*(*nZ)-1); ct++) {
    if(cube[ct]==-1) cube[ct]=(*minValue);
  }
  // rilancia ricorsivamente
  *iterator = *iterator + 1;
  erosion( cube, nX, nY, nZ, mx, my, mz, iterator, minValue );
}


/*
void executeCMDLine( char **stringa, int *strLength ) {
  int i;
  for( i = 0; i< (*strLength) ; i++ ) {
    if( stringa[0][i] == '/' ) stringa[0][i] = '\\';
  }
  int res = system( stringa[0] );
}*/

/*
void c_getInterpolatedSlice(
    int *nx_PT, int *ny_PT,
    int *slice_b_PT, int *slice_t_PT,
    double *IOM_b_PT, double *IOM_t_PT,
    int *nx_CT, int *ny_CT,
    int *slice_CT,
    double *IOM_CT,
    double *image_b_PT, double *image_t_PT,
    double *returnMatrix ) {
  
  struct _c_data *punt;
  struct pointInSpace pis;
  struct pointInSpace *p_pis_b_PT;
  struct pointInSpace *p_pis_t_PT;
  int ct, nx, ny;
  struct DICOM_OrientationMatrix IOM_b_PT_STR, IOM_t_PT_STR, IOM_CT_STR;
  int dimensioniPET,dimensioniCT;
  double valoreCalcolato;
  
  dimensioniPET = (*ny_PT)*(*nx_PT);
  dimensioniCT = (*ny_CT)*(*nx_CT);
  
  return;
  
  //Copia le IOM nelle apposite strutture dati
  IOM_b_PT_STR.a11 = IOM_b_PT[0]; IOM_b_PT_STR.a21 = IOM_b_PT[1]; IOM_b_PT_STR.a31 = IOM_b_PT[2];
  IOM_b_PT_STR.a12 = IOM_b_PT[4]; IOM_b_PT_STR.a22 = IOM_b_PT[5]; IOM_b_PT_STR.a32 = IOM_b_PT[6];
  IOM_b_PT_STR.Sx = IOM_b_PT[12]; IOM_b_PT_STR.Sy = IOM_b_PT[13]; IOM_b_PT_STR.Sz = IOM_b_PT[14];    
  IOM_t_PT_STR.a11 = IOM_t_PT[0]; IOM_t_PT_STR.a21 = IOM_t_PT[1]; IOM_t_PT_STR.a31 = IOM_t_PT[2];
  IOM_t_PT_STR.a12 = IOM_t_PT[4]; IOM_t_PT_STR.a22 = IOM_t_PT[5]; IOM_t_PT_STR.a32 = IOM_t_PT[6];
  IOM_t_PT_STR.Sx = IOM_t_PT[12]; IOM_t_PT_STR.Sy = IOM_t_PT[13]; IOM_t_PT_STR.Sz = IOM_t_PT[14];    
  IOM_CT_STR.a11 = IOM_CT[0]; IOM_CT_STR.a21 = IOM_CT[1]; IOM_CT_STR.a31 = IOM_CT[2];
  IOM_CT_STR.a12 = IOM_CT[4]; IOM_CT_STR.a22 = IOM_CT[5]; IOM_CT_STR.a32 = IOM_CT[6];
  IOM_CT_STR.Sx = IOM_CT[12]; IOM_CT_STR.Sy = IOM_CT[13]; IOM_CT_STR.Sz = IOM_CT[14];    
  
  // Alloca un puntatore alla struttura _c_data
  punt = (struct _c_data *)calloc(1,sizeof(struct _c_data));  if( punt == NULL ) return;  
  // Pulisci la matrice di destinazione
  for(int ct = 0; ct < ( (*nx_CT )*(*ny_CT) - 1 ); ct++ ) returnMatrix[ ct ] = 0;  
  
  p_pis_b_PT = (struct pointInSpace *)calloc( dimensioniPET, sizeof(struct pointInSpace) );  if( p_pis_b_PT == NULL ) return;
  p_pis_t_PT = (struct pointInSpace *)calloc( dimensioniPET, sizeof(struct pointInSpace) );  if( p_pis_t_PT == NULL ) return;    

  // calcola i punti <x,y,z> per tutti i punti delle slice PET e CT
  for( ny=0; ny < *ny_PT; ny++ ) {    
    for( nx=0, ct=0; nx < *nx_PT; nx++, ct++ ) {      
      p_pis_b_PT[ ct ] = get3DPosFromNxNy( nx,  ny, IOM_b_PT_STR );
      p_pis_t_PT[ ct ] = get3DPosFromNxNy( nx,  ny, IOM_t_PT_STR );    
    }  
  }
  
  // ora scorri sulle CT (i punti da calcolare)  
  struct pointInSpace NxNy_b, NxNy_t, TDxyz_b,TDxyz_bnext,TDxyz_t;
  for( ny=0; ny < *ny_CT; ny++ ) {    
    printf("\ny = %d",ny);
    for( nx=0; nx < *nx_CT; nx++ ) {      
      pis = get3DPosFromNxNy( nx,  ny, IOM_CT_STR );

      NxNy_b = getNxNyPos3D( pis.x, pis.y, pis.z, IOM_b_PT_STR );
      NxNy_t = getNxNyPos3D( pis.x, pis.y, pis.z, IOM_t_PT_STR );

      TDxyz_b = get3DPosFromNxNy(  NxNy_b.x, NxNy_b.y, IOM_b_PT_STR);
      TDxyz_t = get3DPosFromNxNy(  NxNy_t.x, NxNy_t.y, IOM_t_PT_STR);
      return;
      TDxyz_bnext = get3DPosFromNxNy(  NxNy_b.x+1, NxNy_b.y+1, IOM_b_PT_STR);
      return;
      printf("\n ");
      printf("\n TDxyz_b.x=%lf, TDxyz_b.y=%lf, TDxyz_b.z=%lf", TDxyz_b.x, TDxyz_b.y, TDxyz_b.z);
      printf("\n TDxyz_t.x=%lf, TDxyz_t.y=%lf, TDxyz_t.z=%lf", TDxyz_t.x, TDxyz_t.y, TDxyz_t.z);
      printf("\n dx1x0=%lf, dy1y0=%lf, dz1z0=%lf", abs(TDxyz_bnext.x - TDxyz_b.x),abs(TDxyz_bnext.y - TDxyz_b.y),abs(TDxyz_t.z - TDxyz_b.z) );
      printf("\n x=%lf, y=%lf, z=%lf", pis.x, pis.y, pis.z );
      
      valoreCalcolato = _c_TrilinearInterpolation(
        image_b_PT[ (int)floor(NxNy_b.y) * (*nx_PT) + (int)floor(NxNy_b.x) ],  //x0y0z0 (sample value)
        image_t_PT[ (int)floor(NxNy_b.y) * (*nx_PT) + (int)floor(NxNy_b.x) ],  //x0y0z1 (sample value)
        image_b_PT[ (int)floor(NxNy_t.y) * (*nx_PT) + (int)floor(NxNy_b.x) ],  //x0y1z0 (sample value)
        image_t_PT[ (int)floor(NxNy_b.y) * (*nx_PT) + (int)floor(NxNy_b.x) ],  //x0y1z1 (sample value)                  
        image_b_PT[ (int)floor(NxNy_b.y) * (*nx_PT) + (int)floor(NxNy_t.x) ],  //x1y0z0 (sample value)
        image_t_PT[ (int)floor(NxNy_b.y) * (*nx_PT) + (int)floor(NxNy_t.x) ],  //x1y0z1 (sample value)
        image_b_PT[ (int)floor(NxNy_t.y) * (*nx_PT) + (int)floor(NxNy_t.x) ],  //x1y1z0 (sample value)                  
        image_t_PT[ (int)floor(NxNy_t.y) * (*nx_PT) + (int)floor(NxNy_t.x) ],  //x1y1z1 (sample value)                  
        TDxyz_b.x, //x0,
        TDxyz_b.y, //y0,
        TDxyz_b.z, //z0,
        abs(TDxyz_bnext.x - TDxyz_b.x), //dx1x0,
        abs(TDxyz_bnext.y - TDxyz_b.y), //dy1y0,
        abs(TDxyz_t.z - TDxyz_b.z), //dz1z0,
        pis.x,pis.y,pis.z
      );
      exit(0);
    }
    exit(0);
  }
  

}
 */