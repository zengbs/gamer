#include<stdio.h>
#include<stdlib.h>
#include<limits.h>
//#include"Prototype.h"

template <typename T>
int Mis_BinarySearch_Real( const T Array[], int Min, int Max, const T Key )
{

// initial check
#  ifdef GAMER_DEBUG
   if ( Min < 0 )    Aux_Error( ERROR_INFO, "incorrect input parameter \"Min (%d) < 0\" !!\n", Min );
   if ( Max <= Min ) Aux_Error( ERROR_INFO, "incorrect input parameters \"Max (%d) <= Min (%d)\" !!\n", Max, Min );
#  endif


// check whether the input key lies outside the target range
   if ( Key <  Array[Min] )   return Min-1;
   if ( Key >= Array[Max] )   return Max;


// binary search
   int Idx = -2;

   while (  ( Idx=(Min+Max)/2 ) != Min  )
   {
      if   ( Array[Idx] > Key )  Max = Idx;
      else                       Min = Idx;
   }


// check whether the found array index is correct
#  ifdef GAMER_DEBUG
   if ( Idx < Min  ||  Idx >= Max )
      Aux_Error( ERROR_INFO, "incorrect output index (Idx %d, Min %d, Max%d) !!\n", Idx, Min, Max );

   if (  Array[Idx] > Key  ||  Array[Idx+1] <= Key )
      Aux_Error( ERROR_INFO, "incorrect output index (Idx %d, ValueL %14.7e, ValueR %14.7e, Key %14.7e) !!\n",
                 Idx, Array[Idx], Array[Idx+1], Key );
#  endif


   return Idx;

} // FUNCTION : Mis_BinarySearch_Real



// explicit template instantiation
template int Mis_BinarySearch_Real <float>  ( const float  Array[], int Min, int Max, const float  Key );
template int Mis_BinarySearch_Real <double> ( const double Array[], int Min, int Max, const double Key );

/*----------------------------------------------------------------------------*/
/*! \fn void*** calloc_3d_array(size_t nt, size_t nr, size_t nc, size_t size)
 *  *  *  \brief Construct 3D array = array[nt][nr][nc]  */
void ***calloc_3d_array (size_t nt, size_t nr, size_t nc, size_t size)
{
  void ***array;
  size_t i, j;

  if ((array = (void ***) calloc (nt, sizeof (void **))) == NULL)
    {
      printf ("[calloc_3d] failed to allocate memory for %d 1st-pointers\n", (int) nt);
      return NULL;
    }

  if ((array[0] = (void **) calloc (nt * nr, sizeof (void *))) == NULL)
    {
      printf ("[calloc_3d] failed to allocate memory for %d 2nd-pointers\n", (int) (nt * nr));
      free ((void *) array);
      return NULL;
    }

  for (i = 1; i < nt; i++)
    {
      array[i] = (void **) ((unsigned char *) array[0] + i * nr * sizeof (void *));
    }

  if ((array[0][0] = (void *) calloc (nt * nr * nc, size)) == NULL)
    {
      printf ("[calloc_3d] failed to alloc. memory (%d X %d X %d of size %d)\n",
      (int) nt, (int) nr, (int) nc, (int) size);
      free ((void *) array[0]);
      free ((void *) array);
      return NULL;
    }
  for (j = 1; j < nr; j++)
    {
      array[0][j] = (void **) ((unsigned char *) array[0][j - 1] + nc * size);
    }
for (i = 1; i < nt; i++)
    {
      array[i][0] = (void **) ((unsigned char *) array[i - 1][0] + nr * nc * size);
      for (j = 1; j < nr; j++)
        {
          array[i][j] = (void **) ((unsigned char *) array[i][j - 1] + nc * size);
        }
    }

  return array;
}

/*! \fn void free_3d_array(void *array)
 *  *  \brief Free memory used by 3D array  */
void free_3d_array(void *array)
{
  void ***ta = (void ***)array;

  free(ta[0][0]);
  free(ta[0]);
  free(array);
}

/* boxsize: the number of total cells excluding ghost cells */

//void GetPeriodicIdx(int i, int *ii, int boxsize, int ghostsize)
//{
//  i>=ghostsize ? *ii = (i-ghostsize)%boxsize: *ii = i-ghostsize+boxsize;
//}


double TrilinearInterpolation(double *FieldAtVertices, double *xyz000, double *dxyz, double *xyz)
{
  double x1, y1, z1, x0, y0, z0, xd, yd, zd, x, y, z;
  double c000, c001, c010, c100, c011, c101, c110, c111, c00, c01, c10, c11, c0, c1, c;

  x0 = xyz000[0];
  y0 = xyz000[1];
  z0 = xyz000[2];

  x1 = xyz000[0] + dxyz[0];
  y1 = xyz000[1] + dxyz[1];
  z1 = xyz000[2] + dxyz[2];
  
  x = xyz[0];
  y = xyz[1];
  z = xyz[2];

  c000 = FieldAtVertices[0];
  c001 = FieldAtVertices[1];
  c010 = FieldAtVertices[2];
  c100 = FieldAtVertices[3];
  c011 = FieldAtVertices[4];
  c101 = FieldAtVertices[5];
  c110 = FieldAtVertices[6];
  c111 = FieldAtVertices[7];

  xd = (x-x0)/(x1-x0);
  yd = (y-y0)/(y1-y0);
  zd = (z-z0)/(z1-z0);

  c00 = c000*(1.0-xd) + c100*xd;
  c01 = c001*(1.0-xd) + c101*xd;
  c10 = c010*(1.0-xd) + c110*xd;
  c11 = c011*(1.0-xd) + c111*xd;

  c0  = c00*(1.0-yd) + c10*yd;
  c1  = c01*(1.0-yd) + c11*yd;

  c = c0*(1.0-zd) + c1*zd;

  return c;
}


double Interpolation_UM_IC( char *FileName, double x, double y, double z )
{


  FILE *pFile;
  long lSize;
  double *buffer;
  size_t result;

  pFile = fopen ( FileName , "rb" );
  if (pFile==NULL) {fputs ("File error",stderr); exit (1);}

  // obtain file size:
  fseek (pFile , 0 , SEEK_END);
  lSize = ftell (pFile);
  rewind (pFile);

  // allocate memory to contain the whole file:
  buffer = (double*) calloc (lSize,sizeof(double));
  if (buffer == NULL) {fputs ("Memory error",stderr); exit (2);}

  // copy the file into the buffer:
  result = fread (buffer,1,lSize,pFile);
  if (result != lSize) {fputs ("Reading error",stderr); exit (3);}

  int HeaderSize = (int)buffer[0];
  printf("HeaderSize=%d\n", HeaderSize);
  int Nx = (int)buffer[1];
  int Ny = (int)buffer[2];
  int Nz = (int)buffer[3];

  if (5*Nx*Ny*Nz > INT_MAX) {printf("integer overflow!!\n"); exit(0);}

  double dx = buffer[4];
  double dy = buffer[5];
  double dz = buffer[6];
  if ( x<0 || x>Nx*dx ) {printf("x is soutside box!\n"); exit(0);}
  if ( y<0 || y>Ny*dy ) {printf("y is soutside box!\n"); exit(0);}
  if ( z<0 || z>Nz*dz ) {printf("z is soutside box!\n"); exit(0);}

  double dxyz[3] = {dx, dy, dz};
 
  int ghostsize = 1;
 
  double ***Rhoo = (double***)calloc_3d_array((size_t)(Nx+2*ghostsize), (size_t)(Ny+2*ghostsize), (size_t)(Nz+2*ghostsize), sizeof(double));
  double ***VelX = (double***)calloc_3d_array((size_t)(Nx+2*ghostsize), (size_t)(Ny+2*ghostsize), (size_t)(Nz+2*ghostsize), sizeof(double));
  double ***VelY = (double***)calloc_3d_array((size_t)(Nx+2*ghostsize), (size_t)(Ny+2*ghostsize), (size_t)(Nz+2*ghostsize), sizeof(double));
  double ***VelZ = (double***)calloc_3d_array((size_t)(Nx+2*ghostsize), (size_t)(Ny+2*ghostsize), (size_t)(Nz+2*ghostsize), sizeof(double));
  double ***Pres = (double***)calloc_3d_array((size_t)(Nx+2*ghostsize), (size_t)(Ny+2*ghostsize), (size_t)(Nz+2*ghostsize), sizeof(double));

  double *X = (double*)calloc((size_t)(Nx+2*ghostsize),sizeof(double));
  double *Y = (double*)calloc((size_t)(Ny+2*ghostsize),sizeof(double));
  double *Z = (double*)calloc((size_t)(Nz+2*ghostsize),sizeof(double));


  for (int i=-ghostsize;i<Nx+ghostsize;i++) X[i+ghostsize] = (0.5+(double)i)*dx;
  for (int i=-ghostsize;i<Ny+ghostsize;i++) Y[i+ghostsize] = (0.5+(double)i)*dy;
  for (int i=-ghostsize;i<Nz+ghostsize;i++) Z[i+ghostsize] = (0.5+(double)i)*dz;
  

  int NX = Nx+2*ghostsize;
  int NY = Ny+2*ghostsize;
  int NZ = Nz+2*ghostsize;


  for (int c=0;c<5*NX*NY*NZ;c++){
    int i, j, k, cc;
    int ii, jj, kk;

    cc = c%(NX*NY*NZ);
    i = (cc - cc%(NY*NZ)) / (NY*NZ);
    j = ((cc - cc%NZ) / NZ) % NY;
    k = cc%NZ;

    if (          0 <= c && c <   NX*NY*NZ ) Rhoo[i][j][k] = buffer[c+HeaderSize];
    if (   NX*NY*NZ <= c && c < 2*NX*NY*NZ ) VelX[i][j][k] = buffer[c+HeaderSize];
    if ( 2*NX*NY*NZ <= c && c < 3*NX*NY*NZ ) VelY[i][j][k] = buffer[c+HeaderSize];
    if ( 3*NX*NY*NZ <= c && c < 4*NX*NY*NZ ) VelZ[i][j][k] = buffer[c+HeaderSize];
    if ( 4*NX*NY*NZ <= c && c < 5*NX*NY*NZ ) Pres[i][j][k] = buffer[c+HeaderSize];
  }


  int Idx = Mis_BinarySearch_Real(X, 0, NX-1, x);
  int Jdx = Mis_BinarySearch_Real(Y, 0, NY-1, y);
  int Kdx = Mis_BinarySearch_Real(Z, 0, NZ-1, z);

  if (Idx<0 || Idx > NX-1){ printf("Idx=%d is out of range!\n", Idx); exit(0); }
  if (Jdx<0 || Jdx > NY-1){ printf("Jdx=%d is out of range!\n", Jdx); exit(0); }
  if (Kdx<0 || Kdx > NZ-1){ printf("Kdx=%d is out of range!\n", Kdx); exit(0); }


  double Vertex000 = Rhoo[Idx  ][Jdx  ][Kdx  ];
  double Vertex001 = Rhoo[Idx  ][Jdx  ][Kdx+1];
  double Vertex010 = Rhoo[Idx  ][Jdx+1][Kdx  ];
  double Vertex100 = Rhoo[Idx+1][Jdx  ][Kdx  ];
  double Vertex011 = Rhoo[Idx  ][Jdx+1][Kdx+1];
  double Vertex101 = Rhoo[Idx+1][Jdx  ][Kdx+1];
  double Vertex110 = Rhoo[Idx+1][Jdx+1][Kdx  ];
  double Vertex111 = Rhoo[Idx+1][Jdx+1][Kdx+1];

  double FieldAtVertices[8]
         = {Vertex000, Vertex001, Vertex010, Vertex100, Vertex011, Vertex101, Vertex110, Vertex111};

  double xyz000[3] = {X[Idx], Y[Jdx], Z[Kdx]};
  double xyz[3] = {x, y, z};

  double InterpolatedField = TrilinearInterpolation(FieldAtVertices, xyz000, dxyz, xyz);

  // free memory
  fclose (pFile);
  free (buffer);
  free_3d_array(Rhoo);
  free_3d_array(VelX);
  free_3d_array(VelY);
  free_3d_array(VelZ);
  free_3d_array(Pres);
  free(X);
  free(Y);
  free(Z);

  return InterpolatedField;
}

int main()
{ 
  double x = 0.1;
  double y = 0.2; 
  double z = 0.3; 
  char FileName[] = "Fluid3D.bin";
  double  FineField;
  FineField = Interpolation_UM_IC( FileName, x, y, z );

  printf("FineField=%e\n", FineField);

  return 0;
}
