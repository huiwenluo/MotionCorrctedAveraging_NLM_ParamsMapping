#include "mex.h"
#include "string.h"
#include "math.h"
#include "stdlib.h"
#include "levmar-2.5/levmar.h"
#include "levmar-2.5/lm.h"


/* Some helper functions */
double dmin(double a, double b);
int imax(int a, int b);
void  createSignal_magn_1R2star(double *shat, double *curs, int m, int nte, void *adata);


#define PI 3.14159265
#define GYRO 42.58
#define MAX_ITER 10

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
  const char **fnames1, **fnames2, **fnames3;       /* pointers to field names */
  const mwSize *dims;
  mxArray    *tmp, *fout, *species;
  char       *pdata=NULL;

  /* Fat water parameters */
  int nx,ny,nte,fnum,nf;
  const mwSize *imDims;
  double *ims, *te, fieldStrength;
  double waterAmp, *fF, *relAmps;
  double *initW, *initF, *initR2, *outW, *outF, *outR2;
  
  /* Internal parameters */
  int kx,ky,kt,kf;
  double *curs, *sfr, *sfi, *swr, *swi, initEst[3], *work, *covar, *adata;



  /* check proper input and output */
  if(nrhs!=3)
    mexErrMsgTxt("Three inputs required: imDataParams, algoParams, initParams.");
  else if(nlhs > 1)
    mexErrMsgTxt("Too many output arguments (outParams needed)");
  else if(!mxIsStruct(prhs[0]) || !mxIsStruct(prhs[1]) || !mxIsStruct(prhs[2]) )
    mexErrMsgTxt("Input must be a structure.");

  /* Get imDataParams */
  fnum = mxGetFieldNumber(prhs[0], "images");
  if( fnum>=0 ){
    ims = mxGetPr(mxGetFieldByNumber(prhs[0], 0, fnum));
    imDims = mxGetDimensions(mxGetFieldByNumber(prhs[0], 0, fnum));
    nx = (int)(imDims[0]);
    ny = (int)(imDims[1]);
    nte = (int)(imDims[2]);
  } else {
    mexErrMsgTxt("Data input error: images not correctly specified");
  }

  fnum = mxGetFieldNumber(prhs[0], "TE");
  if( fnum>=0 ){
    te = mxGetPr(mxGetFieldByNumber(prhs[0], 0, fnum));
    imDims = mxGetDimensions(mxGetFieldByNumber(prhs[0], 0, fnum));
    nte = imax((int)(imDims[0]),(int)(imDims[1]));
  } else {
    mexErrMsgTxt("Data input error: TE not correctly specified");
  }

  fnum = mxGetFieldNumber(prhs[0], "FieldStrength");
  if( fnum>=0 ){
    fieldStrength = (double)(mxGetScalar(mxGetFieldByNumber(prhs[0], 0, fnum)));
      } else {
    fieldStrength = 1.5;
  }

  mexPrintf("%s%d%s%d%s%d%s%f\n", "nx: ", nx, ", ny: ", ny, ", nte :", nte, ", field :", fieldStrength);
  mexPrintf("%s%f%s%f%s%f\n", "im1: ", ims[0], ", te1: ", te[0], ", te2 :", te[1]);
  
  /* Get algoParams */
  fnum = mxGetFieldNumber(prhs[1], "species");
  if( fnum>=0 ){
    species = mxGetFieldByNumber(prhs[1], 0, fnum);
    waterAmp = (double)(mxGetScalar(mxGetField(species, 0, "relAmps")));
    relAmps = mxGetPr(mxGetField(species, 1, "relAmps"));
    fF = mxGetPr(mxGetField(species, 1, "frequency"));
    imDims = mxGetDimensions(mxGetField(species, 1, "relAmps"));
    nf = imax((int)(imDims[0]),(int)(imDims[1]));
  } else {
    mexErrMsgTxt("Algo input error: species not correctly specified");
  }

  /* Get initParams */
  fnum = mxGetFieldNumber(prhs[2], "species");
  if( fnum>=0 ){
    species = mxGetFieldByNumber(prhs[2], 0, fnum);
    initW = mxGetPr(mxGetField(species, 0, "amps"));
    initF = mxGetPr(mxGetField(species, 1, "amps")); 
  } else {
    mexErrMsgTxt("Init input error: init amplitudes not correctly specified");
  }

  fnum = mxGetFieldNumber(prhs[2], "r2starmap");
  if( fnum>=0 ){
    initR2 = mxGetPr(mxGetFieldByNumber(prhs[2], 0, fnum));
  } else {
    mexErrMsgTxt("Init input error: init R2* not correctly specified");
  }


  /* Go process */
  /* Initialize output structure */
  plhs[0] = mxDuplicateArray(prhs[2]);
  fnum = mxGetFieldNumber(plhs[0], "species");
  if( fnum>=0 ){
    species = mxGetFieldByNumber(plhs[0], 0, fnum);
    outW = mxGetPr(mxGetField(species, 0, "amps"));
    outF = mxGetPr(mxGetField(species, 1, "amps")); 
  } else {
    mexErrMsgTxt("Init input error: init amplitudes not correctly specified");
  }

  fnum = mxGetFieldNumber(plhs[0], "r2starmap");
  if( fnum>=0 ){
    initR2 = mxGetPr(mxGetFieldByNumber(plhs[0], 0, fnum));
  } else {
    mexErrMsgTxt("Init input error: init R2* not correctly specified");
  }


  curs = (double *)malloc(nte*sizeof(double));
  sfr = (double *)malloc(nte*sizeof(double));
  sfi = (double *)malloc(nte*sizeof(double));
  swr = (double *)malloc(nte*sizeof(double));
  swi = (double *)malloc(nte*sizeof(double));

  adata = (double *)malloc(5*nte*sizeof(double));

  work = (double *)malloc(nte*3*2*sizeof(double));
  covar = (double *)malloc(3*3*2*sizeof(double));

  /* Initialize water/fat signal models */
  for(kf=0;kf<nf;kf++) {
    fF[kf] = fF[kf]*GYRO*fieldStrength;
  }
  for(kt=0;kt<nte;kt++) {
    swr[kt] = waterAmp;
    swi[kt] = 0.0;
    sfr[kt] = 0.0;
    sfi[kt] = 0.0;
    for(kf=0;kf<nf;kf++) {
      
      sfr[kt] = sfr[kt] + relAmps[kf]*cos(2*PI*te[kt]*fF[kf]);
      sfi[kt] = sfi[kt] + relAmps[kf]*sin(2*PI*te[kt]*fF[kf]);
    }    

    mexPrintf("%s%f%s%f%s%f%s%f\n", "swr: ", swr[kt],", swr: ", swi[kt]," , sfr: ", sfr[kt],", sfi: ", sfi[kt]);


    adata[kt] = te[kt];
    adata[kt+nte] = swr[kt];
    adata[kt+2*nte] = swi[kt];
    adata[kt+3*nte] = sfr[kt];
    adata[kt+4*nte] = sfi[kt];

  }
  
  /* Loop over all pixels */
  for(kx=0;kx<nx;kx++) {
    for(ky=0;ky<ny;ky++) {
      
      /* Get signal at current voxel */
      for(kt=0;kt<nte;kt++) {
	curs[kt] = ims[kx + ky*nx + kt*nx*ny];
      }
      
      initEst[0] = initW[kx + ky*nx];
      initEst[1] = initF[kx + ky*nx];
      initEst[2] = initR2[kx + ky*nx];
      

      /* Do fitting */
      dlevmar_dif(createSignal_magn_1R2star,initEst,curs,3,nte,MAX_ITER,NULL,NULL,work,covar,adata);
	
	}
  }
  
  free(curs);
  free(sfr);
  free(sfi);
  free(swr);
  free(swi);
  free(work);
  free(covar);
  return;
}



void  createSignal_magn_1R2star(double *shat, double *curs, int m, int nte, void *adata) {
  int kt;

  for(kt=0;kt<nte;kt++) {
    shat[kt] = 0.0;
  }
  
  
  
}


double dmin(double a, double b) {
  if( a < b ) {
    return a;
  } else {
    return b;
  }
}

int imax(int a, int b) {
  if( a < b ) {
    return b;
  } else {
    return a;
  }
}
