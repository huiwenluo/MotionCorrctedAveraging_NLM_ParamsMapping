#include "mex.h"
#include "string.h"
#include "math.h"
#include "stdlib.h"
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>

/* Some helper functions */
double dmin(double a, double b);
double complexPhase(double x, double y);
int imax(int a, int b);
int  createSignal_mixed_1R2star_f(const gsl_vector * x, void *data, gsl_vector * f);
int  createSignal_mixed_1R2star_df(const gsl_vector * x, void *data, gsl_matrix * J);
int  createSignal_mixed_1R2star_fdf(const gsl_vector * x, void *data,gsl_vector * f, gsl_matrix * J);
void print_state (size_t iter, gsl_multifit_fdfsolver * s);

#define PI 3.14159265
#define GYRO 42.58
#define MAX_ITER 50


struct data {
  int nte;
  double *te;
  double *cursr;
  double *cursi;
  double *swr;
  double *swi;
  double *sfr;
  double *sfi;
};



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
  double *imsr, *imsi, *te, fieldStrength, clockwise;
  double waterAmp, *fPPM,*fF, *relAmps;
  double *initWr, *initFr, *initWi, *initFi, *initR2, *initFieldmap, *outWr, *outFr, *outWi, *outFi, *outR2, *outFieldmap;
  
  /* Internal parameters */
  int kx,ky,kt,kf;
  double *cursr, *cursi, *sfr, *sfi, *swr, *swi, initEst[6];


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
    imsr = mxGetPr(mxGetFieldByNumber(prhs[0], 0, fnum));
    imsi = mxGetPi(mxGetFieldByNumber(prhs[0], 0, fnum));
    imDims = mxGetDimensions(mxGetFieldByNumber(prhs[0], 0, fnum));
    nx = (int)(imDims[0]);
    ny = (int)(imDims[1]);
    nte = (int)(imDims[4]);
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

  fnum = mxGetFieldNumber(prhs[0], "PrecessionIsClockwise");
  if( fnum>=0 ){
    clockwise = (double)(mxGetScalar(mxGetFieldByNumber(prhs[0], 0, fnum)));
      } else {
    clockwise = 1;
  }

/*     mexPrintf("%s%d%s%d%s%d%s%f\n", "nx: ", nx, ", ny: ", ny, ", nte :", nte, ", field :", fieldStrength); */
/*     mexPrintf("%s%f%s%f%s%f\n", "im1: ", imsr[0], ", te1: ", te[0], ", te2 :", te[1]); */
  
  /* Get algoParams */
  fnum = mxGetFieldNumber(prhs[1], "species");
  if( fnum>=0 ){
    species = mxGetFieldByNumber(prhs[1], 0, fnum);
    waterAmp = (double)(mxGetScalar(mxGetField(species, 0, "relAmps")));
    relAmps = mxGetPr(mxGetField(species, 1, "relAmps"));
    fPPM = mxGetPr(mxGetField(species, 1, "frequency"));
    imDims = mxGetDimensions(mxGetField(species, 1, "relAmps"));
    nf = imax((int)(imDims[0]),(int)(imDims[1]));
  } else {
    mexErrMsgTxt("Algo input error: species not correctly specified");
  }

  /* Get initParams */
  fnum = mxGetFieldNumber(prhs[2], "species");
  if( fnum>=0 ){
    species = mxGetFieldByNumber(prhs[2], 0, fnum);
    initWr = mxGetPr(mxGetField(species, 0, "amps"));
    initFr = mxGetPr(mxGetField(species, 1, "amps")); 
    initWi = mxGetPi(mxGetField(species, 0, "amps"));
    initFi = mxGetPi(mxGetField(species, 1, "amps")); 
  } else {
    mexErrMsgTxt("Init input error: init amplitudes not correctly specified");
  }


/*   mexPrintf("%s%f%s%f%s%f\n", "initwi1: ", initWi[0], ", te1: ", te[0], ", te2 :", te[1]);  */


  fnum = mxGetFieldNumber(prhs[2], "r2starmap");
  if( fnum>=0 ){
    initR2 = mxGetPr(mxGetFieldByNumber(prhs[2], 0, fnum));
  } else {
    mexErrMsgTxt("Init input error: init R2* not correctly specified");
  }

  fnum = mxGetFieldNumber(prhs[2], "fieldmap");
  if( fnum>=0 ){
    initFieldmap = mxGetPr(mxGetFieldByNumber(prhs[2], 0, fnum));
  } else {
    mexErrMsgTxt("Init input error: init R2* not correctly specified");
  }


  /* Go process */
  /* Initialize output structure */
  plhs[0] = mxDuplicateArray(prhs[2]);
  fnum = mxGetFieldNumber(plhs[0], "species");
  if( fnum>=0 ){
    species = mxGetFieldByNumber(plhs[0], 0, fnum);
    outWr = mxGetPr(mxGetField(species, 0, "amps"));
    outFr = mxGetPr(mxGetField(species, 1, "amps")); 
    outWi = mxGetPi(mxGetField(species, 0, "amps"));
    outFi = mxGetPi(mxGetField(species, 1, "amps")); 
  } else {
    mexErrMsgTxt("Init input error: init amplitudes not correctly specified");
  }

  fnum = mxGetFieldNumber(plhs[0], "r2starmap");
  if( fnum>=0 ){
    outR2 = mxGetPr(mxGetFieldByNumber(plhs[0], 0, fnum));
  } else {
    mexErrMsgTxt("Init input error: init R2* not correctly specified");
  }

  fnum = mxGetFieldNumber(plhs[0], "fieldmap");
  if( fnum>=0 ){
    outFieldmap = mxGetPr(mxGetFieldByNumber(plhs[0], 0, fnum));
  } else {
    mexErrMsgTxt("Init input error: init fieldmap not correctly specified");
  }

  cursr = (double *)malloc(nte*sizeof(double));
  cursi = (double *)malloc(nte*sizeof(double));
  sfr = (double *)malloc(nte*sizeof(double));
  sfi = (double *)malloc(nte*sizeof(double));
  swr = (double *)malloc(nte*sizeof(double));
  swi = (double *)malloc(nte*sizeof(double));

  /* Initialize gsl stuff   */
  const gsl_multifit_fdfsolver_type *T;
  gsl_multifit_fdfsolver *s;

  int status;
  size_t i, iter = 0;

  const size_t n = nte*2 - 1;
  const size_t p = 5;

  gsl_matrix *covar = gsl_matrix_alloc (p, p);

  struct data d;

  d.nte = nte;
  d.cursr = cursr;
  d.cursi = cursi;
  d.te = te;
  d.swr = swr;
  d.swi = swi;
  d.sfr = sfr;
  d.sfi = sfi;
  
  
  gsl_multifit_function_fdf f; 
  gsl_vector_view x = gsl_vector_view_array (initEst, p);
  const gsl_rng_type * type;
  gsl_rng * r;

  gsl_rng_env_setup();

  type = gsl_rng_default;
  r = gsl_rng_alloc (type);

  f.f = &createSignal_mixed_1R2star_f;
  f.df = &createSignal_mixed_1R2star_df;
  f.fdf = &createSignal_mixed_1R2star_fdf;
  f.n = n;
  f.p = p;
  f.params = &d;

  T = gsl_multifit_fdfsolver_lmsder;
  s = gsl_multifit_fdfsolver_alloc (T, n, p);


/*   mexPrintf("%s%f%s%d%s%d\n", "initwi1: ", initWi[0], ", n : ", n, ", p :", p);  */

  /* Initialize water/fat signal models */
  fF = (double *)malloc(nf*sizeof(double));
  for(kf=0;kf<nf;kf++) {
    fF[kf] = fPPM[kf]*GYRO*fieldStrength;
    /*    mexPrintf("%s%f%s%f%s%f\n", "fieldSgrength: " , fieldStrength, "gyro: " , GYRO, "fF: " , fF[kf]);*/
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

/*         mexPrintf("%s%d%s%f%s%f%s%f%s%f\n", "nf: " , nf, "swr: ", swr[kt],", swi: ", swi[kt]," , sfr: ", sfr[kt],", sfi: ", sfi[kt]); */


  }
  

  double curPhi, curAmpW, curAmpF;

  /* Loop over all pixels */
  for(kx=0;kx<nx;kx++) {
    for(ky=0;ky<ny;ky++) {
      
      /* Get signal at current voxel */
      if(clockwise>0) {
	for(kt=0;kt<nte;kt++) {
	  cursr[kt] = imsr[kx + ky*nx + kt*nx*ny];
	  cursi[kt] = imsi[kx + ky*nx + kt*nx*ny];
	}
      } else {
	for(kt=0;kt<nte;kt++) {
	  cursr[kt] = imsr[kx + ky*nx + kt*nx*ny];
	  cursi[kt] = -imsi[kx + ky*nx + kt*nx*ny];
	}
      }
      
      curAmpW = sqrt(initWr[kx + ky*nx]*initWr[kx + ky*nx] + initWi[kx + ky*nx]*initWi[kx + ky*nx]);
      curAmpF = sqrt(initFr[kx + ky*nx]*initFr[kx + ky*nx] + initFi[kx + ky*nx]*initFi[kx + ky*nx]);
      
      if(curAmpW>curAmpF) {
	curPhi = complexPhase(initWr[kx + ky*nx],initWi[kx + ky*nx]);
      } else {
	curPhi = complexPhase(initFr[kx + ky*nx],initFi[kx + ky*nx]);
      }

      curPhi = complexPhase(initWr[kx + ky*nx]+initFr[kx + ky*nx],initWi[kx + ky*nx]+initFi[kx + ky*nx]);



      initEst[0] = curAmpW; 
/*       initEst[1] = initWi[kx + ky*nx]; */
      initEst[1] = curAmpF;
      initEst[2] = curPhi;
      initEst[3] = initR2[kx + ky*nx];
      initEst[4] = initFieldmap[kx + ky*nx];
      
/*       mexPrintf("%s%f%s%d%s%d\n", "initwi1ph: ", curPhi, ", kx : ", kx, ", ky :", ky);   */

      /* Do fitting */
      x = gsl_vector_view_array (initEst, p);
      gsl_multifit_fdfsolver_set (s, &f, &x.vector);
      iter = 0;
/*       print_state (iter, s);  */
      do
	{
	  iter++;
	  status = gsl_multifit_fdfsolver_iterate (s);
	  
/*  	  printf ("status = %s\n", gsl_strerror (status)); */
/*  	  print_state (iter, s);  */
	  
	  if (status)
	    break;
	  
	  status = gsl_multifit_test_delta (s->dx, s->x,1e-4, 1e-4);
	}
      while (status == GSL_CONTINUE && iter < MAX_ITER);


      

      curAmpW = gsl_vector_get (s->x, 0);
      curAmpF = gsl_vector_get (s->x, 1);
      curPhi = gsl_vector_get (s->x, 2);


      outWr[kx + ky*nx] = curAmpW*cos(curPhi);
      outWi[kx + ky*nx] = curAmpW*sin(curPhi);
      outFr[kx + ky*nx] = curAmpF*cos(curPhi);
      outFi[kx + ky*nx] = curAmpF*sin(curPhi);
      outR2[kx + ky*nx] = gsl_vector_get (s->x, 3); 
      outFieldmap[kx + ky*nx] = gsl_vector_get (s->x, 4); 
      
    }
  }
  
  free(cursr);
  free(cursi);
  free(sfr);
  free(sfi);
  free(swr);
  free(swi);
  free(fF);
  gsl_multifit_fdfsolver_free (s);
  return;
}

void
print_state (size_t iter, gsl_multifit_fdfsolver * s)
{
  printf ("iter: %3u x = % 15.8f % 15.8f % 15.8f  % 15.8f % 15.8f "
          "|f(x)| = %g\n",
          iter,
          gsl_vector_get (s->x, 0), 
          gsl_vector_get (s->x, 1),
          gsl_vector_get (s->x, 2), 
          gsl_vector_get (s->x, 3), 
          gsl_vector_get (s->x, 4),
          gsl_blas_dnrm2 (s->f));
}


int  createSignal_mixed_1R2star_f(const gsl_vector * x, void *data, gsl_vector * f) {
  size_t kt;
  double shatr, shati, shat;

  int nte = ((struct data *)data)->nte;
  double *cursr = ((struct data *)data)->cursr;
  double *cursi = ((struct data *)data)->cursi;
  double *te = ((struct data *)data)->te;
  double *swr = ((struct data *)data)->swr;
  double *swi = ((struct data *)data)->swi;
  double *sfr = ((struct data *)data)->sfr;
  double *sfi = ((struct data *)data)->sfi;

  double W = gsl_vector_get (x, 0);
  double F = gsl_vector_get (x, 1);
  double phi = gsl_vector_get (x, 2);
  double r2 = gsl_vector_get (x, 3);
  double fieldmap = gsl_vector_get (x, 4);
  
/*   mexPrintf("%s%f%s%f\n", "Wr: ", Wr, ", Wi: ", Wi ); */
/*   mexPrintf("%s%f%s%f\n", "Fr: ", Fr, ", Fi: ", Fi ); */


/*   for(kt=0;kt<nte;kt++) { */
/*     mexPrintf("%s%f%s%f%s%f%s%f\n", "swr: ", swr[kt],", swi: ", swi[kt]," , sfr: ", sfr[kt],", sfi: ", sfi[kt]); */

/*   } */

  int NUM_MAGN = 1;

  for(kt=0;kt<NUM_MAGN;kt++) {
    shat = exp(-te[kt]*r2)*sqrt((W*swr[kt] + F*sfr[kt])*(W*swr[kt] + F*sfr[kt]) + (W*swi[kt] + F*sfi[kt])*(W*swi[kt] + F*sfi[kt]));
    gsl_vector_set (f, kt,  shat - sqrt(cursr[kt]*cursr[kt] + cursi[kt]*cursi[kt]));
  }  

  for(kt=NUM_MAGN;kt<nte;kt++) {

    shatr = cos(phi + 2*PI*fieldmap*te[kt])*exp(-te[kt]*r2)*(W*swr[kt] + F*sfr[kt]) - sin(phi + 2*PI*fieldmap*te[kt])*exp(-te[kt]*r2)*(W*swi[kt] + F*sfi[kt]);
    shati = sin(phi + 2*PI*fieldmap*te[kt])*exp(-te[kt]*r2)*(W*swr[kt] + F*sfr[kt]) + cos(phi + 2*PI*fieldmap*te[kt])*exp(-te[kt]*r2)*(W*swi[kt] + F*sfi[kt]);
    
    gsl_vector_set (f, kt,  shatr - cursr[kt]);
    gsl_vector_set (f, kt+nte-NUM_MAGN,  shati - cursi[kt]);

/*     mexPrintf("%s%f%s%f\n", "sr: ", cursr[kt], ", shatr: ", shatr); */
/*     mexPrintf("%s%f%s%f\n", "si: ", cursi[kt], ", shati: ", shati); */
  }

  return GSL_SUCCESS;
}


/*
  k1 = real(sum(RA.*exp(j*2*pi*FS.*T),2));
  k2 = imag(sum(RA.*exp(j*2*pi*FS.*T),2));

  J = zeros(3,N);

  
  J = [ds1(:) , ds2(:), ds3(:)+ds4(:)];
*/

int  createSignal_mixed_1R2star_df(const gsl_vector * x, void *data, gsl_matrix * J) {
  size_t kt,kp;
  double curJ1,curJ2,curJ3,shatr,shat;
  double curJ4,curJ5,shati;
  
  int nte = ((struct data *)data)->nte;
  double *te = ((struct data *)data)->te;
  double *swr = ((struct data *)data)->swr;
  double *swi = ((struct data *)data)->swi;
  double *sfr = ((struct data *)data)->sfr;
  double *sfi = ((struct data *)data)->sfi;

  double W = gsl_vector_get (x, 0);
  double F = gsl_vector_get (x, 1);
  double phi = gsl_vector_get (x, 2);
  double r2 = gsl_vector_get (x, 3);
  double fieldmap = gsl_vector_get (x, 4);

  double expr2, sinfm, cosfm;

  int NUM_MAGN = 1;


  for(kt=0;kt<NUM_MAGN;kt++) {

    shat = sqrt((W*swr[kt] + F*sfr[kt])*(W*swr[kt] + F*sfr[kt]) + (W*swi[kt] + F*sfi[kt])*(W*swi[kt] + F*sfi[kt]));

    curJ1 = exp(-te[kt]*r2)*(W*swr[kt] + F*sfr[kt])/shat;
    gsl_matrix_set (J, kt, 0, curJ1); 

    curJ2 = exp(-te[kt]*r2)*(F*(sfr[kt]*sfr[kt] + sfi[kt]*sfi[kt]) + W*sfr[kt])/shat;
    gsl_matrix_set (J, kt, 1, curJ2); 
      
    gsl_matrix_set (J, kt, 2, 0.0); 

    curJ4 = exp(-te[kt]*r2)*(-W*W*te[kt] - W*F*sfr[kt]*te[kt] -F*F*te[kt]*(sfr[kt]*sfr[kt] + sfi[kt]*sfi[kt]) - W*F*sfr[kt]*te[kt])/shat;
    gsl_matrix_set (J, kt, 3, curJ4); 

    gsl_matrix_set (J, kt, 4, 0.0); 

  }



  for(kt=NUM_MAGN;kt<nte;kt++) {
    
    expr2 = exp(-te[kt]*r2);
    sinfm = sin(phi + 2*PI*fieldmap*te[kt]);
    cosfm = cos(phi + 2*PI*fieldmap*te[kt]);


    shatr=cosfm*expr2*(W*swr[kt] + F*sfr[kt]) - sinfm*expr2*(W*swi[kt] + F*sfi[kt]);
    shati=sinfm*expr2*(W*swr[kt] + F*sfr[kt]) + cosfm*expr2*(W*swi[kt] + F*sfi[kt]);

    curJ1 = cosfm*expr2*swr[kt] - sinfm*expr2*swi[kt];
    gsl_matrix_set (J, kt, 0, curJ1); 
    curJ1 = sinfm*expr2*swr[kt] + cosfm*expr2*swi[kt];
    gsl_matrix_set (J, kt+nte-NUM_MAGN, 0, curJ1); 

    curJ2 = cosfm*expr2*sfr[kt] - sinfm*expr2*sfi[kt];
    gsl_matrix_set (J, kt, 1, curJ2); 
    curJ2 = sinfm*expr2*sfr[kt] + cosfm*expr2*sfi[kt];
    gsl_matrix_set (J, kt+nte-NUM_MAGN, 1, curJ2); 

    curJ3 = -sinfm*expr2*(W*swr[kt] + F*sfr[kt]) - cosfm*expr2*(W*swi[kt] + F*sfi[kt]);
    gsl_matrix_set (J, kt, 2, curJ3); 
    curJ3 = cosfm*expr2*(W*swr[kt] + F*sfr[kt]) - sinfm*expr2*(W*swi[kt] + F*sfi[kt]);
    gsl_matrix_set (J, kt+nte-NUM_MAGN, 2, curJ3); 

    curJ4 = -te[kt]*shatr;
    gsl_matrix_set (J, kt, 3, curJ4); 
    curJ4 = -te[kt]*shati;
    gsl_matrix_set (J, kt+nte-NUM_MAGN, 3, curJ4); 
      
    curJ5 = 2*PI*te[kt]*(-sinfm*expr2*(W*swr[kt] + F*sfr[kt]) - cosfm*expr2*(W*swi[kt] + F*sfi[kt]));
    gsl_matrix_set (J, kt, 4, curJ5); 
    curJ5 =  2*PI*te[kt]*(cosfm*expr2*(W*swr[kt] + F*sfr[kt]) - sinfm*expr2*(W*swi[kt] + F*sfi[kt]));
    gsl_matrix_set (J, kt+nte-NUM_MAGN, 4, curJ5); 
      
    /*    mexPrintf("%s%f%s%f%s%f\n", "j1: ", curJ1, ", j2: ", curJ2 , ", j3: ", curJ3);*/

  }



/*   mexPrintf("%s\n", "J matrix = "); */
/*   for (kt=0;kt<2*nte-NUM_MAGN;kt++) { */
/*     mexPrintf("%s%f%s%f%s%f%s%f%s%f\n", " ", gsl_matrix_get(J,kt,0), ", ", gsl_matrix_get(J,kt,1), ", ", gsl_matrix_get(J,kt,2), ", ", gsl_matrix_get(J,kt,3), ", ", gsl_matrix_get(J,kt,4)); */
/*   } */
/*    mexPrintf("\n\n"); */

  return GSL_SUCCESS;
}

int  createSignal_mixed_1R2star_fdf(const gsl_vector * x, void *data,gsl_vector * f, gsl_matrix * J) {
  createSignal_mixed_1R2star_f (x, data, f);
  createSignal_mixed_1R2star_df (x, data, J);
  return GSL_SUCCESS;
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


double complexPhase(double x, double y) {
  double curPhi = 0.0; 


  if(x>0) {
    curPhi = atan(y/x);
  } else if(x<0) {
    if(y>=0) {
      curPhi = atan(y/x) + PI;
    } else {
      curPhi = atan(y/x) - PI;
    }
    
  } else if(x==0) {
    
    if(y>=0) {
      curPhi = PI/2;
    } else {
      curPhi = -PI/2;
    }
    
  }


  return curPhi;
}


