/* === Doxygen Comment ======================================= */
/*! 
 * \file          rootify.cxx
 * \copyright     (c) 2003 SNIFS-Supernova Factory Experiment
 * \date          Wed Aug  6 17:44:35 2003
 * \author        Emmanuel Gangler <e.gangler@ipnl.in2p3.fr>
 * \version       0.0
 * \brief         
 *                
 * $Id$
 **/
/* =========================================================== */

/* ----- lib includes ------------------------------ */
#include "gsl/gsl_wavelet.h"

/* ----- image includes ------------------------------ */
#include "imagesnifs.hxx"
#include "catorfile.hxx"
//#include "section.hxx"
//#include "analyser.hxx"


/* ----- wavelettransform ------------------------------ */
void wavelettransform(double* datain, double* dataout,int n, int level) {
  
  gsl_wavelet *w;
  if (level==1)
    w = gsl_wavelet_alloc (gsl_wavelet_haar_centered, 2);
  else
    w = gsl_wavelet_alloc (gsl_wavelet_daubechies_centered, level*2);

  gsl_wavelet_workspace *work = gsl_wavelet_workspace_alloc (n);
  for (int i=0;i<n;i++) {
    dataout[i]=datain[i];
  }
  
  gsl_wavelet_transform_forward (w, dataout, 1, n, work);  

  gsl_wavelet_free (w);
  gsl_wavelet_workspace_free(work);
  
}

void waveletinversetransform(double* datain, double* dataout,int n, int level) {
  
  gsl_wavelet *w;
  if (level==1)
    w = gsl_wavelet_alloc (gsl_wavelet_haar_centered, 2);
  else
    w = gsl_wavelet_alloc (gsl_wavelet_daubechies_centered, level*2);

  gsl_wavelet_workspace *work = gsl_wavelet_workspace_alloc (n);
  for (int i=0;i<n;i++) {
    dataout[i]=datain[i];
  }
  
  gsl_wavelet_transform_inverse (w, dataout, 1, n, work);  

  gsl_wavelet_free (w);
  gsl_wavelet_workspace_free(work);
  
}


/*------------ main ----------------------*/
int main(int argc, char **argv) {

  char **argval, **arglabel;
  ImageSnifs* in,*out;
  char inName[lg_name+1], outName[lg_name+1];
  int level;

  set_arglist("-in none -out none -level 1 -threshold");
  init_session(argv,argc,&arglabel,&argval);

  CatOrFile inCat(argval[0]);
  CatOrFile outCat(argval[1]);
  get_argval(2,"%d", &level);

  /* loop on catalog */
  while(inCat.NextFile(inName) && outCat.NextFile(outName)) {
    printf("Opening now %s\n",inName);
    in = new ImageSnifs(inName);
    out = new ImageSnifs(*in,outName);

    int y=in->Ny(),maxY=1;
    // only power of 2 ...
    while ((y=y>>1)) 
      maxY=maxY<<1;
    double *data = new double[maxY];

    // initialize wavelet
    gsl_wavelet *w;
    if (level==1)
      w = gsl_wavelet_alloc (gsl_wavelet_haar_centered, 2);
    else
      w = gsl_wavelet_alloc (gsl_wavelet_daubechies_centered, level*2);
    gsl_wavelet_workspace *work = gsl_wavelet_workspace_alloc (maxY);

    // transform image
    for (int i=0;i<in->Nx();i++) {
      for (int j=0;j<maxY;j++) {
        data[j]=in->RdFrame(i,j);
      }
      gsl_wavelet_transform_forward (w, data, 1, maxY, work);  
      if (is_set(argval[3])) {
        int j=0;
        int l=3500;
        int l2=l/2+2048;
        int l4=l/4+1024;
        int l8=l/4+2048;
                for (;j<1024;j++) {
                    if (fabs(data[j])<250) data[j]=0;
         }
        for (;j<maxY;j++) {
          if (fabs(data[j])<2000) data[j]=0;
        }
        gsl_wavelet_transform_inverse (w, data, 1, maxY, work);
      }
      
      for (int j=0;j<maxY;j++) {
        out->WrFrame(i,j,data[j]);
      }
    }
    delete[] data;
    gsl_wavelet_free (w);
    gsl_wavelet_workspace_free(work);
     
    delete in;
    delete out;
  }
  
  exit_session(0);

  return(0);
}



