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
#include "gsl/gsl_wavelet2d.h"

/* ----- image includes ------------------------------ */
#include "imagesnifs.hxx"
#include "catorfile.hxx"
#include "section.hxx"
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

/*---------- To be not lost ------------------------------*/
#ifdef OLD_CODE
// This is taken fomr somewhere in the main
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
#endif

int tomaxpow2(int x) {
    int maxX=1;
    // only power of 2 ...
    while ((x=x>>1)) 
      maxX=maxX<<1;
    return maxX;
}

void dwt_horizontal(ImageSnifs* in, ImageSnifs* out,int level){
    int maxX=tomaxpow2(in->Nx());

    double *data = new double[maxX];

    // initialize wavelet
    gsl_wavelet *w;
    if (level==1)
      w = gsl_wavelet_alloc (gsl_wavelet_haar_centered, 2);
    else
      w = gsl_wavelet_alloc (gsl_wavelet_daubechies_centered, level*2);
    gsl_wavelet_workspace *work = gsl_wavelet_workspace_alloc (maxX);

    // transform image
    for (int j=0;j<in->Ny();j++) {
      for (int i=0;i<maxX;i++) {
        data[i]=in->RdFrame(i,j);
      }
    
      gsl_wavelet_transform_forward (w, data, 1, maxX, work);  
      
      for (int i=0;i<maxX;i++) {
        out->WrFrame(i,j,data[i]);
      }
    }

    delete[] data;
    gsl_wavelet_free (w);
    gsl_wavelet_workspace_free(work);
}

void dwt_vertical(ImageSnifs* in, ImageSnifs* out,int level){
    int maxY=tomaxpow2(in->Ny());

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
      
      for (int j=0;j<in->Ny();j++) {
        out->WrFrame(i,j,data[j]);
      }
    }

    delete[] data;
    gsl_wavelet_free (w);
    gsl_wavelet_workspace_free(work);
}

void dwt2d(ImageSnifs* in, ImageSnifs* out,int level){
// needs a square matrix we don't have...
// order doesn't matter :-D
dwt_horizontal(in,out,level);
dwt_vertical(out,out,level);
}

void dwtThreshold(ImageSnifs* out) {
  Section S(1,out->Nx(),1,out->Ny());
  
  double level;
  //if (out->Variance()) {
  //double rms;
  // int nout;
    
    // wait for the new GSL ...

    //ImageAnalyser ana(out->Variance(),&S);
    //ana.SigmaClippedInfo(10,&level,&rms,&nout);
    //level *=3;
  //}
  //else 
    print_warning("dwtThreshold : can't compute the threshold\n");
  int nout=out->AbsThreshold(level);
  printf("Removed %d values for level %f\n",nout,level);
}


/*------------ main ----------------------2*/
int main(int argc, char **argv) {

  char **argval, **arglabel;
  ImageSnifs* in,*out;
  char inName[lg_name+1], outName[lg_name+1];
  int level;

  set_arglist("-in none -out none -level 1 -threshold -d2|hori|vert");
  init_session(argv,argc,&arglabel,&argval);

  CatOrFile inCat(argval[0]);
  CatOrFile outCat(argval[1]);
  get_argval(2,"%d", &level);

  /* loop on catalog */
  while(inCat.NextFile(inName) && outCat.NextFile(outName)) {
    printf("Opening now %s\n",inName);
    in = new ImageSnifs(inName);
    out = new ImageSnifs(*in,outName,0,1);

    if (!strcmp(arglabel[4],"-d2"))
      dwt2d(in, out,level);
    else if (!strcmp(arglabel[4],"-hori"))
      dwt_horizontal(in, out,level);
    else if (!strcmp(arglabel[4],"-vert"))
      dwt_vertical(in, out,level);
    else 
      print_error("You must supply one of -2d -hori or -vert option");

    if (is_true(argval[3]))
      dwtThreshold(out);
      
    delete in;
    delete out;
  }
  
  exit_session(0);

  return(0);
}



