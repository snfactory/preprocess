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
#define GSL_DISABLE_DEPRECATED
#include "gsl/gsl_wavelet.h"
#include "gsl/gsl_wavelet2d.h"

/* ----- image includes ------------------------------ */
#include "imagesnifs.hxx"
#include "catorfile.hxx"
#include "section.hxx"
#include "analyser.hxx"


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

void dwt2d(ImageSnifs* in, ImageSnifs* out,int level,double cut){
// needs a square matrix we don't have...
// embed it in a larger square matrix
  int maxX=tomaxpow2(in->Nx());
  if (maxX<in->Nx())
    maxX*=2;
  int maxY=tomaxpow2(in->Ny());
  if (maxY<in->Ny())
    maxY*=2;
  if (maxX<maxY)
    maxX=maxY;
  else
    maxY=maxX;
  int X0=(maxX-in->Nx())/2;
  int Y0=(maxY-in->Ny())/2;


  // Build image
  double *data = new double [maxX*maxY];

    for (int i=0;i<in->Nx();i++) {
      for (int j=0;j<in->Ny();j++) {
        data[(i+X0)*maxX + j+Y0]=in->RdFrame(i,j);
      }
    }
    for (int i=0;i<X0;i++) {
      for (int j=0;j<in->Ny();j++) {
        data[(X0-i-1)*maxX + j+Y0]=data[(i+X0)*maxX + j+Y0];
      }
    }
    for (int i=in->Nx()+X0;i<maxX;i++) {
      for (int j=0;j<in->Ny();j++) {
        data[i*maxX + j+Y0]=data[((in->Nx()+X0)*2-i-1)*maxX + j+Y0];
      }
    }
    for (int i=0;i<maxX;i++) {
      for (int j=0;j<Y0;j++) {
        data[i*maxX + Y0-j-1]=data[i*maxX + j+Y0];
      }
    }
    for (int i=0;i<maxX;i++) {
      for (int j=in->Ny()+Y0;j<maxY;j++) {
        data[i*maxX + j]=data[i*maxX + (in->Ny()+Y0)*2-j-1];
      }
    }

    // initialize wavelet
    gsl_wavelet *w;
    if (level==1)
      w = gsl_wavelet_alloc (gsl_wavelet_haar_centered, 2);
    else
      w = gsl_wavelet_alloc (gsl_wavelet_daubechies_centered, level*2);
    gsl_wavelet_workspace *work = gsl_wavelet_workspace_alloc (maxY);

    gsl_wavelet2d_transform_forward (w, data, maxX, maxX,maxY, work);  

    if (cut>0) {
      int nremoved=0;
      for (int i=0;i<maxX;i++) 
        for (int j=0;j<maxY;j++) 
          if (fabs(data[i*maxX+j])<cut) {
            nremoved++;
            data[i*maxX+j]=0;
          }
      printf("removed %f%% of data\n",nremoved*1.0/maxX/maxY*100);
    }
    
    gsl_wavelet2d_transform_inverse (w, data, maxX, maxX,maxY, work);  
    

        //out->DeleteFrame();
        //out->CreateFrame("image",maxX,maxY);
    for (int i=0;i<in->Nx();i++)
      for (int j=0;j<in->Ny();j++) {
        out->WrFrame(i,j,data[(i+X0)*maxX + (j+Y0)]);
      }



    gsl_wavelet_free (w);
    gsl_wavelet_workspace_free(work);
    //dwt_horizontal(in,out,level);
    //dwt_vertical(out,out,level);}
    
}

double dwtThreshold(ImageSnifs* out) {
  Section S(1,out->Nx(),1,out->Ny());
  
  if (out->Variance()) {
    ImageAnalyser ana(out->Variance(),&S);
    return ana.Quantile(0.5);
  }
  return 0;
}


/*------------ main ----------------------2*/
int main(int argc, char **argv) {

  char **argval, **arglabel;
  ImageSnifs* in,*out;
  char inName[lg_name+1], outName[lg_name+1];
  int level;

  set_arglist("-in none -out none -level 1 -threshold 0 -d2|hori|vert");
  init_session(argv,argc,&arglabel,&argval);

  CatOrFile inCat(argval[0]);
  CatOrFile outCat(argval[1]);
  get_argval(2,"%d", &level);

  double cut;
  get_argval(3,"%lf",&cut);

  /* loop on catalog */
  while(inCat.NextFile(inName) && outCat.NextFile(outName)) {
    printf("Opening now %s\n",inName);
    in = new ImageSnifs(inName);
    out = new ImageSnifs(*in,outName,0,1);
    double variance = dwtThreshold(out);
    printf("Variance estimated %f \n",variance);

    if (!strcmp(arglabel[4],"-d2"))
      dwt2d(in, out,level,cut*sqrt(variance));
    else if (!strcmp(arglabel[4],"-hori"))
      dwt_horizontal(in, out,level);
    else if (!strcmp(arglabel[4],"-vert"))
      dwt_vertical(in, out,level);
    else 
      print_error("You must supply one of -d2 -hori or -vert option");

    delete in;
    delete out;
  }
  
  exit_session(0);

  return(0);
}



