/* === Doxygen Comment ======================================= */
/*! 
 * \file          test_wavelets.cxx
 * \copyright     (c) 2003 SNIFS-Supernova Factory Experiment
 * \date          Wed Aug  6 17:44:35 2003
 * \author        Emmanuel Gangler <e.gangler@ipnl.in2p3.fr>
 * \version       0.0
 * \brief         
 *                
 * $Id$
 **/
/* =========================================================== */

/* ----- IFU includes ------------------------------ */
#include <IFU_io.h>
#include <IFU_math.h>
#include "gsl/gsl_wavelet.h"
#include "gsl/gsl_statistics.h"
#include "gsl/gsl_sort.h"

#include "TH1.h"
#include "TFile.h"

#include "utils.h"

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

void readdata(double* data,char* name ,int nmax) {
  int i;
  FILE *f = fopen (name, "r");
  for (i = 0; i < nmax; i++)
    {
      fscanf (f, "%lg", &data[i]);
    }
  fclose (f);
}

void writedata(double* data, char* name, int nmax) {
  int i;
  FILE *f = fopen (name, "w");
  for (i = 0; i < nmax; i++)
    {
      fprintf (f, "%lg", data[i]);
    }
  fclose (f);
}

void Data2Histo(double* data, int n, char* histName) {

  TH1F * hist = new TH1F(histName,histName,n,-0.5, n-0.5);
  for (int i=0;i<n;i++) {
    hist->SetBinContent(i+1,data[i]);
  }
  hist->Write();
  delete hist;

}

void DataSquareHisto(double* data, int n, char* histName,char* IntegName) {

  double datasq[n];
  for (int i=0;i<n;i++) {
    datasq[i]=data[i]*data[i];
  }
  ut_sort_ascending(datasq,n);
  Data2Histo(datasq,n,histName);

  // integral histogram
  for (int i=1;i<n;i++){
    datasq[i]+=datasq[i-1];
  }
  
  Data2Histo(datasq,n,IntegName);
  
}

void DataSquareCut(double* data, int n, double efrac, double* cut) {

  double datasq[n];
  for (int i=0;i<n;i++) {
    datasq[i]=data[i]*data[i];
  }
  ut_sort_ascending(datasq,n);

  // integral - constant term
  double sum=0;
  for (int i=0;i<n-1;i++){
    sum+=datasq[i];
  }
   
  // look for the threshold
  double integ=0;
  for (int i=0;i<n-1;i++){
    integ+=datasq[i];
    if (integ/sum > efrac) {
      *cut = sqrt(datasq[i]);
      return;
    }
  }
}

void DataThreshold(double* data, int n, double cut) {

  for (int i=0;i<n-1;i++) {
    if (fabs(data[i])<cut)
      data[i]=0;
  }
}

void DataHisto(double* data, int n, char* histName) {
  int i;
  double low=ut_big_value,up=-ut_big_value;
  
  for (i=0;i<n;i++) {
    if (data[i] < low)
      low=data[i];
    if (data[i] > up)
      up=data[i];
  }

  double mean = ut_mean(data,n);
  double sigma = gsl_stats_sd_m(data,1,n,mean);
  double width = 3.49 * pow(n,-1/3) * sigma;
  width *=0.04;
  int nbin = int( (up-low)/width+3 );

  TH1F * hist = new TH1F(histName,histName,nbin,low-2*width,up+2*width);
  for (int i=0;i<n;i++) {
    hist->Fill(data[i]);
  }

  hist->Write();
  delete hist;
}


/*------------ main ----------------------*/
int main(int argc, char **argv) {

  int  n = 256;
  double *data = (double* ) malloc (n * sizeof (double));
  double *dataout = (double*) malloc (n * sizeof (double));
  //  double *abscoeff = malloc (n * sizeof (double));
  //size_t *p = malloc (n * sizeof (size_t));

  TFile *out;

  char **argval, **arglabel;
  set_arglist("-out none -efrac 0.01");
  init_session(argv,argc,&arglabel,&argval);

  out = new TFile(argval[0],"RECREATE");
  double efrac;
  get_argval(1,"%lf", &efrac);

  for (int i = 0; i < n; i++)    {
      data[i]=1+sin(i*1.0/20);
    }

  Data2Histo(data,n,"RawData");

  for (int level=1;level<10;level++) {
    wavelettransform(data,dataout,n,level);
    char name[80];
    sprintf(name,"Level%d",level);
    Data2Histo(dataout,n,name);
    sprintf(name,"Coefs%d",level);
    DataHisto(dataout,n,name);
    
    char iname[80];
    sprintf(name,"Square%d",level);
    sprintf(iname,"Integral%d",level);
    DataSquareHisto(dataout,n,name,iname);

    double cut;
    DataSquareCut(dataout,n,efrac,&cut);
    DataThreshold(dataout,n,cut);

    sprintf(name,"TSquare%d",level);
    sprintf(iname,"TIntegral%d",level);
    DataSquareHisto(dataout,n,name,iname);

    waveletinversetransform(dataout,dataout,n,level);
    sprintf(name,"Rec%d",level);
    Data2Histo(dataout,n,name);

  }
  
  //    gsl_wavelet_transform_inverse (w, data, 1, n, work);

  out->Purge();
  out->Close();
  delete out;
}

