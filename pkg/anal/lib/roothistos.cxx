
#include "gsl/gsl_statistics.h"

#include "image.hxx"
#include "section.hxx"
#include "utils.h"
#include "roothistos.hxx"

/*--- gsl routines ---*/
#include "gsl/gsl_fft_real.h"
#include "gsl/gsl_fft_halfcomplex.h"

/*--- ROOT includes ---*/
#include "TH1F.h"
#include "TProfile.h"
#include "TMatrix.h"

/* --------------- HistoSetMinMax ------------------- */
void RootAnalyser::HistoSetMinMax(TH1 * hist) {
  double min = hist->GetMinimum();
  double max = hist->GetMaximum();
  double center = (max+min)/2;
  hist->SetMaximum(center+(max-center)*1.1);
  hist->SetMinimum(center-(max-center)*1.1);
}


/* --------------- FillHistoLine ------------------- */
void RootAnalyser::FillHistoLine() {

  int npix = fSec->XLength()*fSec->YLength();
  char histName[lg_name+1];
  sprintf(histName,"Line%s",fSec->Name());
  TH1F * histo = new TH1F(histName,"Pixel values",npix,1-0.5,npix+0.5);
  
  for (int j=fSec->YFirst();j<fSec->YLast();j++)
    for (int i=fSec->XFirst();i<fSec->XLast();i++) {
      int ipix = i -fSec->XFirst() + fSec->XLength()*(j-fSec->YFirst());
      histo->SetBinContent(ipix+1,fImage->RdFrame(i,j));
    }
  HistoSetMinMax(histo);
  histo->Write();
  delete histo;
}

/* ----- HorizontalProfile ------------------------- */
void RootAnalyser::HorizontalProfile(float sigma) {
  
  char histName[lg_name+1];
  sprintf(histName,"HoriPrfl%s",fSec->Name());
  TProfile * prof = new TProfile(histName,"Horizontal Profile",fSec->XLength(),fSec->XFirst()-0.5,fSec->XLast()-0.5);
  
  int i,j,remain;
  double *over, *overptr;
  int nover = fSec->YLength();

  over = new double[nover];

  for (i=fSec->XFirst();i<fSec->XLast();i++) {
    for (j=fSec->YFirst();j<fSec->YLast();j++) 
      over[j-fSec->YFirst()] = fImage->RdFrame(i,j);

    /* purify the overscan from spurious */
    overptr=over;
    remain = nover;
    if (sigma>0) 
      ut_trunc_sigma_unknown(&overptr,&remain,sigma);
    
    for (j=0;j<remain;j++) {
      prof->Fill(i,overptr[j]);
    }
  }
  HistoSetMinMax(prof);
  prof->Write();
  delete prof;
  delete[]over;
  
}

/* ----- VerticalProfile ------------------------- */
void RootAnalyser::VerticalProfile(float sigma){
  
  char histName[lg_name+1];
  sprintf(histName,"VertPrfl%s",fSec->Name());
  TProfile * prof = new TProfile(histName,"Vertical Profile",fSec->YLength(),fSec->YFirst()-0.5,fSec->YLast()-0.5);
  
  int i,j,remain;
  double *over,*overptr;
  int nover = fSec->XLength();
  over = new double[nover];

  for (j=fSec->YFirst();j<fSec->YLast();j++) {
    for (i=fSec->XFirst();i<fSec->XLast();i++) 
      over[i-fSec->XFirst()] = fImage->RdFrame(i,j);

    /* purify the overscan from spurious */
    remain=nover;
    overptr=over;
    if (sigma>0) 
     ut_trunc_sigma_unknown(&overptr,&remain,sigma);
    
    for (i=0;i<remain;i++) {
      prof->Fill(j,overptr[i]);
    }
  }
  HistoSetMinMax(prof);
  prof->Write();
  delete prof;
  delete[] over; 
}


/*-------------- OddEvenVerticalProfile ----------------------------- */
void RootAnalyser::OddEvenVerticalProfile(float sigma){
  
  char histName[lg_name+1];
  sprintf(histName,"OeVertPrfl%s",fSec->Name());
  TH1D * prof = new TH1D(histName,"Odd-Even Vertical profile",fSec->YLength(),fSec->YFirst()-0.5,fSec->YLast()-0.5);
  
  int i,j,remainOdd,remainEven;
  double *odd,*oddptr;
  double *even,*evenptr;

  odd = new double[(fSec->XLength()+1)/2];
  even = new double[fSec->XLength()/2];

  for (j=fSec->YFirst();j<fSec->YLast();j++) {
    for (i=fSec->XFirst();i<fSec->XLast();i+=2) { /* 'odd' */ 
      odd[(i-fSec->XFirst())/2] = fImage->RdFrame(i,j);
      if (i+1<fSec->XLast()) {
        even[(i-fSec->XFirst())/2] = fImage->RdFrame(i+1,j);
      }
    }
  
    /* purify from spurious */
    remainOdd=(fSec->XLength()+1)/2;
    remainEven=(fSec->XLength())/2;
    oddptr=odd;
    evenptr=even;

    if (sigma>0) {
      ut_trunc_sigma_unknown(&oddptr,&remainOdd,sigma);
      ut_trunc_sigma_unknown(&evenptr,&remainEven,sigma);
    }

    /* fill the profile */
    double oddMean=ut_mean(oddptr,remainOdd);
    double evenMean=ut_mean(evenptr,remainEven);
    double err = sqrt(gsl_stats_variance_m(oddptr,1,remainOdd,oddMean)/remainOdd + gsl_stats_variance_m(evenptr,1,remainEven,evenMean)/remainEven);
    int bin = prof->GetXaxis()->FindBin(j);
    prof->SetBinContent(bin,oddMean-evenMean);
    prof->SetBinError(bin,err);
  }
  HistoSetMinMax(prof);
  prof->Write();
  delete prof;
  delete[] odd;
  delete[] even;
}

/* ----------- OverscanError ---------------- */
void RootAnalyser::OverscanError(float sigma){
  
  int i,j,remain;
  double *over,*overptr; 
  int nover = fSec->XLength();
  float offset,var;

  over = new double[nover];

  char histName[lg_name+1];
  sprintf(histName,"OvscErr%s",fSec->Name());
  TH1F * histo = new TH1F(histName,"Error from Oversvcan",100,0,0);

  /* then fill the histos*/
  for (j=fSec->YFirst();j<fSec->YLast();j++) {
    for (i=fSec->XFirst();i<fSec->XLast();i++) 
      over[i-fSec->XFirst()] = fImage->RdFrame(i,j);

    
    /* purify the overscan from spurious */
    remain=nover;
    overptr=over;
    if (sigma>0) 
      ut_trunc_sigma_unknown(&overptr,&remain,sigma);

    offset = ut_mean(overptr,remain);
    var=gsl_stats_variance_m(overptr,1,remain,offset);

    for (i=0;i<remain;i++) {
      float val = (overptr[i]-offset)*sqrt((remain+1.0)/remain);
      histo->Fill(val);
    }
  }
  
  histo->Write();
  
  delete histo;
  delete[] over;
}

/*-------------- HistoData ---------------------------- */
void RootAnalyser::HistoData(){

  double min,max,center;
  fImage->MinMax(fSec,&min, &max);
  center = (min+max)/2;
  max = center + (max-center)*1.01;
  min = center - (max-center)*1.01;

  char histName[lg_name+1];
  sprintf(histName,"Histo%s",fSec->Name());
  TH1F * histo = new TH1F(histName,"Data values",1000,min,max);

  for (int iy = fSec->YFirst(); iy<fSec->YLast();iy++)
    for (int ix = fSec->XFirst(); ix<fSec->XLast();ix++)
      histo->Fill(fImage->RdFrame(ix,iy));
    
  histo->Write();
  delete histo;
}

/*-------------- HistoData ---------------------------- */
void RootAnalyser::MatrixData(){

  double min,max;
  fImage->MinMax(fSec,&min, &max);

  TMatrix * mat = new TMatrix(fSec->XLength(),fSec->YLength());
  for (int iy = fSec->YFirst(); iy<fSec->YLast();iy++)
    for (int ix = fSec->XFirst(); ix<fSec->XLast();ix++)
      (*mat)(ix-fSec->XFirst(),iy-fSec->YFirst())= fImage->RdFrame(ix,iy);

  char histName[lg_name+1];
  sprintf(histName,"Matrix%s",fSec->Name());
  mat->Write(histName);
  delete mat;
}

/*-------------- hf ------------------------------------------- */
void RootAnalyser::HighFrequency()
{
  double rapporth,rapportv,minh,minv,mintot,maxh,maxv,maxtot;

  minh = ut_big_value;
  minv = ut_big_value;
  maxh = -ut_big_value;
  maxv = -ut_big_value;

  // first compute histogram bounds
  for (int j=fSec->YFirst(); j<fSec->YLast();j++)
    for (int i=fSec->XFirst(); i<fSec->XLast();i++){
      if (j<fSec->YLast()-1) {
        rapporth = fImage->RdFrame(i,j)/fImage->RdFrame(i,j+1);
        if (rapporth<minh) 
          minh = rapporth;
        else if (rapporth>maxh) 
          maxh = rapporth; 
      }
      if (i<fSec->XLast()-1) {
        rapportv = fImage->RdFrame(i,j)/fImage->RdFrame(i+1,j);
        if (rapportv<minv) 
          minv = rapportv;
        else if (rapportv>maxv) 
          maxv = rapportv;
      }
    }
 
  if (minh<minv){
    mintot = minh;
  }
  else mintot = minv;
  if (maxh<maxv){
    maxtot = maxv;
  }
  else maxtot = maxh;

  char histName[lg_name+1];
  sprintf(histName,"HfRatioH%s",fSec->Name());
  TH1F* histoh = new TH1F(histName,"Horizontal ratio",1000,minh-1,maxh+1);
  sprintf(histName,"HfRatioV%s",fSec->Name());
  TH1F* histov = new TH1F(histName,"Vertical ratio",1000,minv-1,maxv+1);
  sprintf(histName,"HfRatio%s",fSec->Name());
  TH1F* histotot = new TH1F(histName,"Neighbours Ratio",1000,mintot-1,maxtot+1);
    

  for (int j=fSec->YFirst(); j<fSec->YLast();j++){
    for (int i=fSec->XFirst(); i<fSec->XLast();i++){
      if (j<fSec->YLast()-1) {
        rapporth = fImage->RdFrame(i,j)/fImage->RdFrame(i,j+1);
        histoh->Fill(rapporth);
        histotot->Fill(rapporth);
      }
      if (i<fSec->XLast()-1) {
        rapportv = fImage->RdFrame(i,j)/fImage->RdFrame(i+1,j);
        histov->Fill(rapportv);
        histotot->Fill(rapportv);
      }
    }
  }
  
  histoh->Write();
  histov->Write();
  histotot->Write();
  
  delete histoh;
  delete histov;
  delete histotot;
  
}

/*-------------- fft-it ------------------*/
void RootAnalyser::Fft()
{

  int length = fSec->XLength();
  int nlines=fSec->YLength();
  
  double *line = new double[length];
  double *fft = new double[length/2+1];
  double *autocorr = new double[length];
  double *autocorrwork = new double[length];
  double *powerspec = new double[length/2+1];
  double *powerspec2 = new double[length/2+1];

  for (int i=0;i<length/2+1;i++) {
    fft[i]=0;
    powerspec2[i] = 0;
    autocorr[i]=0;
  }
  for (int i=length/2+1;i<length;i++) {
    autocorr[i]=0;
  }
  

  //  double seuils[6]={-0.1,0.2,0.4,0.6,0.8,0.99};
  //int nseuil=0, maxseuil=6;

  char histName[lg_name+1];
  sprintf(histName,"Fft%s",fSec->Name());
  TH1D *hfft = new TH1D(histName,"Sum of FFTs",length/2+1,-0.5/length, (length/2)*1.0/length + 0.5/length);
  sprintf(histName,"AutoCov%s",fSec->Name());
  TH1D *hauto = new TH1D(histName,"Autocovariance",length,0,length);
  sprintf(histName,"Power%s",fSec->Name());
  TH1D *hpower = new TH1D(histName,"Power spectrum",length/2+1,-0.5/length, (length/2)*1.0/length + 0.5/length);
  sprintf(histName,"Power2%s",fSec->Name());
  TH1D *hpower2 = new TH1D(histName,"Power spectrum 2",length/2+1,-0.5/length, (length/2)*1.0/length + 0.5/length);

  gsl_fft_real_wavetable * real;
  gsl_fft_real_workspace * work;

  work = gsl_fft_real_workspace_alloc (length);
  real = gsl_fft_real_wavetable_alloc (length);

  for (int j=fSec->YFirst();j<fSec->YLast();j++) {
    if (!((j-fSec->YFirst())%1000))
      printf("Fft : processing line %d\n",j);
    
    /* fill a line */
    for (int i=fSec->XFirst();i<fSec->XLast();i++)
      line[i-fSec->XFirst()] =fImage->RdFrame(i,j);

    /* take its autocorrelation function now since fft destroys data*/
    ut_autocovariance(line,autocorrwork,length);
    for (int i=0;i<length;i++)
      autocorr[i] += autocorrwork[i]/nlines;

    /* take its fft */
    gsl_fft_real_transform (line, 1, length, real, work);

    /* extract the direct fft in some cases */
    //TH1D* hffttmp=0;    
    //if (nseuil<maxseuil && (j-yin)*1.0/(yout-yin) > seuils[nseuil]) {
    //  nseuil++;
    //  char name[80],title[80];
    //  sprintf (name,"fft%d",j);
    //  sprintf (title,"Fft of line %d",j);
    //  hffttmp = new TH1D(name,title,length/2+1,0,length/2+1);
    //  }

    /* increment fft knowledge */
    fft[0] += fabs(line[0])/nlines;
    powerspec2[0] += line[0]*line[0]/nlines;
    //    if (hffttmp)
    //  hffttmp->SetBinContent(0,fabs(line[0]));
    int i=1;
    for (;i<length-1;i+=2) {
      double module2 = line[i]*line[i]+line[i+1]*line[i+1];
      fft[i/2+1] += sqrt(module2)/nlines;
      powerspec2[i/2+1] += module2/nlines;
      //      if (hffttmp)
      //  hffttmp->SetBinContent(i/2+1,module);
    }
    if (i<length) {
      fft[i/2+1] += fabs(line[i])/nlines;    
      powerspec2[i/2+1] += line[i]*line[i]/nlines;    
      //      if (hffttmp)
      //  hffttmp->SetBinContent(i/2+1,fabs(line[i]));
    }

    //    if (hffttmp) {
    //  hffttmp->Write();
    //  delete hffttmp;
    //}
    
  }

  for (int i=0;i<length/2+1;i++) {
    hfft->SetBinContent(i+1,fft[i]);
    hpower2->SetBinContent(i+1,powerspec2[i]);
  }
  
  for (int i=0;i<length;i++)
    hauto->SetBinContent(i+1,autocorr[i]);
  

  /* take the fft of autocorrelation function */
  gsl_fft_real_transform (autocorr, 1, length, real, work);
  gsl_fft_real_wavetable_free (real);
  gsl_fft_real_workspace_free (work);

  powerspec[0]=fabs(autocorr[0]);
  int i=1;
  for (;i<length-1;i+=2) {
    powerspec[i/2+1]=sqrt(autocorr[i]*autocorr[i] 
                      + autocorr[i+1]*autocorr[i+1]);
  }
  if (i<length)
    powerspec[i/2+1]=fabs(autocorr[i]);
  
  for (int i=0;i<length/2+1;i++)
    hpower->SetBinContent(i+1,powerspec[i]);


  hfft->Write();
  hauto->Write();
  hpower->Write();
  hpower2->Write();
  delete hfft;
  delete hauto;
  delete hpower;
  delete[] line;
  delete[] fft;
  delete[] autocorr;
  delete[] autocorrwork;
  delete[] powerspec;
  delete[] powerspec2;
  
}

