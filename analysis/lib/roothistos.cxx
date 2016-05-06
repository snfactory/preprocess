
#include <vector>

#include "gsl/gsl_statistics.h"

#include "image.hxx"
#include "section.hxx"
#include "utils.h"
#include "kombinator.hxx"
#include "roothistos.hxx"

/*--- gsl routines ---*/
#include "gsl/gsl_fft_real.h"
#include "gsl/gsl_fft_halfcomplex.h"

/*--- ROOT includes ---*/
#include "TH1F.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TMatrix.h"

/* --------------- Constructor ------------------- */
RootAnalyser::RootAnalyser(ImageSimple* I, Section* S) {
  SetImage(I);
  SetSection(S);
}

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
void RootAnalyser::HorizontalProfile(float sigma, ImageSnifs* selection) {
  
  char histName[lg_name+1];
  sprintf(histName,"HoriPrfl%s",fSec->Name());
  TH1D*  prof = new TH1D(histName,"Horizontal Profile",fSec->XLength(),fSec->X1()-0.5,fSec->X2()+0.5);
  
  int i,j,remain;
  double *over, *overptr;
  int nover = fSec->YLength();
  over = new double[nover];
  
    
  vector<double> vals;
  vector<double> vars;
  KGauss K=KGauss(sigma);
  K.SetBrute(1);
  double mean,var,rms;

  for (i=fSec->XFirst();i<fSec->XLast();i++) {
    if (fImage->Variance() && sigma>0) {
      vals.clear();
      vars.clear();
      for (j=fSec->YFirst();j<fSec->YLast();j++) 
        if ( ! selection || selection->RdFrame(i,j)==1) {
          vals.push_back(fImage->RdFrame(i,j));
          vars.push_back(fImage->Variance()->RdFrame(i,j));
        }
      
      K.Kombine(&vals,&vars,&mean,&var);
      rms=sqrt(var);
    } else {
      int iover=0;
      for (j=fSec->YFirst();j<fSec->YLast();j++) 
        if ( ! selection || selection->RdFrame(i,j)==1) {
          over[iover++] = fImage->RdFrame(i,j);
        }
      overptr=over;
      remain = iover;
      if (sigma>0) 
        ut_trunc_sigma_unknown(&overptr,&remain,sigma);
      mean=gsl_stats_mean(overptr,1,remain);
      rms=gsl_stats_sd_m(overptr,1,remain,mean);
    }
    prof->SetBinContent(i+1-fSec->XFirst(),mean);
    prof->SetBinError(i+1-fSec->XFirst(),rms); 
  }
  HistoSetMinMax(prof);
  prof->Write();
  delete prof;
  delete[]over;
  
}

/* ----- HorizontalProfile ------------------------- */
void RootAnalyser::HorizontalMode( ImageSnifs* selection) {
  
  char histName[lg_name+1];
  sprintf(histName,"HoriMode%s",fSec->Name());
  TH1D * hmode = new TH1D(histName,"Horizontal Profile",fSec->XLength(),fSec->X1()-0.5,fSec->X2()+0.5);
  
  int i,j;
  double *over;
  int nover = fSec->YLength();

  over = new double[nover];

  for (i=fSec->XFirst();i<fSec->XLast();i++) {
    int iover=0;
    for (j=fSec->YFirst();j<fSec->YLast();j++) 
      if ( ! selection || selection->RdFrame(i,j)==1)
        over[iover++] = fImage->RdFrame(i,j);

    if (iover>0)
    /* purify the overscan from spurious */
    hmode->SetBinContent(i-fSec->XFirst()+1,ut_mode2(over,iover));

  }
  HistoSetMinMax(hmode);
  hmode->Write();
  delete hmode;
  delete[]over;
  
}

/* ----- VerticalProfile ------------------------- */
void RootAnalyser::VerticalProfile(float sigma,ImageSnifs* selection){
  
  char histName[lg_name+1];
  sprintf(histName,"VertPrfl%s",fSec->Name());
  TProfile * prof = new TProfile(histName,"Vertical Profile",fSec->YLength(),fSec->Y1()-0.5,fSec->Y2()+0.5);
  
  int i,j,remain;
  double *over,*overptr;
  int nover = fSec->XLength();
  over = new double[nover];

  for (j=fSec->YFirst();j<fSec->YLast();j++) {
    int iover=0;
    for (i=fSec->XFirst();i<fSec->XLast();i++) 
      if ( ! selection || selection->RdFrame(i,j)==1)
        over[iover++] = fImage->RdFrame(i,j);

    /* purify the overscan from spurious */
    remain=iover;
    overptr=over;
    if (sigma>0) 
     ut_trunc_sigma_unknown(&overptr,&remain,sigma);
    
    for (i=0;i<remain;i++) {
      prof->Fill(j+1,overptr[i]);
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
  TH1D * prof = new TH1D(histName,"Odd-Even Vertical profile",fSec->YLength(),fSec->Y1()-0.5,fSec->Y2()-0.5);
  
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
void RootAnalyser::HistoData(double start, double end,int nbins, ImageSnifs* selection){

  char histName[lg_name+1];
  sprintf(histName,"Histo%s",fSec->Name());
  TH1F * histo = HistoDataBuild(histName,nbins,start,end,selection);

  histo->Write();
  delete histo;
}

/*-------------- HistoDataBuild ---------------------------- */
TH1F* RootAnalyser::HistoDataBuild(char* histName, int nbins,double start, double end,ImageSnifs* selection){

  double min=ut_big_value,max=-ut_big_value,center;
  if (start==end) {
  for (int iy = fSec->YFirst(); iy<fSec->YLast();iy++)
    for (int ix = fSec->XFirst(); ix<fSec->XLast();ix++)
      if ( ! selection || selection->RdFrame(ix,iy)==1) {
        if (fImage->RdFrame(ix,iy)<min)
          min=fImage->RdFrame(ix,iy);
        if (fImage->RdFrame(ix,iy)>max)
          max=fImage->RdFrame(ix,iy);
        center = (min+max)/2;
      }
    max = center + (max-center)*1.01;
    min = center - (max-center)*1.01;
  } else {
    max=end;
    min=start;
  }
  

  TH1F * histo = new TH1F(histName,"Data values",nbins,min,max);

  for (int iy = fSec->YFirst(); iy<fSec->YLast();iy++)
    for (int ix = fSec->XFirst(); ix<fSec->XLast();ix++)
      if ( ! selection || selection->RdFrame(ix,iy)==1)
        histo->Fill(fImage->RdFrame(ix,iy));
    
  return histo;
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
  TH1F* histoh = new TH1F(histName,"Horizontal ratio",10000,0.5,2);
  sprintf(histName,"HfRatioV%s",fSec->Name());
  TH1F* histov = new TH1F(histName,"Vertical ratio",10000,0.5,2);
  sprintf(histName,"HfRatio%s",fSec->Name());
  TH1F* histotot = new TH1F(histName,"Neighbours Ratio",10000,0.5,2);
  sprintf(histName,"HfDiff%s",fSec->Name());
  TH1F* histodiff = new TH1F(histName,"Neighbours diff",10000,-5000,5000);
    

  for (int j=fSec->YFirst(); j<fSec->YLast();j++){
    for (int i=fSec->XFirst(); i<fSec->XLast();i++){
      if (j<fSec->YLast()-1) {
        rapporth = fImage->RdFrame(i,j)/fImage->RdFrame(i,j+1);
        histoh->Fill(rapporth);
        histotot->Fill(rapporth);
        histodiff->Fill(fImage->RdFrame(i,j)-fImage->RdFrame(i,j+1));
      }
      if (i<fSec->XLast()-1) {
        rapportv = fImage->RdFrame(i,j)/fImage->RdFrame(i+1,j);
        histov->Fill(rapportv);
        histotot->Fill(rapportv);
        histodiff->Fill(fImage->RdFrame(i,j)-fImage->RdFrame(i+1,j));
      }
    }
  }
  
  histoh->Write();
  histov->Write();
  histotot->Write();
  histodiff->Write();

  delete histoh;
  delete histov;
  delete histotot;
  delete histodiff;
  
}

/*-------------- hf ------------------------------------------- */
void RootAnalyser::VertHighFrequencyProf(double min, double max, int nbin)
{

  char histName[lg_name+1];
  sprintf(histName,"HfVDiff%s",fSec->Name());
  TProfile* prof = new TProfile(histName,"Vertical Diff vs. value",nbin,min,max);

  for (int j=fSec->YFirst(); j<fSec->YLast();j++){
    for (int i=fSec->XFirst(); i<fSec->XLast();i++){
      if (i<fSec->XLast()-1) {
        double diffv = fImage->RdFrame(i,j) - fImage->RdFrame(i+1,j);
        double sumv = (fImage->RdFrame(i,j) + fImage->RdFrame(i+1,j))/2;
        prof->Fill(sumv,diffv);
      }
    }
  }
  
  prof->Write();
  delete prof;
  
}

/*-------------- impulse response ------------------------------------------- */
void RootAnalyser::ImpulseProfile(int minx, int maxx, int nbiny, double miny, double maxy)
{

  char histName[lg_name+1];
  int nbinx=maxx-minx+1;
  sprintf(histName,"Impulse%s",fSec->Name());
  TProfile2D* prof = new TProfile2D(histName,"Impulse response",nbinx,minx-0.5,maxx+0.5,nbiny,miny,maxy);
  int maxData=20;
  double data[nbinx][nbiny][maxData];
  int n[nbinx][nbiny];
  
  for(int i=0;i<nbinx;i++)
    for(int j=0;j<nbiny;j++)
      n[i][j]=0;

  for (int j=fSec->YFirst(); j<fSec->YLast();j++){
    for (int i=fSec->XFirst()-miny; i<fSec->XLast()-maxy-1;i++){
      // 1-D laplacian detection of impulse signals
      //double laplacian= fImage->RdFrame(i,j)- fImage->RdFrame(i-1,j)/2- fImage->RdFrame(i+1,j)/2;
      //double lapvar=fImage->Variance()->RdFrame(i,j)- fImage->RdFrame(i-1,j)/4- fImage->RdFrame(i+1,j)/4;
      if (fImage->RdFrame(i,j)>fImage->RdFrame(i+1,j)*2 && fImage->RdFrame(i,j)>fImage->RdFrame(i-1,j)*2 && fImage->RdFrame(i,j)>10 ){
        double y=log10(fImage->RdFrame(i,j));
        for (int ix=minx;ix<maxx+1;ix++) {
          int binx=prof->GetXaxis()->FindFixBin(ix)-1;
          int biny=prof->GetYaxis()->FindFixBin(y)-1;
          if (binx<0 || biny<0 || binx>=nbinx || biny>=nbiny)
            continue;
          if (n[binx][biny]==maxData) {
            double z=ut_median(data[binx][biny],n[binx][biny]);
            prof->Fill(ix,y,z);
            n[binx][biny]=0;
          }
          data[binx][biny][n[binx][biny]]=fImage->RdFrame(i+ix,j);
          n[binx][biny]++;
        }
      }
    }
  }
  /*
  for(int i=0;i<nbinx;i++)
    for(int j=0;j<nbiny;j++)
      if (n[i][j]>0) {
        double z=ut_median(data[i][j],n[i][j]);
        prof->Fill(prof->GetXaxis()->GetBinCenter(i+1),
                   prof->GetYaxis()->GetBinCenter(j+1),z);
      }
  */

  
  prof->Write();
  delete prof;
  
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

/*-------------- ADC Bits --------------------------- */
void RootAnalyser::ADCBits(int NBits){

  char histName[lg_name+1];
  sprintf(histName,"ADC%s",fSec->Name());
  TProfile * hADC = new TProfile(histName,"ADC Bits value",16,-0.5,15.5);
  
  for (int iy = fSec->YFirst(); iy<fSec->YLast();iy++)
    for (int ix = fSec->XFirst(); ix<fSec->XLast();ix++) {
      int val = (int) fImage->RdFrame(ix,iy);
      for (int bit=0;bit<16;bit++) {
        hADC->Fill(bit,(val>>bit)&1);
      }
    }
  
  hADC->Write(histName);
  delete hADC;
}
