/* === Doxygen Comment ======================================= */
/*! 
 * \file          analyser.cxx
 * \copyright     (c) 2003 SNIFS-Supernova Factory Experiment
 * \date          Wed Aug  6 17:44:35 2003
 * \author        Emmanuel Gangler <e.gangler@ipnl.in2p3.fr>
 * \version       0.0
 * \brief         
 *                
 * $Id$
 **/
/* =========================================================== */

#include "gsl/gsl_statistics.h"

#include "section.hxx"
#include "image.hxx"
#include "imagesnifs.hxx"
#include "analyser.hxx"
#include "utils.h"

/* ##### ImageAnalyser ################################################# */

/* ===== constructor/destructor ======================================= */
/* ----- void constructor  ----------------------------------*/
ImageAnalyser::ImageAnalyser(){
  fVal=0;
  fSec=0;
  fImage=0;
  fFftLength=0;
}

/* ----- constructor Image, Sec ----------------------------------*/
ImageAnalyser::ImageAnalyser(ImageSimple * Image, Section *Sec){
  fVal=0;
  fSec=0;
  fFftLength=0;
  SetImage(Image);
  SetSection(Sec);
}

/* ----- constructor Image, Sec ----------------------------------*/
ImageAnalyser::ImageAnalyser(ImageSnifs * Image, Section *Sec) {
  fVal=0;
  fSec=0;
  fFftLength=0;
  SetImage(Image);
  SetSection(Sec);
  // ImageAnalyser(Image,Sec); -> does not work ...
}

/* ----- destructor --------------------------------------------*/
ImageAnalyser::~ImageAnalyser(){
  if (fVal) {
    delete[] fVal;
  }
  if (fFftLength) {
     gsl_fft_real_wavetable_free (fReal);
     gsl_fft_real_workspace_free (fWork);
     fFftLength=0;
  }
}


/* ===== setters ====================================================== */
/* ----- SetSection ----------------------------------------------------*/
void ImageAnalyser::SetSection(Section* Sec){
  // Set only the section if it belongs to the image.
  // (now handled by FillVal)
  fSec=Sec;
  ResetVal();
}

/* ----- SetImage ----------------------------------------------------*/
void ImageAnalyser::SetImage(ImageSimple * Image){
  // Set only the image if it fits the section.

  fImage=Image;
  ResetVal();
}

/* ----- SetImage ----------------------------------------------------*/
void ImageAnalyser::ResetVal(){
  // Set only the image if it fits the section.
  if (fVal)
    delete[] fVal;
  fVal=0;
  fNVal=0;
  fMeanDone=0;
  fVarDone=0;
  fIsSorted=0;
}


/* ----- SetImage ----------------------------------------------------*/
void ImageAnalyser::FillVal(){
  // Set only the image if it fits the section.
  ResetVal();

  if (fSec && fSec->XFirst()>=0 && fSec->YFirst()>=0 && fSec->XLast()<=fImage->Nx() && fSec->YLast()<=fImage->Ny()) {
    fNVal = fSec->XLength()*fSec->YLength();
    fVal = new double[fNVal];
    for (int iy=fSec->YFirst();iy<fSec->YLast();iy++)
      for (int ix=fSec->XFirst();ix<fSec->XLast();ix++) {
        fVal[ix-fSec->XFirst()+(iy-fSec->YFirst())*fSec->XLength()] 
            = fImage->RdFrame(ix,iy);
      }
  } else {
    print_warning("Analyser::FillVal Section does not match image");
  }
  
}

/* ===== Analysis ======================================= */
/* ----- MeanLevel -------------------------------------- */
double ImageAnalyser::MeanLevel(){
  if (!Val())
    FillVal();
  if (fMeanDone)
    return fMean;
  fMeanDone=1;
  return fMean = ut_mean(fVal,fNVal);
}

/* ----- MeanLevel -------------------------------------- */
double ImageAnalyser::StatsVariance(){
  if (fVarDone)
    return fVar;
  
  fVarDone=1;
  double mean = MeanLevel();    
  return fVar = gsl_stats_variance_m(fVal,1,fNVal,mean);
}

/* ----- MeanMapVariance ------------------------------ */
double ImageAnalyser::MeanMapVariance() {
  // Returns the variance computed from the variance map.

  // force crash ?
  //if (!fImage->Variance()) {
  //  return -1;
  //}
  double weight=0;
  for (int j=fSec->YFirst(); j< fSec->YLast();j++) 
    for (int i=fSec->XFirst(); i< fSec->XLast();i++) {
      weight += 1/fImage->Variance()->RdFrame(i,j);
    }
  return 1/weight;
}
  

/* ----- NPixOut -------------------------------------- */
int ImageAnalyser::NPixOut(double SigmaCut){
  double mean = MeanLevel();
  double sigma = sqrt(StatsVariance());
  int npix=0;
  for (int i=0;i<fNVal;i++) {
    if (fabs(fVal[i]-mean)/sigma > SigmaCut)
      npix++;
  }
  return npix;
}

/* ----- NPixOver -------------------------------------- */
int ImageAnalyser::NPixOver(double Limit){
  int npix=0;
  if (!Val())
    FillVal();
  for (int i=0;i<fNVal;i++) {
    if (fVal[i] > Limit)
      npix++;
  }
  return npix;
}

/* ----- OutPixMean -------------------------------------- */
double ImageAnalyser::OutPixMean(double SigmaCut){
  double mean = MeanLevel();
  double sigma = sqrt(StatsVariance());
  int npix=0;
  double meanOut=0;
  for (int i=0;i<fNVal;i++) {
    if (fabs(fVal[i]-mean)/sigma > SigmaCut) {
      npix++;
      meanOut+=fVal[i];
    }
  }
  if (npix)
    return meanOut/npix;
  else return 0;
}

/* ----- Quantile -------------------------------------- */
double ImageAnalyser::Quantile(double Q){
  // returns the value which realizes a fraction Q of the data under it
  // Q = 0.5 corresponds to the median

  // Note that the sorting is rather expensive !

  if (!Val())
    FillVal();
  if (!fIsSorted) {
    ut_sort_ascending(fVal,fNVal);
    fIsSorted=1;
  }
  
  return gsl_stats_quantile_from_sorted_data(fVal,1,fNVal,Q);
 
}

/* ----- SigmaClippedInfo -------------------------------------- */
void ImageAnalyser::SigmaClippedInfo(double SigmaCut,double * Mean, double * Rms, int * Nout){

  if (!Val())
    FillVal();
  if (!fIsSorted) {
    ut_sort_ascending(fVal,fNVal);
    fIsSorted=1;
  }

  double * valRef = fVal;
  int nValRef = fNVal;
  ut_trunc_sigma_unknown_fast_sorted(&valRef,&nValRef,SigmaCut);
  
  *Mean = ut_mean(valRef,nValRef);
  *Rms = gsl_stats_sd_m(valRef,1,nValRef,*Mean);
  *Nout = fNVal - nValRef;
}

/* ----- LineFft -------------------------------------- */
ImageSnifs * ImageAnalyser::LineFft(char* outName){
  
  // prepare the gsl fft
  if (fFftLength && fFftLength!=fSec->XLength()) {
     gsl_fft_real_wavetable_free (fReal);
     gsl_fft_real_workspace_free (fWork);
     fFftLength=0;
  }
  fFftLength=fSec->XLength();
  fWork = gsl_fft_real_workspace_alloc (fSec->XLength());
  fReal = gsl_fft_real_wavetable_alloc (fSec->XLength());
  double ffts[fSec->XLength()];

  // builds the output image
  ImageSnifs* out = new ImageSnifs();
  out->CreateFrame(outName,fSec->XLength(),fSec->YLength());
  out->Image()->ImportHeader(fImage);
  out->DeleteDesc("DATASEC");
  out->DeleteDesc("BIASSEC");
  for (int j=fSec->YFirst();j<fSec->YLast();j++) {

    for (int i=fSec->XFirst();i<fSec->XLast();i++)
      ffts[i-fSec->XFirst()] =fImage->RdFrame(i,j);
    
    gsl_fft_real_transform (ffts, 1, fSec->XLength(), fReal, fWork);
    for (int i=fSec->XFirst();i<fSec->XLast();i++) 
      out->WrFrame(i-fSec->XFirst(),j-fSec->YFirst(),ffts[i-fSec->XFirst()]);  
  }
  return out;
}
