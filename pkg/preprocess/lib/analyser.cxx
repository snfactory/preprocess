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
  if (fVal)
    delete[] fVal;

  if (fImage && Sec->XFirst()>=0 && Sec->YFirst()>=0 && Sec->XLast()<=fImage->Nx() && Sec->YLast()<=fImage->Ny()) {
    fSec = Sec;
    fNVal = fSec->XLength()*fSec->YLength();
    fVal = new double[fNVal];
    for (int iy=fSec->YFirst();iy<fSec->YLast();iy++)
      for (int ix=fSec->XFirst();ix<fSec->XLast();ix++) {
        fVal[ix-fSec->XFirst()+(iy-fSec->YFirst())*fSec->XLength()] 
            = fImage->RdFrame(ix,iy);
      }
  }
  else {
    fSec=Sec;
    fVal=0;
    fNVal=0;
  }
}

/* ----- SetImage ----------------------------------------------------*/
void ImageAnalyser::SetImage(ImageSimple * Image){
  // Set only the image if it fits the section.
  if (fVal)
    delete[] fVal;

  if (fSec && fSec->XFirst()>=0 && fSec->YFirst()>=0 && fSec->XLast()<=Image->Nx() && fSec->YLast()<=Image->Ny()) {
    fImage = Image;
    fNVal = fSec->XLength()*fSec->YLength();
    fVal = new double[fNVal];
    for (int iy=fSec->YFirst();iy<fSec->YLast();iy++)
      for (int ix=fSec->XFirst();ix<fSec->XLast();ix++) {
        fVal[ix-fSec->XFirst()+(iy-fSec->YFirst())*fSec->XLength()] 
            = fImage->RdFrame(ix,iy);
      }
  }
  else {
    fImage=Image;
    fVal=0;
    fNVal=0;
  }
}

/* ===== Analysis ======================================= */
/* ----- MeanLevel -------------------------------------- */
double ImageAnalyser::MeanLevel(){
  return ut_mean(fVal,fNVal);
}

/* ----- MeanLevel -------------------------------------- */
double ImageAnalyser::StatsVariance(){
  return gsl_stats_variance(fVal,1,fNVal);
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

/* ----- OutPixMean -------------------------------------- */
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
