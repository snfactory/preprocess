/* === Doxygen Comment ======================================= */
/*! 
 * \file          filter.cxx
 * \copyright     (c) 2004 SNIFS-Supernova Factory Experiment
 * \date          Wed Aug  6 17:44:35 2003
 * \author        Emmanuel Gangler <e.gangler@ipnl.in2p3.fr>
 * \version       0.0
 * \brief         
 *                
 * $Id$
 **/
/* =========================================================== */

/* ----- local includes ----- */
#include "filter.hxx"
#include "image.hxx"
#include "section.hxx"
#include "utils.h"
#include "analyser.hxx"

/* ##### Filter ################################################# */

void ImageFilter::Filter() {
  int i,j;
  Section S;
  
  switch (fBound) {

  case kNoData :
    for (j=0;j<fYsize;j++)
      for (i=0;i<fInput->Ny();i++)
        WriteNoData(i,j);
    for (j=fYsize;j<fInput->Ny()-fYsize;j++) {
      for (i=0;i<fXsize;i++ ) 
        WriteNoData(i,j);
      for (i=fXsize;i<fInput->Nx()-fXsize;i++) {
        S.SetXFirst(i-fXsize);
        S.SetXLast(i+fXsize);
        S.SetYFirst(j-fYsize);
        S.SetYLast(j+fYsize);
        Filter(i,j,&S);
      }
      for (i=fInput->Nx()-fXsize;i<fInput->Nx();i++)
        WriteNoData(i,j);
    }
    for (j=fInput->Ny()-fYsize;j<fInput->Ny();j++)
      for (i=0;i<fInput->Ny();i++)
        WriteNoData(i,j);
    break;

  case kShrinks :
    for (j=0;j<fInput->Ny();j++) {
      S.SetYFirst(j - fYsize);
      S.SetYLast(j + fYsize+1);
      if (S.YFirst()<0)
        S.SetYFirst(0);
      if (S.YLast()>fInput->Ny())
        S.SetYLast(fInput->Ny());
      for (i=0;i<fInput->Ny();i++) {
        S.SetXFirst(i - fXsize);
        S.SetXLast(i + fXsize+1);
        if (S.XFirst()<0)
          S.SetXFirst(0);
        if (S.YLast()>fInput->Ny())
          S.SetYLast(fInput->Ny());
        Filter(i,j,&S);
      }
    }
    
  case kSymetric :
    for (j=0;j<fInput->Ny();j++) {
      S.SetYFirst(j - fYsize);
      S.SetYLast(j + fYsize+1);
      if (S.YFirst()<0) {
        S.SetYFirst(0);
        S.SetYLast(2*j+1);
      }
      if (S.YLast()>fInput->Ny()) {
        S.SetYLast(fInput->Ny());
        S.SetYFirst(- fInput->Ny() + 2*j + 1);
      }
      for (i=0;i<fInput->Ny();i++) {
        S.SetXFirst(i - fXsize);
        S.SetXLast(i + fXsize+1);
        if (S.XFirst()<0) {
          S.SetXFirst(0);
          S.SetYLast(2*i+1);
        }
        if (S.YLast()>fInput->Ny()) {
          S.SetYLast(fInput->Ny());
          S.SetYFirst(-fInput->Nx() + 2*i +1);
        }
        Filter(i,j,&S);
      }
    }
    break;
  } 
}


/* ----- WriteNoData ------------------------------ */
void ImageFilter::WriteNoData(int I, int J) 
{
  fOutput->WrFrame(I,J,fNoData);
  if (fOutput->Variance())
    fOutput->Variance()->WrFrame(I,J,fNoDataVar);
}




/* ##### ImageFilterHF ################################################# */


/* ----- constructor -------------------------------------------------- */
ImageFilterHF::ImageFilterHF(int Xsize, int Ysize, Bound_t B) {
  fXsize=Xsize;
  fYsize=Ysize;
  fBound = B;
  fThreshold = -ut_big_value;
  fSignificance=0;
  fAnal=new ImageAnalyser();
}

/* ----- destructor ------------------------------ */
ImageFilterHF::~ImageFilterHF() {
  delete fAnal;
}

/* ----- SetInputImage ---------------------------------------- */
void ImageFilterHF::SetInputImage(ImageSimple* I){
  fInput=I;
  fAnal->SetImage(I);
}


/* ----- Filter -------------------------------------------------- */
void ImageFilterHF::Filter(int i, int j, Section* S) {
  if (fInput->RdFrame(i,j)<fThreshold) {
    WriteNoData(i,j);
    return;
  }  
  fAnal->SetSection(S);
  double median = fAnal->Quantile(0.5);
  // this is the variance of the mean, with actual point excluded, not the variance of the median.
  double pointvar = fInput->Variance()->RdFrame(i,j);
  double var = 1/(1/fAnal->MeanMapVariance() - 1/pointvar);
  int secSize = S->XLength()*S->YLength();
  // now rough variance of the median :
  var *= 3 * secSize / (secSize + 2 );
  
  // we consider 2 cases : 
  // the central point = the median. In that caes, we have no correction
  // to perform (suspicion variation of signal over the window.)
  double signif=fabs( fInput->RdFrame(i,j) - median) / sqrt(var + pointvar);
  if (signif < fSignificance)
    WriteNoData(i,j);
  else {
    if (median != 0) {
      double val=fInput->RdFrame(i,j)/median;
      fOutput->WrFrame(i,j,val);
      fOutput->Variance()->WrFrame(i,j,(var* val*val + pointvar)/ (median*median));
    }
    else {
      WriteNoData(i,j);
    }
  }
}

