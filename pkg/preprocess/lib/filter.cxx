
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
//#include "IFU_io.h"

/* ##### Filter ################################################# */

void ImageFilter::Filter() {
  int i,j;
  Section S;

  reset_print_progress();
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
        print_progress("Filtering",100*(i+j*fInput->Nx())/(fInput->Ny()*fInput->Nx()*1.0),1.0);      
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
      for (i=0;i<fInput->Nx();i++) {
        S.SetXFirst(i - fXsize);
        S.SetXLast(i + fXsize+1);
        if (S.XFirst()<0)
          S.SetXFirst(0);
        if (S.XLast()>fInput->Nx())
          S.SetXLast(fInput->Nx());
        print_progress("Filtering",100*(i+j*fInput->Nx())/(fInput->Ny()*fInput->Nx()*1.0),1.0);      
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
      for (i=0;i<fInput->Nx();i++) {
        S.SetXFirst(i - fXsize);
        S.SetXLast(i + fXsize+1);
        if (S.XFirst()<0) {
          S.SetXFirst(0);
          S.SetXLast(2*i+1);
        }
        if (S.XLast()>fInput->Nx()) {
          S.SetXLast(fInput->Nx());
          S.SetXFirst(-fInput->Nx() + 2*i +1);
        }
        print_progress("Filtering",100*(i+j*fInput->Nx())/(fInput->Ny()*fInput->Nx()*1.0),1.0);      
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
  // variance of the mean
  double var = fAnal->MeanMapVariance();
  int secSize = S->XLength()*S->YLength();
  // now rough (guess) variance of the median :
  var *= 3.0 * secSize / (secSize + 2.0 );
  
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


/* ##### ImageFilterMedian ################################################# */


/* ----- constructor -------------------------------------------------- */
ImageFilterMedian::ImageFilterMedian(int Xsize, int Ysize, Bound_t B) {
  fXsize=Xsize;
  fYsize=Ysize;
  fBound = B;
  fAnal=new ImageAnalyser();
}

/* ----- destructor ------------------------------ */
ImageFilterMedian::~ImageFilterMedian() {
  delete fAnal;
}

/* ----- SetInputImage ---------------------------------------- */
void ImageFilterMedian::SetInputImage(ImageSimple* I){
  fInput=I;
  fAnal->SetImage(I);
}


/* ----- Filter -------------------------------------------------- */
void ImageFilterMedian::Filter(int i, int j, Section* S) {

  fAnal->SetSection(S);
  double median = fAnal->Quantile(0.5);
  // this is the variance of the mean, with actual point excluded, not the variance of the median.
  fOutput->WrFrame(i,j,median);

  if (fInput->Variance() && fOutput->Variance()) {
    
    // variance of the mean
    double var = fAnal->MeanMapVariance();
    int secSize = S->XLength()*S->YLength();
    // now rough (guess) variance of the median :
    var *= 3.0 * secSize / (secSize + 2.0 );
  
    fOutput->Variance()->WrFrame(i,j,var);
  }
}

/* ##### ImageFilterSigmaClip ################################################# */


/* ----- constructor -------------------------------------------------- */
ImageFilterSigmaClip::ImageFilterSigmaClip(int Xsize, int Ysize, Bound_t B) {
  fXsize=Xsize;
  fYsize=Ysize;
  fBound = B;
  fAnal=new ImageAnalyser();
}

/* ----- destructor ------------------------------ */
ImageFilterSigmaClip::~ImageFilterSigmaClip() {
  delete fAnal;
}

/* ----- SetInputImage ---------------------------------------- */
void ImageFilterSigmaClip::SetInputImage(ImageSimple* I){
  fInput=I;
  fAnal->SetImage(I);
}


/* ----- Filter -------------------------------------------------- */
void ImageFilterSigmaClip::Filter(int i, int j, Section* S) {


  double mean,rms;
  int nout;
  fAnal->SetSection(S);

  if (fInput->Variance()) 
    fAnal->SigmaClippedInfoVarKnown(fSigma,&mean,&rms);
  else
    fAnal->SigmaClippedInfo(fSigma,&mean,&rms,&nout);

  fOutput->WrFrame(i,j,mean);
  if ( fOutput->Variance()) {
    fOutput->Variance()->WrFrame(i,j,rms*rms);
  }
}

/* ##### ImageFilterMax ################################################# */


/* ----- constructor -------------------------------------------------- */
ImageFilterMax::ImageFilterMax(int Xsize, int Ysize, Bound_t B) {
  fXsize=Xsize;
  fYsize=Ysize;
  fBound = B;
}

/* ----- destructor ------------------------------ */
ImageFilterMax::~ImageFilterMax() {
}

/* ----- SetInputImage ---------------------------------------- */
//void ImageFilterMax::SetInputImage(ImageSimple* I){
//  fInput=I;
//}


/* ----- Filter -------------------------------------------------- */
void ImageFilterMax::Filter(int i, int j, Section* S) {

  double max = -ut_big_value;
  int ix0,iy0;
  for (int iy=S->YFirst();iy<S->YLast();iy++) {
    for (int ix=S->XFirst();ix<S->XLast();ix++) {
      //      if (ix<0 || ix>=fInput->Nx() || iy<0 || iy > fInput->Ny())
      //  print_error("ImageFilterMax::Filter i,j = %d,%d",ix,iy );
      
      if (fInput->RdFrame(ix,iy)>max) {
        max = fInput->RdFrame(ix,iy);
        ix0=ix;
        iy0=iy;
      }
    }
  }

  fOutput->WrFrame(i,j,max);
  // very crude estimate now ...
  if (fInput->Variance()&&fOutput->Variance()) {
    fOutput->Variance()->WrFrame(i,j,
           fInput->Variance()->RdFrame(ix0,iy0));
  }
}

/* ##### ImageFilterLaplacian ################################################# */


/* ----- constructor -------------------------------------------------- */
ImageFilterLaplacian::ImageFilterLaplacian(int Xsize, int Ysize, Bound_t B) {
  fXsize=1;
  fYsize=1;
  fBound = kNoData;
}

/* ----- destructor ------------------------------ */
ImageFilterLaplacian::~ImageFilterLaplacian() {
}

/* ----- SetInputImage ---------------------------------------- */
//void ImageFilterMax::SetInputImage(ImageSimple* I){
//  fInput=I;
//


/* ----- Filter -------------------------------------------------- */
void ImageFilterLaplacian::Filter() {
  double noise[2];
  if (!fInput->Variance())
    print_error("ImageFilterLaplacian::Filter works only on a slow preprocessed file");
  fInput->RdDesc("RDNOISE",DOUBLE,2,noise);
  if (noise[1]>noise[0])
    fCut=noise[1];
  else
    fCut=noise[0];
  

  for (int j=1;j<fInput->Ny()-1;j++){
    for (int i=1;i<fInput->Nx()-1;i++)
      Filter(i,j,0);
  }
  for (int j=0;j<fInput->Ny();j++){
    WriteNoData(0,j);
    WriteNoData(fInput->Nx()-1,j);
  }
  for (int i=0;i<fInput->Nx();i++) {
    WriteNoData(i,0);
    WriteNoData(i,fInput->Ny()-1);
  }
}


/* ----- Filter -------------------------------------------------- */
void ImageFilterLaplacian::Filter(int i0, int j0, Section* S) {

  double nsigma=3.0;

  if (fInput->RdFrame(i0,j0)<nsigma*fCut){
    fOutput->WrFrame(i0,j0,0);
    return;
  }
  
  // i-1,i+1 and j-1, j+1 must have been checked before
  double laplacian=fInput->RdFrame(i0,j0) -0.25*
    ( fInput->RdFrame(i0-1,j0) + fInput->RdFrame(i0+1,j0)
    + fInput->RdFrame(i0,j0-1) + fInput->RdFrame(i0,j0+1));
  ImageSimple* var=fInput->Variance();
  double varlap=var->RdFrame(i0,j0) + 0.0625*
    (var->RdFrame(i0-1,j0) + var->RdFrame(i0+1,j0)
     + var->RdFrame(i0,j0-1) + var->RdFrame(i0,j0+1));
  if (varlap<ut_big_value*0.001 && fabs(laplacian)<3*sqrt(varlap)){
    fOutput->WrFrame(i0,j0,0);
    return;
  }
  
  double buffer[9];
  for (int j=0;j<3;j++){
    for (int i=0;i<3;i++) {
      buffer[i+j*3]=fInput->RdFrame(i0-1+i,j0-1+j);
    }
    
  }
  double median = ut_median(buffer,9);
  if (median<fCut)
    median=fCut;
  fOutput->WrFrame(i0,j0,laplacian/median);
}

/* ##### ImageFilterLaplacian ################################################# */


/* ----- constructor -------------------------------------------------- */
ImageFilterRemCosmic::ImageFilterRemCosmic(int Xsize, int Ysize, Bound_t B):ImageFilterLaplacian(Xsize,Ysize,B) {
  double fRatio=2.0;//2.0 is safe for arcs ... insufficient for R dark

}

/* ----- destructor ------------------------------ */
ImageFilterRemCosmic::~ImageFilterRemCosmic() {
}

/* ----- SetInputImage ---------------------------------------- */
//void ImageFilterMax::SetInputImage(ImageSimple* I){
//  fInput=I;
//

/* ----- Filter -------------------------------------------------- */
void ImageFilterRemCosmic::Filter() {
  double noise[2];
  if (!fInput->Variance())
    print_error("ImageFilterRemCosmic::Filter needs variance");
  fInput->RdDesc("RDNOISE",DOUBLE,2,noise);
  if (noise[1]>noise[0])
    fCut=noise[1];
  else
    fCut=noise[0];
  

  for (int j=1;j<fInput->Ny()-1;j++){
    for (int i=1;i<fInput->Nx()-1;i++) {
      fRemoved=0;
      Filter(i,j,0);
      if (fRemoved) {
        if (i>1)
          i-=2;
        if (j>1)
          j--;
      }
    }
  }
}

/* ----- Filter -------------------------------------------------- */
void ImageFilterRemCosmic::Filter(int i0, int j0, Section* S) {

  double nsigma=3.0;

  if (fOutput->RdFrame(i0,j0)<nsigma*fCut){
    return;
  }
  
  // i-1,i+1 and j-1, j+1 must have been checked before
  double laplacian=fOutput->RdFrame(i0,j0) -0.25*
    ( fOutput->RdFrame(i0-1,j0) + fOutput->RdFrame(i0+1,j0)
    + fOutput->RdFrame(i0,j0-1) + fOutput->RdFrame(i0,j0+1));
  ImageSimple* var=fOutput->Variance();
  double varlap=var->RdFrame(i0,j0) + 0.0625*
    (var->RdFrame(i0-1,j0) + var->RdFrame(i0+1,j0)
     + var->RdFrame(i0,j0-1) + var->RdFrame(i0,j0+1));
  if (varlap<ut_big_value*0.001 && fabs(laplacian)<3*sqrt(varlap)){
    return;
  }
  
  int maxi=3;
  int maxj=3;
  double buffer[maxi*maxj];
  for (int j=0;j<maxj;j++){
    for (int i=0;i<maxi;i++) {
      buffer[i+j*maxi]=fOutput->RdFrame(i0-maxi/2+i,j0-maxj/2+j);
    }
    
  }
  double median = ut_median(buffer,maxi*maxj);
  if (median<fCut)
    median=fCut;
  if (laplacian/median>fRatio) {
    int maxj=9;
    double buffer[maxi*maxj];
    for (int j=0;j<maxj;j++){
      for (int i=0;i<maxi;i++) {
        buffer[i+j*maxi]=fOutput->RdFrame(i0-maxi/2+i,j0-maxj/2+j);
      }
      median = ut_median(buffer,maxi*maxj);
    }
    fOutput->WrFrame(i0,j0,median);
    fOutput->Variance()->WrFrame(i0,j0,ut_big_value);
    fRemoved=1;
  }
  
}


