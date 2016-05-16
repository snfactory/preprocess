/* === Doxygen Comment ======================================= */
/*! 
 * \file          overscan.cxx
 * \copyright     (c) 2004 SNIFS Collaboration
 * \date          Wed Aug  6 17:44:35 2003
 * \author        Emmanuel Gangler <e.gangler@ipnl.in2p3.fr>
 * \version       0.0
 * \brief         
 *                
 * $Id$
 **/
/* =========================================================== */

/* ----- IFU include ----- */
// IFU_io.h because SPECTRUM is needed by IFU_math.h
#include "IFU_io.h"
#include "IFU_math.h"
#include "gsl/gsl_fit.h"

/* ----- local includes -----*/
#include "image.hxx"
#include "imagesnifs.hxx"
#include "bichip.hxx"
#include "overscan.hxx"
#include "section.hxx"
#include "algocams.hxx"
#include "utils.h"

/* ##### OVERSCAN BASE ################################################# */

/* The big idea behind the overscan correction is that there is an image,
   with a section (BiasSec) which contains the data necessary to correct the 
   overscan, and another region (SubSec) which is the section where the 
   correction has to ba applied */

/* ===== constructor/destructor ======================================= */

/* ----- OverscanBase -------------------------------------------------- */
OverscanBase::OverscanBase() {

  fImage=0;
  fBiasSec=0;
  fSubSec=0;
  fLineZero=0;
  fLineLength=0;
  fNlines=5;
}

/* ----- ~OverscanBase ---------------------------------------- */
OverscanBase::~OverscanBase() {
  SetBiasSec(0);
  SetSubSec(0);
}


/* ===== setters ======================================= */

/* ----- SetBiasSec ---------------------------------------- */
void OverscanBase::SetBiasSec(Section* Sec) {
  if (fBiasSec)
    delete fBiasSec;
  fBiasSec=Sec;
}


/* ----- SetSubSec ---------------------------------------- */
void OverscanBase::SetSubSec(Section* Sec) {
  if (fSubSec)
    delete fSubSec;
  fSubSec=Sec;
}

/* ===== utilities methods ======================================= */

/* ----- Substract ------------------------------------------------ */

void OverscanBase::Substract(double *values, double * var) {
  return SubstractRamp(values,var);
}

/* ----- Substract Ramp ------------------------------------------------ */

void OverscanBase::SubstractRamp(double *ov, double * var) {
  // Substracts the overscan from fSubSec section of the image, 
  // taking the datas from the supplied window (that if values are those
  // computed from ComputeLines)

  // The assumption is that the overscan is recorder AFTER the data, 
  // and not before

  // 2 loops for in-place substraction

  // hack in order to minimize debugging
  Section* S = fBiasSec;
  // hack for inlining speed
  int SYFirst = S->YFirst();
  
  // Substract the overscan - makes a ramp between left and right
  // overscans. Left overscan = right line before. The fixed point
  // for the ramp are the middle ov overscan strip
  // A small extrapolation is made on one line after the fixed point
  // As it amounts to about 14 points, it is considered as a negligible
  // approximation 
  for (int col=fSubSec->YFirst();col<fSubSec->YLast();col++) {
    double ovGuessAdd, ovGuessMul;
    // result will be (Add + line * Mul )
    // no data before - be careful
    if (col<SYFirst) {
      ovGuessAdd=ov[0];
      ovGuessMul=0;
    }
    // no data after - no info
    else if (col>=S->YLast()) {
      ovGuessAdd = ov[S->YLength()-1];
      ovGuessMul = 0;
    }
    // first line : extrapolate the slope (best guess)
    else if (col==SYFirst) {
      // fLineZero = (S->X1()+S->X2())/2.0+1+Nx()
      // fLineLength = Nx()
      ovGuessMul = (ov[SYFirst+1]-ov[SYFirst]) / fLineLength ;
      ovGuessAdd = (ov[SYFirst+1]-ov[SYFirst])
        * ( - fLineZero +fLineLength ) / fLineLength 
        + 2*ov[SYFirst] - ov[SYFirst+1];;
    }
    // Ok, do it.
    else {
      ovGuessAdd = (ov[col-SYFirst]-ov[col-SYFirst-1])
        * ( - fLineZero +fLineLength) / fLineLength 
        + ov[col-SYFirst-1];
      ovGuessMul = (ov[col-SYFirst]-ov[col-SYFirst-1]) / fLineLength ;
    }
    for (int line=fSubSec->XFirst();line<fSubSec->XLast();line++) {
      double val = fImage->RdFrame(line,col) -ovGuessAdd - ovGuessMul * line;
      fImage->WrFrame(line,col,val);   
    } // for line
  }

  // And now the effect on the variance
  if (fImage->Variance()) {
    for (int col=fSubSec->YFirst();col<fSubSec->YLast();col++) {
      double addedVarAdd, addedVarMul;
      // no data before - be careful
      if (col<=SYFirst) {
        addedVarAdd=var[0];
        addedVarMul=0;
      }
      // no data after - no info
      else if (col>=S->YLast()) {
        addedVarAdd = var[S->YLength()-1];
        addedVarMul=0;
      }
      // Ok, do it. As the values may be correlated, this is perhaps just wrong
      // So we make an educated guess, in case of full correlation
      // in the case of ImproveMedian
      else {
        addedVarAdd = (var[col-SYFirst]-var[col-SYFirst-1])
          * ( - fLineZero +fLineLength) / fLineLength 
          + var[col-SYFirst-1];
        addedVarMul = (var[col-SYFirst]-var[col-SYFirst-1]) / fLineLength ;
      }
      for (int line=fSubSec->XFirst();line<fSubSec->XLast();line++) {
        double val = fImage->Variance()->RdFrame(line,col) + addedVarAdd + addedVarMul * line;
        fImage->Variance()->WrFrame(line,col,val);   
      } // for line
    }
  } //if variance
}

/* ----- ComputeLines -------------------------------------------------- */

void OverscanBase::ComputeLines(double* ov, double * var) {

  // Computes the value of the overscan for 1 line
  // Various methods exists. 
  // The overscan pattern can not be treated as a bias frame.
  // Polynomial fit was proven not to be satisfactory.
  // The sauron method was : median of the mean of the lines
  // Mean of the line is not fully adequate as we saw that the first 
  // lines of the 
  // overscan are dangerous to use for the blueCCD, As they amount for 1% of
  // last data column. But median of the line needs a removed odd-even effect
  // if present.
  // Median of n lines is perhaps adequate in case of sudden bumps. Variants 
  // of the Lee filter may perhaps also be applied. 

  // for the moment : median window of 5(fNlines parameter) lines, 
  // mean of the line.
  // the size of the window is restricted near the edges of the section in order  // to have only symmetric windows : it really helps for the median.

  // These comments corresponds to the base algorithms.
  // Theu may be overwritten by derived classes ...

  ComputeLinesMean(ov, var);

  ImproveLinesMedian(ov, var);
  
}

/* ----- ComputeLinesMean -------------------------------------------------- */

void OverscanBase::ComputeLinesMean(double* ov, double * var) {

  // Computes the value of the overscan for 1 line
  // Various methods exists. 
  // The overscan pattern can not be treated as a bias frame.
  // Polynomial fit was proven not to be satisfactory.
  // The sauron method was : median of the mean of the lines
  // Mean of the line is not fully adequate as we saw that the first 
  // lines of the 
  // overscan are dangerous to use for the blueCCD, As they amount for 1% of
  // last data column. But median of the line needs a removed odd-even effect
  // if present.
  // Median of n lines is perhaps adequate in case of sudden bumps. Variants 
  // of the Lee filter may perhaps also be applied. 

  // for the moment : median window of 5(fNlines parameter) lines, 
  // mean of the line.
  // the size of the window is restricted near the edges of the section in order  // to have only symmetric windows : it really helps for the median.

  // Small hacks to minimize debugging
  Section * S = fBiasSec;
  double * ovColRaw = ov;

  //
  // First the value line per line
  //

  double ovLine[S->XLength()];
  int col;

  for (col=S->YFirst();col<S->YLast();col++) {

    for (int line=S->XFirst();line<S->XLast();line++)
      ovLine[line-S->XFirst()] = fImage->RdFrame(line,col);

    // implies the overscan is pure from spurious...
    // but we do not get enough data in order to know
    // what happens in the beginning of the overscan strip.
    // so the best would be the truncated mean with known variance method
    // -> which cuts ?
    ovColRaw[col-S->YFirst()] = ut_mean(ovLine,S->XLength());
  }

  // The associated variance :
  if (fImage->Variance()) {
    double rms;
    // not really clean. Variance effect should be stored in an array
    fImage->RdDesc("RDNOISE",DOUBLE,1,&rms);
    double linevar = rms*rms / S->XLength();
    for (int col=S->YFirst();col<S->YLast();col++) {
      var[col-S->YFirst()] = linevar;
    }
  }

}

/* ----- ComputeLinesMean -------------------------------------------------- */

void OverscanBase::ComputeLinesMode(double* ov, double* var) {

  // This routine computes the mode of the current line.
  // This is the best automatic 0 approximation in case there are no
  // pure overscan region.

  Section * S = fBiasSec;

  double * ovColRaw = ov;
  double ovLine[S->XLength()];
  int col;

  for (col=fBiasSec->YFirst();col<fBiasSec->YLast();col++) {

    for (int line=fBiasSec->XFirst();line<fBiasSec->XLast();line++)
      ovLine[line-fBiasSec->XFirst()] = fImage->RdFrame(line,col);

    // implies the overscan is pure from spurious...
    // but we do not get enough data in order to know
    // what happens in the beginning of the overscan strip.
    // so the best would be the truncated mean with known variance method
    // -> which cuts ?
    ovColRaw[col-fBiasSec->YFirst()] = ut_mode(ovLine,fBiasSec->XLength());
  }

  // The associated variance : no idea -> put 0 !!! (alternatively, a gaussian fit
  // around the mode could be used... but time consuming ! )
  if (fImage->Variance()) {
    for (int col=S->YFirst();col<S->YLast();col++) {
      var[col-S->YFirst()] = 0;
    }
  }
}

/* ----- ImproveLinesMedian -------------------------------------------------- */

void OverscanBase::ImproveLinesMedian(double* ov, double * var) {

  // Note that the ov has to be computed first with a ComputeLines
  //
  // This method tries to improve the estimate using a window in which
  // the median of the values is computed.
  // This has the disavantage to correlate the lines errors - but this 
  // should be a small effect if the overscan is computed with enough values.

  // Small hack to minimize debugging
  Section * S = fBiasSec;

  // whith this choice, this makes a theoretical error mean,median of 
  // less than
  // 1/10 of the readout error for our CCDs, which is affordable !
  // except of course at the beginning, but there the systematics
  // effects are unavoidable due to the strip structure (essentially at 
  // the beginning) 
  double buffer[fNlines*2+1];
  int index[fNlines*2+1];
  double ovColRaw[S->YLength()];

  /* copy of ov */
  memcpy (ovColRaw,ov,sizeof(double)*S->YLength());
  
  for (int col=S->YFirst() ; col<S->YLast() ; col++) {
    int nLines = fNlines;
    // we restrict the window near the edges in order to get a symetric
    // window around the value. Bad overscan guess else occurs.
    if (col-S->YFirst() < fNlines) 
      nLines = col-S->YFirst();
    if (S->YLast()-1-col < nLines)
      nLines = S->YLast()-1-col;    

    // now the loop for the median
    for (int iCol=col-nLines; iCol<=col+nLines;iCol++) {
      buffer[iCol-col+nLines] = ovColRaw[iCol];
    }
    ov[col-S->YFirst()]=median(buffer,nLines*2+1,index);
  } // for col

 // And now the effect on the variance :
 // only accurate if all lines have the same variance ! (which is the case
 // for normal preprocessing)
 // or put the right computation here ...
  if (fImage->Variance()) {
    double linevar = var[0];
    
    for (int col=S->YFirst();col<S->YLast();col++) {
      for (int line=S->XFirst();line<S->XLast();line++) {
        // no data before or after : we don't know anything
        // choice is to get the first parameter
        // we could have put the variance to infinity

        // nlines of the 2 windows
        int nLines = fNlines;
        if (col-S->YFirst() < fNlines) 
          nLines = col-S->YFirst();
        if (S->YLast()-1-col < nLines)
          nLines = S->YLast()-1-col;
        // the math are computed as if it were the median of
        // a square distribution (which is false, and only a good
        // approximation with big N)
        // As the lines are fully correlated, we shall 
        // add also a small something due to the interpolation
        // but the error we do here is nevertheless negligible
        // with respect to readout error.
        // formula is lineval / 3 / (N + 2)
        var[col-S->YFirst()] = linevar / (2*nLines+3) * 3 ; // N = 2n + 1
      } // for line
    }
  } //if variance
}

/* ----- SubstractOddEven -------------------------------------------------- */

void OverscanBase::SubstractOddEven(double * param, double sigcut) {

  //   Odd-even substraction works as follows
  //   1) compute odd-even for each line in the overscan strip
  // 2) interpolate a linear fit
  // 3) substract it
  // NOTE : if there is no odd-even effect, the total error added shall be 
  // negligible
  //
  // param is a return value of the fit parameters (2 values), 
  // sigcut the level of oulier removal (in sigmas)
  
  // Hacks for debugging minimization
  Section * S = fBiasSec;

  // fill the array
  double oddEven[S->XLength()/2];
  double array[S->YLength()];
  double weight[S->YLength()];
  int i,j,remain;
  
  double sigma=0;
  if (sigcut) // rustic outlier detection.
    sigma = fImage->SectionRms(S,ut_fraction_sigcut(1.0/(S->XLength()-0.5)));

  for (j= S->YFirst() ;j< S->YLast();j++) {
    /* fill the intermediate odd-even array to get rid of cosmics */
    for (i= S->XFirst() ; i<S->XLast()-1 ; i+=2) /* 'odd' */
      oddEven[(i-S->XFirst())/2] = fImage->RdFrame(i,j)-fImage->RdFrame(i+1,j);
    
    remain = S->XLength()/2;

    double * oddEvenRet = oddEven;
    if (sigcut>0)
      ut_trunc_sigma_known(&oddEvenRet,&remain,sigma,sigcut);
    
    /* now compute the mean */
    array[j- S->YFirst() ]=0;
    for(i=0;i<remain;i++)
      array[j- S->YFirst() ] += oddEvenRet[i];
    array[j- S->YFirst() ] /= remain;  
    weight[j- S->YFirst() ] = remain;
  }

  // linear fit to the array (with wrapper to gsl_fit_wlinear)
  double X[ S->YLength() ];
  double cov00,cov01,cov11,chi2;
  for (i=0;i< S->YLength() ;i++)
    X[i]=i;
  gsl_fit_wlinear (X, 1, weight, 1, array, 1, S->YLength(), 
                   param, param+1, &cov00, &cov01, &cov11, 
                   &chi2);
  
  // perform the substraction (on sub frame)
  double sub,val;
  for (j=fSubSec->YFirst();j<fSubSec->YLast();j++) {
    sub = (param[0] + (j- S->YFirst() )*param[1])/2;
    for ( i= fSubSec->XFirst() +  abs( fSubSec->XFirst() - S->XFirst() )%2;
         i<fSubSec->XLast() ; i+=2) {
      val = fImage->RdFrame(i,j) -sub;
      fImage->WrFrame(i,j,val);
    }
    for ( i=  fSubSec->XFirst() +  abs( fSubSec->XFirst() - S->XFirst() +1)%2;
          i<S->XLast() ; i+=2) {
      val = fImage->RdFrame(i,j) +sub;
      fImage->WrFrame(i,j,val);
    }
  }

  // Effect on the variance is negligible.

}


/* ##### OVERSCAN SNIFS ################################################# */

/* Version of the overscan substraction for a SNIFS Detcom-formatted image 
Otcome images should be first detcom-formatted for this to be used */

/* ===== constructor/destructor ======================================= */

OverscanSnifs::OverscanSnifs(): OverscanBase() {
}

/* ===== Methods ======================================== */

/* ----- SetImage ----------------------------------------*/
void OverscanSnifs::SetImage(ImageSnifs* Image) 
  /* sets the image, and compute the relevant sections */
{
  OverscanBase::SetImage(Image->Image());
  // fImage = Image; // beware of the tricky hiding !
  SetBiasSec( Image->Algo()->SafeOverscanStrip(Image) );
  SetSubSec( new Section (1,Image->Nx(),1,Image->Ny()) );
  SetLineZero((fBiasSec->X1()+fBiasSec->X2())/2.0-1);
  SetLineLength(Image->Nx());
}

/* ----- Correct ----------------------------------------*/
void OverscanSnifs::Correct(BiChipSnifs* BiChip) {
  for (int i=0;i<BiChip->NChips();i++) {
    Correct(BiChip->Chip(i));
  }
}


/* ----- Correct ----------------------------------------*/
void OverscanSnifs::Correct(ImageSnifs* Image) {

  SetImage(Image);
  OddEvenCorrect();
  // after odd-even correction : removal of the odd-even improves the noise
  AddOverscanVariance(); 
  SubstractOffset();  

}

/* ----- OverscanNoise ---------------------------------------- */
double OverscanSnifs::OverscanNoise() {

  return fImage->SectionRms(fBiasSec,0);
}

/* ----- AddOverscanNoise ---------------------------------------- */
void OverscanSnifs::AddOverscanVariance() {

  int ovscNoise;

  // check
  if (Image()->ParanoMode() && fImage->Variance()) {
    if (fImage->Variance()->RdIfDesc("OVSCNOIS",INT, 1, &ovscNoise) > 0 
        && ovscNoise ) {
      print_error(" OverscanSnifs::AddOverscanVariance already done in %s",fImage->Name());
      return;
    }
  }
  

  // work
  double rms = OverscanNoise();
  fImage->WrDesc("RDNOISE",DOUBLE,1,&rms);
  if (fImage->Variance()) {
    
    fImage->Variance()->Add(rms*rms);
    // status for the check
    ovscNoise=1;
    fImage->Variance()->WrDesc("OVSCNOIS",INT,1,&ovscNoise);
  }
  
}

/* ----- OddEvenCorrect ---------------------------------------- */
void OverscanSnifs::OddEvenCorrect() {

  // odd-even substraction
  //
  double param[2];
  param[0]=param[1]=0;
  
  if (fOddEven) {
    if (Image()->ParanoMode()) {
      // check fcalsses !
    
      // check type of CCD (is operation allowed ?)

      // check overscan was not already substracted
    // returns 0 if not done !
      if ( fImage->RdIfDesc("OEPARAM",DOUBLE, 2, param) > 0 ){
	print_error("OverscanSnifs::OddEvenCorrect Odd-Even already substracted for %s\n",fImage->Name());
	print_error("Nothing Done\n");
	return;
      }
    }

    SubstractOddEven(param,0);
  }

  fImage->WrDesc("OEPARAM",DOUBLE, 2, param);

}



/* ----- SubstractOffset ----------------------------------------*/
void OverscanSnifs::SubstractOffset() {

  int overscanDone;

  if (Image()->ParanoMode()) {
    // check fcalsses !

    // check overscan was not already substracted
    if ( fImage->RdIfDesc("OVSCDONE",INT, 1, &overscanDone) >= 0 
         && overscanDone ){
      print_error("Overscan already substracted for %s\n",fImage->Name());
      print_error("Nothing Done\n");
      return;
    }
  }

  //
  // Get the real overscan work done (core of the computation)
  //
  double* ov= new double[fBiasSec->YLength()];
  double* var= new double[fBiasSec->YLength()];
  ComputeLines(ov,var);
  Substract(ov,var);

  // back to keyword business
  overscanDone=1;
  fImage->WrDesc("OVSCDONE",INT,1,&overscanDone);
  // set the biasframe if allowed
  char obstype[lg_name+1];
  fImage->RdDesc("OBSTYPE",CHAR,lg_name+1,obstype);
  if (! strcmp(obstype,"BIAS") || Image()->GetFClass() == RAW_BIAS ) {
    int biasFrame=1;
    fImage->WrDesc("BIASFRAM",INT,1,&biasFrame);
  }
  if (Image()->GetFClass() == RAW_BIAS)
    Image()->SetFClass(BIAS_FRAME);
  
  // for the posterity : write the OVSC median value
  double ovscMed=ut_median(ov,fBiasSec->YLength());
  fImage->WrDesc("OVSCMED",DOUBLE,1,&ovscMed);

  // last : get the max value of the substracted offset
  double max=0;
  for (int i=0;i<fBiasSec->YLength();i++) {
    if (ov[i]>max)
      max=ov[i];
  }
  if ((2*ov[0]-ov[1])>max)
    max=2*ov[0]-ov[1];
  fImage->WrDesc("OVSCMAX",DOUBLE,1,&max);
  delete[] ov;
  delete[] var;
  

}

/* ##### OVERSCAN NONE ################################################# */


/* ===== Methods ======================================== */

/* ----- SetImage ----------------------------------------*/
void OverscanNone::SetImage(ImageSnifs* Image) 
{ // force subsequent crashes if misused !
  OverscanBase::SetImage(0);
  SetBiasSec( 0 );
  SetSubSec( 0 );
}


/* ----- Substract ------------------------------------------------ */

void OverscanNone::Substract(double *values, double * var) {
  return ;
}

/* ----- ComputeLines -------------------------------------------------- */

void OverscanNone::ComputeLines(double* ov, double * var) {
  return;
}

/* ----- OverscanNoise -------------------------------------------------- */

double OverscanNone::OverscanNoise() {
  return 0;
}

/* ----- ComputeLines -------------------------------------------------- */

void OverscanNone::OddEvenCorrect() {
  return;
}

/* ##### OVERSCAN FROM DATA ############################################## */

/* Version of the overscan for a non-

Version of the overscan substraction for a SNIFS Detcom-formatted image 
Otcome images should be first detcom-formatted for this to be used */

/* ===== constructor/destructor ======================================= */


/* ===== Methods ======================================== */

/* ----- SetImage ----------------------------------------*/
void OverscanFromData::SetImage(ImageSnifs* Image) 
  /* sets the image, and compute the relevant sections */
{
  OverscanBase::SetImage(Image->Image());
  // take only the datasec
  char dataSecString[lg_name+1];
  Image->RdDesc("DATASEC",CHAR,lg_name+1,dataSecString);
  SetSubSec( new Section (dataSecString ) );
  SetBiasSec(   new Section (dataSecString )); // same region
  SetLineZero((fBiasSec->X1()+fBiasSec->X2())/2.0-1);
  SetLineLength(Image->Nx());
}


/* ----- ComputeLines -------------------------------------------------- */

void OverscanFromData::ComputeLines(double* ov, double* var) {

  ComputeLinesMode(ov, var);

  // perhaps not the best choice, but this is a rescue class anyway
  ImproveLinesMedian(ov, var);

}

/* ----- OverscanNoise ---------------------------------------- */
double OverscanFromData::OverscanNoise() {

  // no computation done on the noise, as it shall be wrong.
  return 0;
}

/* ----- OverscanNoise ---------------------------------------- */
void OverscanFromData::OddEvenCorrect() {

  // no odd-even correction
  double param[2];
  param[0]=param[1]=0;
  fImage->WrDesc("OEPARAM",DOUBLE, 2, param);

  return ;
}
