/* === Doxygen Comment ======================================= */
/*! 
 * \file          analyser.hxx
 * \copyright     (c)  2003 SNIFS-Supernova Factory Experiment
 * \date          Wed Aug  6 18:32:01 2003
 * \author        Emmanuel Gangler <e.gangler@ipnl.in2p3.fr>
 * \version       0.0
 * \brief         
 *                
 * $Id$
 **/
/* =========================================================== */

#ifndef ANALYSER_H
#define ANALYSER_H
/*
The analyser contains various tools to analyse an image.
Both the image and the section are not actual members,
but just references.
Note that the Image shall not be modified while analyzed.
*/

class ImageSimple;
class Section;
#include "gsl/gsl_fft_real.h"
#include "gsl/gsl_fft_halfcomplex.h"
#include "imagesnifs.hxx"
/* ===== IMAGE ANALYSER ============================== */

class ImageAnalyser {
public :
  
  ImageAnalyser();
  ImageAnalyser(ImageSimple * Image, Section *Sec);
  ImageAnalyser(ImageSnifs * Image, Section *Sec);
  ~ImageAnalyser();
  
  
  /* ----- Setters ---------------------------------------- */
  void SetImage(ImageSimple * Image);
  void SetImage(ImageSnifs * Image) {SetImage(Image->Image());}
  void ResetVal();
  void SetSection(Section * Section); 
  void FillVal();

  /* ----- Getters ---------------------------------------- */
  double *Val() {return fVal;}

  /* ----- Analysis ---------------------------------------- */
  double MeanLevel();
  double StatsVariance();
  double MeanMapVariance();
  int NPixOut(double SigmaCut);
  int NPixOver(double Value);
  double OutPixMean(double SigmaCut);
  double Quantile(double Q);
  void SigmaClippedInfo(double SigmaCut,double * Mean, double * Rms, int * Nout);
  ImageSnifs* LineFft(char* outName);

  
protected :
  ImageSimple * fImage;
  Section * fSec;
  double * fVal;
  int fNVal;

  int fMeanDone;
  double fMean;
  int fVarDone;
  double fVar;
  int fIsSorted;
   
  gsl_fft_real_wavetable * fReal;
  gsl_fft_real_workspace * fWork;
  int fFftLength;
};

#endif
