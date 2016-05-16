/* === Doxygen Comment ======================================= */
/*! 
 * \file          overscan.hxx
 * \copyright     (c) 2004 SNIFS Collaboration
 * \date          Wed Aug  6 18:32:01 2003
 * \author        Emmanuel Gangler <e.gangler@ipnl.in2p3.fr>
 * \version       0.0
 * \brief         
 *                
 * $Id$
 **/
/* =========================================================== */

#ifndef OVERSCAN_H
#define OVERSCAN_H

/* Overscan is the interface to various possible overscan substraction 
   algotithms. 
   The assumptions are that there is an overscan region (defined as a section)
   of the image that may be used to substrcat the 0 in some region of the image
   (defined as a section, not necessarily the datasec - so the SubSec name).

   The main methods are the determination of the overscan level
   And the substraction of it.
   The overscan substraction is an operation of the image to itself

*/

/* ----- forward definitions ----- */
class ImageSimple;
class BiChipSnifs;

/* ----- local includes ----- */
#include "imagesnifs.hxx"


/* ===== OverscanBase ========================================*/

class OverscanBase {

public:

  OverscanBase();
  virtual ~OverscanBase();

  /* setters */ 
  virtual void SetImage(ImageSimple* image) {fImage = image;};
  virtual void SetBiasSec(Section* Sec);
  virtual void SetSubSec(Section* Sec);
  virtual void SetLineZero(double Zero) {fLineZero = Zero;};
  virtual void SetLineLength(double Length) {fLineLength = Length;};

  /* no interface method ... */
  // The interface has to be specified with the true image class, 
  // i.e. ImageSnifs
  // virtual void Correct()=0; // should remain like this
  
  /*utilities methods*/

  virtual void ComputeLinesMean(double * values, double * var);
  virtual void ComputeLinesMode(double * values, double * var);
  virtual void ImproveLinesMedian(double* values, double * var);
  virtual void ComputeLines(double * values, double * var);
  virtual void Substract(double* values, double * var);
  virtual void SubstractRamp(double* values, double * var);
  virtual void SubstractOddEven(double * param, double sigcut);

protected:
  ImageSimple * fImage;
  Section* fSubSec;
  Section* fBiasSec;
  // the physical position on the line (in subsec frame) 
  // where the overscan is effectively computed (used by SubsctractRamp)
  double fLineZero; 
  // The physical number of pixels read through the amplifier between 2 lines
  // (used by SubstractRamp)
  double fLineLength;

  // some parameters
  int fNlines; // number of lines used to reduce the noise
};


/* ===== OverscanSnifs ========================================*/

/* The basic overscan substractor for SNIFS images = the one interfaced
by Preprocessor*/

class OverscanSnifs : public OverscanBase {

public:  
  OverscanSnifs();
  virtual ~OverscanSnifs(){}

  // setters and useful getters
  void SetOddEven(int IsTrue) {fOddEven = IsTrue;}
  ImageSnifs* Image() {return dynamic_cast<ImageSnifs *>(fImage);}
  
  
  virtual void SetImage(ImageSnifs* Image) ;
  virtual void Correct(ImageSnifs* Image);
  virtual void Correct(BiChipSnifs* BiChip);
  virtual double OverscanNoise();
  virtual void AddOverscanVariance();
  virtual void OddEvenCorrect();
  virtual void SubstractOffset();
  
protected :

  int fOddEven;

};

/* ===== OverscanNone ========================================*/

/* this is non-operating overscan, to keep the client code clean */

class OverscanNone : public OverscanSnifs {

public:  
  OverscanNone(){}  
  virtual ~OverscanNone(){}

  virtual void SetImage(ImageSnifs* Image) ;
  virtual void ComputeLines(double * values, double *var);
  virtual void Substract(double* values, double * var);
  virtual double OverscanNoise();
  virtual void OddEvenCorrect();

  // Correct has the same interface as OverscanSnifs

};

/* ===== OverscanFromData ========================================*/

/* this is an overscan estimator from the backgroud of the image itself, 
in case we do not have an overscan region. 
This is very dangerous to use if there is not enough pixel values near the true '0' */

class OverscanFromData : public OverscanSnifs {

public:  
  OverscanFromData() {}
  
  virtual ~OverscanFromData(){}

  virtual void SetImage(ImageSnifs* Image) ;
  virtual void ComputeLines(double * values, double *var);
  virtual double OverscanNoise();
  virtual void OddEvenCorrect();
};

#endif
