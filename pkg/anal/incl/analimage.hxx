/* === Doxygen Comment ======================================= */
/*! 
 * \file          analimage.hxx
 * \copyright     (c)  2003 SNIFS-Supernova Factory Experiment
 * \date          Wed Aug  6 18:32:01 2003
 * \author        Emmanuel Gangler <e.gangler@ipnl.in2p3.fr>
 * \version       0.0
 * \brief         
 *                
 * $Id$
 **/
/* =========================================================== */

#ifndef ANALIMAGE_H
#define ANALIMAGE_H

/* ----- ROOT includes ------------------------------ */
#include "Rtypes.h"
#include "TObject.h"
class TTree;
class TFile;

/* ----- image includes ------------------------------ */
class ImageSnifs;


/* ===== IMAGE SIGNATURE ============================== */

class ImageSignature : public TObject {
  public :
  
    ImageSignature() {}

  protected :
  
  friend class AnalImageSignature;

  int fImageType;
  double fExpTime;
  double fJulDate;
  double fRobustMean;
  double fRobustRMS;
  double fQuant99;
  double fQuant999; 
  int fSatu;
  double fFocus;
  double fOverscanLevel;
  double fNoise;
  double fGain;
  double fColor;

public:
  ClassDef(ImageSignature,1)



};

/* ===== STEERING IMAGE SIGNATURE ============================== */

class AnalImageSignature {
  public :
  
    AnalImageSignature(TFile* Rfile=0);
  
    ~AnalImageSignature();

    void SetImage(ImageSnifs* Image);

    void ImageType(char * TypeString); // Image type name
    void FillImageSignature();
  //void ApplyCutsToImageSignature();
    void StoreImageSignature();
  
  protected :
  
  ImageSnifs * fImage; //image
  ImageSignature* fSignature; // image signature
  //ImageCuts* fCuts;
  TFile* fFile;
  TTree* fTree;
};

#endif

