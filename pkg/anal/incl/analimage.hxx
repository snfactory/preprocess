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

const int ImageSignatureLineLength = 256;

/* ===== PrintInfo ============================== */

class PrintInfo {
public:
  virtual void FillHelp(char* HelpString)=0;
  virtual void FillLine(char* LineString)=0;
  virtual void FillName(char* NameString)=0;
  virtual void FillValue(char* ValueString)=0;
};

template <class T> class PrintInfoType : public PrintInfo {
public :

  PrintInfoType(const char* Help, const char* Name, const char* Format, int length, T* toprint);
  void FillHelp(char* HelpString);
  void FillLine(char* LineString);
  void FillName(char* NameString);
  void FillValue(char* ValueString);
  
protected:
  char fHelpInfo[lg_name+1];
  char fName[lg_name+1];
  char fFormat[lg_name+1];
  unsigned int fLength;
  T* fInfo;
  
};



/* ===== IMAGE SIGNATURE ============================== */

class ImageSignature : public TObject {
  public :
  
    ImageSignature();
    ~ImageSignature();
    void Reset();
    void FillExternal(int External); 
    void Fill(ImageSnifs* Preprocessed);
    void SetPrintItems(char* ItemsDescr) ;
    void PrintBlurb();
    void PrintHeader();
    void PrintContent();
    void PrintTrailer();

  protected :
  
  //friend class AnalImageSignature;

  // for custom print
  static const int fkNitems;
  PrintInfo** fPrintInfo;
  int* fIsAlive;

  char fName[80];
  int fFClass;
  int fExternal;
  int fChannel;
  double fExpTime;
  double fJulDate;

  double fRobustMean;
  double fRobustRMS;
  int fNOutliers;
  int fSatu;
  
  double fQuant99;
  double fQuant999; 
  double fQuant9999;
  double fFocus;
  //  double fOverscanLevel;
  double fNoise;
  double fGainMed;
  double fGainQuant;

public:
  ClassDef(ImageSignature,1)

};

/* ===== STEERING IMAGE SIGNATURE ============================== */

// The image shall be preprocessed ...

class AnalImageSignature {
  public :
  
    AnalImageSignature(TFile* Rfile=0);
  
    ~AnalImageSignature();

  //   void SetImage(ImageSnifs* Image);
  //  ImageSnifs* Image(){return fImage;}
    ImageSignature * Signature() {return fSignature;}
  

  void MakeTree();
  //  void ImageType(char * TypeString); // Image type name
  // void FillImageSignature(ImageSnifs* Preprocessed);
  //void ApplyCutsToImageSignature();
    void StoreImageSignature();
  
  protected :
  
  ImageSnifs * fImage; //image
  ImageSignature* fSignature; //! image signature
  //ImageCuts* fCuts;
  TFile* fFile;
  TTree* fTree;
};

#endif

