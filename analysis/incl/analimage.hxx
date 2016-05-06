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

/* ----- standard includes ------------------------------ */
#include <vector>
using namespace std;

/* ----- ROOT includes ------------------------------ */
#include "Rtypes.h"
#include "TObject.h"
class TTree;
class TFile;

/* ----- image includes ------------------------------ */
class ImageSnifs;

/* ----- IFU includes ------------------------------ */
// rootcint does not cooperate with IFU_io.h
const int ImageSignatureLineLength = 256;
// to placate dict.cxx
#ifndef lg_name
#define lg_name 80L
#endif

/* ----- local includes ------------------------------ */
//#include "signaturecut.hxx"
class SignatureCut;



/* ===== PrintInfo ============================== */

class PrintInfo {
public:
  virtual ~PrintInfo() {}
  virtual void FillHelp(char* HelpString)=0;
  virtual void FillLine(char* LineString)=0;
  virtual void FillName(char* NameString)=0;
  virtual void FillValue(char* ValueString)=0;
  virtual char* GetName()=0;
  virtual double GetValue()=0;
  
};

template <class T> class PrintInfoType : public PrintInfo {
public :

  PrintInfoType(const char* Help, const char* Name, const char* Format, int length, T* toprint);
  virtual void FillHelp(char* HelpString);
  virtual void FillLine(char* LineString);
  virtual void FillName(char* NameString);
  virtual void FillValue(char* ValueString);
  
  virtual char* GetName() {return fOrgName;}
  virtual double GetValue();
  // {return *((double*) fInfo);}
  

protected:
  char fHelpInfo[lg_name+1];
  char fOrgName[lg_name+1];
  char fName[lg_name+1];
  char fFormat[lg_name+1];
  unsigned int fLength;
  T* fInfo;
  
};

template<> class PrintInfoType<char> : public PrintInfoType<int> {
public :

  PrintInfoType(const char* Help, const char* Name, const char* Format, int length, char* toprint);
  virtual void FillValue(char* ValueString);
  virtual double GetValue();  

  //protected:
  char* fInfo; // hides the typed <int> fInfo
  
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
    // cuts
    void ParseCutsFile(char* file);
    SignatureCut* ParseCutsLine(char* line);
    void ApplyCuts();
    void ResetCuts();
  

  protected :
  
  //friend class AnalImageSignature;

  // for custom print
  static const int fkNitems;
  PrintInfo** fPrintInfo; //!
  int* fIsAlive;
  // for selection
  vector<SignatureCut *> fSigCuts; //!

  char fName[80];
  int fFClass;
  int fExternal;
  int fChannel;
  double fExpTime;
  double fJulDate;

  double fRobustMean;
  double fRobustRMS;
  int fNUnder;
  int fNOver;
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

