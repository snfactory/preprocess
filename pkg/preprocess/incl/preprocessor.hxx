/* === Doxygen Comment ======================================= */
/*! 
 * \file          preprocessor.hxx
 * \copyright     (c) 2003 CRAL-Observatoire de Lyon
 * \date          Wed Aug  6 18:32:01 2003
 * \author        Emmanuel Gangler <e.gangler@ipnl.in2p3.fr>
 * \version       0.0
 * \brief         
 *                
 * $Id$
 **/
/* =========================================================== */

#ifndef PREPROCESSOR_H
#define PREPROCESSOR_H

/* ----- Local includes and definitions ----- */
class BiChipSnifs;
class ImageSnifs;
class OverscanSnifs;

#include "algocams.hxx"
#include "iomethod.hxx"

/*   
 Preprocessor is a steering for the preprocessing.
 Ideally, all should be plugged in there, but this is still evolutive ...
 This is of course a SNIFS preprocessor.
(i.e. if you have non-SNIFS images, you have to rewrite...)

 Besides a collection of mains, there is some utilities, and algorithm switches.

 The only implementation yet is the otcom/detcom sniffer
*/

/* ===== PREPROCESSOR ======================================== */

class Preprocessor {
  public :
   
    // Constructors/Destructors
    Preprocessor( );
  
  ~Preprocessor();

    // setters and getters
    void SetIoMethod(IoMethod_t mode)  {fMode = mode;  }
    OverscanSnifs* Overscan() {return fOverscan;}
    void SetOverscanAuto(ImageSnifs* Image);
    void SetFastMode(int IsFast) {fFast = IsFast;}
    int FastMode() {return fFast;}
    void SetAllImage(int IsAll) {fAllImage = IsAll;}
    int GetAllImage() {return fAllImage;}
  

    // method

    BiChipSnifs* BuildRawBiChip(char* name,char* outName="");
    BiChipSnifs* PreprocessBias(char* name,char* outName="");

    ImageSnifs * PreprocessAssemble(char* name, char* outName, BiChipSnifs* bias);
    ImageSnifs * PreprocessDark(char* name, char* outName,BiChipSnifs* bias);
    ImageSnifs * PreprocessFlat(char* name, char* outName,BiChipSnifs* bias, ImageSnifs* dark);
    ImageSnifs* Preprocess(char* name, char* OutName,BiChipSnifs *bias,ImageSnifs *dark,ImageSnifs* flat) ;
  
    
  
  
  protected:
    
    BiChipSnifs* fBiChip;
    ImageSnifs* fImage;
    IoMethod_t fMode;
    OverscanSnifs * fOverscan;
    OverscanSnifs * fOverscanSnifs;
    OverscanSnifs * fOverscanRescue;
    int fFast; // if =1 : only basic algorithms applied
    int fAllImage; // =1 => all image =0 => only 2 first chips
};


#endif
