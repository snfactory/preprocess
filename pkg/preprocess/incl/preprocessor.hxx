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
#include "algocams.hxx"
#include "iomethod.hxx"

/*   
 Preprocessor is a steering for the preprocessing.
 Ideally, all should be plugged in there, but this is still evolutive ...

 Consider it as a collection of mains.

 The only implementation yet is the otcom/detcom sniffer
*/

/* ===== PREPROCESSOR ======================================== */

class Preprocessor {
  public :
   
    // Constructors/Destructors
    Preprocessor( );
  
    ~Preprocessor() {}

    // setters
    void SetIoMethod(IoMethod_t mode)  {fMode = mode;  }

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
  
  
};


#endif
