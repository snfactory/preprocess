/* === Doxygen Comment ======================================= */
/*! 
 * \file          imagestack.hxx
 * \copyright     (c) 2003 CRAL-Observatoire de Lyon
 * \date          Wed Aug  6 18:32:01 2003
 * \author        Emmanuel Gangler <e.gangler@ipnl.in2p3.fr>
 * \version       0.0
 * \brief         
 *                
 * $Id$
 **/
/* =========================================================== */

#ifndef IMAGESTACKSNIFS_H
#define IMAGESTACKSNIFS_H
/*   
     ImageStackSnifs is the SNIFS processor for ImageStack.
*/

#include <vector>
using namespace std;

class CatOrFile;
class Kombinator;
class KombinatorFit;
class ValueGetter;

#include "imagesnifs.hxx"
#include "bichip.hxx"

/* ===== ImageStackSnifs ======================================== */

class ImageStackSnifs  {
  public :

    ImageStackSnifs(int NLines=1){fNLinesMem = NLines;}
    ImageStackSnifs(CatOrFile * Cat, char * Mode="I", int NLines = 5000);
    ~ImageStackSnifs();
  
    ImageSnifs* Kombine(char *outName,Kombinator* k);
    void AddImage(ImageSnifs* Image) {Image->SetNLines(fNLinesMem);fImages.push_back(Image);}
  

  protected :
    vector<ImageSnifs *> fImages;
    int fNLinesMem;

};

/* ===== BiChipStackSnifs ======================================== */

class BiChipStackSnifs  {
  public :

    BiChipStackSnifs(int NLines=1) {fNLinesMem=NLines;}
    BiChipStackSnifs(CatOrFile * Cat, char * Mode="I",int NLines=1);
    ~BiChipStackSnifs();

    void AddBiChip(BiChipSnifs* BiChip) {BiChip->SetNLines(fNLinesMem);fBiChips.push_back(BiChip);}
  
  //    BiChipStackSnifs* PreprocessBias(CatOrFile * Out,int NLines=1); 
  //  ImageStackSnifs* PreprocessDark(CatOrFile * Out, BiChipSnifs * Bias, int NLines=1);
    BiChipSnifs* Kombine(char* OutName,Kombinator* K);
    void KombineFit(BiChipSnifs** out, char** OutName,KombinatorFit* K,ValueGetter *V);
  //BiChipSnifs* MakeBiasFrame(CatOrFile * tmpOut, char* BiasName,double SigCut);

  protected :
    vector<BiChipSnifs *> fBiChips;
    // number of resident lines in memory when the whole image is not needed
    int fNLinesMem;

};

#endif
