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
class ImageSnifs;
class BiChipSnifs;

/* ===== ImageStackSnifs ======================================== */

class ImageStackSnifs  {
  public :

    ImageStackSnifs(CatOrFile * Cat, char * Mode="I");

  protected :
    vector<ImageSnifs *> fImages;

};

/* ===== BiChipStackSnifs ======================================== */

class BiChipStackSnifs  {
  public :

    BiChipStackSnifs(int NLines=1) {fNLinesMem=NLines;}
    BiChipStackSnifs(CatOrFile * Cat, char * Mode="I",int NLines=1);
    ~BiChipStackSnifs();

    void AddBiChip(BiChipSnifs* BiChip) {BiChip->SetNLines(fNLinesMem);fBiChips.push_back(BiChip);}
  
    BiChipStackSnifs* PreprocessBias(CatOrFile * tmpOut,int NLines=1); 
    BiChipSnifs* KombineGauss(char* OutName,double SigCut);
    BiChipSnifs* MakeBiasFrame(CatOrFile * tmpOut, char* BiasName,double SigCut);

  protected :
    vector<BiChipSnifs *> fBiChips;
    // number of resident lines in memory when the whole image is not needed
    int fNLinesMem;

};

#endif
