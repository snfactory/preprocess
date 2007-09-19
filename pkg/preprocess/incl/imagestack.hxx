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

#ifndef IMAGESTACK_H
#define IMAGESTACK_H
/*   
     Image stack is a base for handling a stack of images.
     The only implemented purpose is now the combination of identical
     images. But more can come in the future.
*/

class ImageSimple;
class CatOrFile;
class Kombinator;
class KombinatorFit;

#include <vector>
using namespace std;

#include "valuegetter.hxx"

/* ===== ImageStack ======================================== */

class ImageStack  {
  public :

    ImageStack();
    ImageStack(vector<ImageSimple*> Images);
    ImageStack(CatOrFile * Cat, char * Mode="I");
    ~ImageStack();
  
    // "normal" kombinator
    void SetKombinator(Kombinator * Kombinator) {fKombinator = Kombinator;}
    Kombinator * GetKombinator() {return fKombinator;}
    void Kombine(ImageSimple *ToFill, int FillsVarOut=1, int UpdateInitialVar=0);
    // kombinator with fit added
    void SetKombinatorFit(KombinatorFit * Kombinator) {fKombinatorFit = Kombinator;}
    KombinatorFit * GetKombinatorFit() {return fKombinatorFit;}
    void KombineFit(ImageSimple **ToFill, int FillsVarOut=1, int UpdateInitialVar=0);
    // and the trickiest getvalue (part of kombinatorfit)
    double GetValue(ImageSimple* image) 
       {return fValueGetter->GetValue(image);};
    void SetValueGetter(ValueGetter * Getter) {fValueGetter = Getter;}

    int Nx(){return fNx;}
    int Ny(){return fNy;}
    vector<ImageSimple *> * List() {return &fImageList;}
  

  protected :
    // protected because Nx and Ny shall be correctly set

    vector<ImageSimple *> fImageList; // the stack owns the images is fIsOwner is set.
    Kombinator * fKombinator;
    KombinatorFit * fKombinatorFit;
    ValueGetter * fValueGetter;

    int fIsOwner;
    int fNx;
    int fNy;

};

#endif
