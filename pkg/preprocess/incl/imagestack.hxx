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
#include <vector>
using namespace std;

/* ===== ImageStack ======================================== */

class ImageStack  {
  public :

    ImageStack();
    ImageStack(vector<ImageSimple*> Images);
    ImageStack(CatOrFile * Cat, char * Mode="I");
    ~ImageStack();
  
    void SetKombinator(Kombinator * Kombinator) {fKombinator = Kombinator;}
    Kombinator * GetKombinator() {return fKombinator;}
    void Kombine(ImageSimple *ToFill, int FillsVarOut=1, int UpdateInitialVar=0);
    int Nx(){return fNx;}
    int Ny(){return fNy;}
  

  protected :
    // protected because Nx and Ny shall be correctly set
    vector<ImageSimple *> * List() {return &fImageList;}

    vector<ImageSimple *> fImageList; // the stack owns the images is fIsOwner is set.
    Kombinator * fKombinator;
    int fIsOwner;
    int fNx;
    int fNy;

};

#endif
