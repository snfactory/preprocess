/* === Doxygen Comment ======================================= */
/*! 
 * \file          image.hxx
 * \copyright     (c) 2003 CRAL-Observatoire de Lyon
 * \date          Wed Aug  6 18:32:01 2003
 * \author        Emmanuel Gangler <e.gangler@ipnl.in2p3.fr>
 * \version       0.0
 * \brief         
 *                
 * $Id$
 **/
/* =========================================================== */

#ifndef IMAGESNIFS_H
#define IMAGESNIFS_H
/*   
 ImageSnifs : methods which are compliant with te SNIFS header
   This contains
    - explicit tasks for the reduction
    - all checks for paranoiac mode
 The dependency with respect to the image is a 'Has a', in order to 
 allow for multiple image opening styles.   
*/

/* ----- local includes and definitions ----- */
#include "image.hxx"
class AlgoCams;
class DarkModel;

/* ===== CHANNEL ======================================== */
/* the right comparison shall not be channel = kBlue
  but rather channel & kBlue*/

enum Channel_t {kUnknown=0, kBlueChannel = 1, kRedChannel = 2, kPhotometric = 4, kGuiding = 8};
//const int kNChannel = 4;

/* ===== IMAGE SNIFS ======================================== */

class ImageSnifs :public ImageSimple {
  public :
   
    // Constructors/Destructors

    ImageSnifs(IoMethod_t Method=kIoPlain,int MParam=0);
    ImageSnifs(const ImageSnifs &image,char* newname,short newtype = 0,int copydata=0,IoMethod_t Method=kIoPlain,int MParam=0);
    ImageSnifs(char* name, char* mode="Input",IoMethod_t Method=kIoPlain,int MParam=0);
    ~ImageSnifs();

    // ImageSimple access (kept to allow the switch from 'is a' to 'has a') 
    // if necessary one day
    const ImageSimple * Image() const {return this;}
    ImageSimple * Image() {return this;}
    // ShortCuts
    
  int Nx() const { return Image()->Nx();}
  int Ny() const { return Image()->Ny();}
  char* Name() const { return Image()->Name(); }
  void WrFrame(int line, int col, double value) 
      {Image()->WrFrame(line,col,value);}
  double RdFrame(int line, int col) const {return Image()->RdFrame(line,col);}
  int WrDesc(char* Descr, short Type, int NbElements, const void* Values) 
      {return Image()->WrDesc(Descr,Type,NbElements,Values);}
  int RdDesc(char* Descr, short Type, int NbElements, void* Values) const
      {return Image()->RdDesc(Descr,Type,NbElements,Values);}
  int RdIfDesc(char* Descr, short Type, int NbElements, void* Values) const
      {return Image()->RdIfDesc(Descr,Type,NbElements,Values);}
  int DeleteDesc(char* Descr) {return Image()->DeleteDesc(Descr);}
  int OpenFrame(char *name, char *mode="Input") 
      {return Image()->OpenFrame(name,mode);}
  int CloseFrame() {return Image()->CloseFrame();}
  int DeleteFrame() {return Image()->DeleteFrame();}
  int CreateFrame(char *name,int nx, int ny, short Type=FLOAT ){return Image()->CreateFrame(name,nx,ny,Type);}

  void ImportHeader(ImageSnifs * From)
      { Image()->ImportHeader(From->Image());}
  void ImportSection(ImageSnifs * From, Section* Sec, int X1Start, int Y1Start,int XDir=1,int YDir=1,double ZScale=1 ) 
      { Image()->ImportSection(From->Image(),Sec,X1Start,Y1Start,XDir,YDir,ZScale);}
  void SetVarianceFrame(ImageSimple* Var){Image()->SetVarianceFrame(Var);}
    ImageSimple* Variance() const {return Image()->Variance();}  



    // Tools
    ImageSnifs* BuildSubImage(Section* Sec,char* Name);

    // Algorithms
    
    // in-place overscan substraction
    // void SubstractOverscan();
    // in-place odd-even substrction
    void OddEvenCorrect();

    void UpdateFClass();
    void SubstractBiasModel(DarkModel* Model);
    void SubstractBias(ImageSnifs* Bias);
    void SubstractDarkModel(DarkModel* Model);
    void SubstractDark(ImageSnifs* Dark);
    void ApplyFlat(ImageSnifs* Flat);

    void BuildFlat();
    void BuildDark();
    void AddPoissonNoise();
    void HandleSaturation();
    void HandleCosmetics();
    void CheatCosmetics();
    void SpecialRedCosmetics();
    void CustomFlat();
  

    // Processing management
    void CreateVarianceFrame(char* name="");
    void AddOverscanVariance();
  //    void Assembled2Dark();
  //  void Assembled2Flat(ImageSnifs *dark);
  //  void Assembled2Preprocessed(ImageSnifs *dark,ImageSnifs *flat);
    int CanBeStackedWith(ImageSnifs* image);
    int HasOverscan();
  
  

    // photometric tools
    void FilterSectionFill(Section * Sec,int Filter);
    void SplitMultiFilter(char* NameRecipee, ImageSnifs** Subs);
  

    // Hacks

    void HackFitsKeywords();

    // Utilities and setters/getters

    // Parano mode
    bool ParanoMode() const { return fParano; }
    void SetParanoMode(bool Parano) { fParano=Parano; }
    // FClass
    int GetFClass();
    void SetFClass(int);
    // algorithms
    void SetAlgo(char*); // algo is coded in fits keywords
    AlgoCams* Algo();
    // Channel  
    int GetChannel();
    void SetChannel(int Channel);
    // Sections
    Section * DataSec();
    Section * BiasSec();
    // Other utilities
    double RdNoise();
  
    

  protected :
    // bacause it is too specialized
    int IdenticalPreprocDesc(ImageSnifs * ToCheck, char* Desc, int* val1, int* val2);
  

    bool fParano;
  //ImageSimple* fImage;
    AlgoCams * fAlgo;
    
};


#endif
