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

#ifndef IMAGE_H
#define IMAGE_H
/*
 ImageSimple : methods that do not rely on precise fits keywords
   This contains 
    - wrappers to IFU methods
    - wrappers for methods without i/o
    - basic image manipulation
   
 ImageSnifs : methods which are compliant with te SNIFS header
   This contains
    - explicit tasks for the reduction
    - all checks for paranoiac mode
    
*/

#include "IFU_io.h"
class Section;


/* ===== IMAGE SIMPLE ============================== */

class ImageSimple {
  public :

    ImageSimple();
    ImageSimple(const ImageSimple &image,char* newname,short newtype = 0,int copydata=0);
    ImageSimple(char* name, char* mode="Input");
  
    ~ImageSimple();

    // Wrappers to IMAGE2D content
    int Nx() const { return fFrame->nx;}
    int Ny() const { return fFrame->ny;}
    char* Name() const { return fFrame->name; }

    // Wrappers to IMAGE2D methods (and very simple methods)
    void WrFrame(int line, int col, double value);
    void WrFrame(int line, int col, long value);
    double RdFrame(int line, int col) const;
    int WrDesc(char* Descr, short Type, int NbElements, void* Values);
    int RdDesc(char* Descr, short Type, int NbElements, void* Values) const;
    int RdIfDesc(char* Descr, short Type, int NbElements, void* Values) const;
    int DeleteDesc(char* Descr);
    int OpenFrame(char *name, char *mode="Input");
    int CloseFrame();
    int DeleteFrame();
    // for create_frame from an existing image, consider copy contructor
    int CreateFrame(char *name,int nx, int ny, short Type=FLOAT );

    // Utility Methods
    int Inside(int Xc,int Yc);
    void ImportHeader(ImageSimple * From){ CP_non_std_desc(From->fFrame,fFrame);}
    void ImportSectionFrame(ImageSimple * From, Section* Sec, int X1Start, int Y1Start,int XDir=1,int YDir=1,double ZScale=1 );
    void ImportSection(ImageSimple * From, Section* Sec, int X1Start, int Y1Start,int XDir=1,int YDir=1,double ZScale=1 );
    void Add(ImageSimple* ToAdd, double Scale=1);
    void Add(double Constant);
    void Scale(double Scale);
    void SetTo(double Value);
    void CutLow(double LowBound);
    void Divide(ImageSimple* Denom);
    double MeanValue(Section* Sec,int step=1) const;

    // variance settings
    void SetVarianceFrame(ImageSimple* Var){fVariance = Var;}
    ImageSimple* Variance() const {return fVariance;}  

    // Algorithms

    // Overscan
    void SubstractOverscan(Section* Sec);
    void ComputeLinesOverscan(Section* Sec,double * values);
    double OverscanRms(Section* Sec,double sigcut=0);

    // odd-even
    void OddEvenCorrect(Section* Sec,double*param, double sigcut=0);
    //
    void AddPoissonNoise();
    void HandleSaturation(double Level);
  
  
  protected :
    IMAGE2D *fFrame;
    bool fLoaded;
    ImageSimple* fVariance;
    static const int nLinesDefault; // parameter for the overscan

};



/* ===== IMAGE SNIFS ======================================== */

class ImageSnifs : public ImageSimple {
  public :
   
    // Constructors/Destructors

    ImageSnifs();
    ImageSnifs(const ImageSnifs &image,char* newname,short newtype = 0,int copydata=0);
    ImageSnifs(char* name, char* mode="Input");
    ~ImageSnifs();

    // Tools
    ImageSnifs* BuildSubImage(Section* Sec,char* Name);

    // Algorithms
    
    // in-place overscan substraction
    void SubstractOverscan();
    // in-place odd-even substrction
    void OddEvenCorrect();

    void SubstractBias(ImageSnifs* Bias);
    void SubstractDark(ImageSnifs* Dark);
    void ApplyFlat(ImageSnifs* Flat);

    void BuildFlat();
    void AddPoissonNoise();
    void HandleSaturation();
  
    // Processing management
    void CreateVarianceFrame(char* name="");
    void AddOverscanVariance();

    // photometric tools
    void FilterSectionFill(Section * Sec,int Filter);
    void SplitMultiFilter(char* NameRecipee, ImageSnifs** Subs);
  

    // Hacks

    void HackFitsKeywords();

    // Utilities

    bool ParanoMode() const { return fParano; }
    void SetParanoMode(bool Parano) { fParano=Parano; }
  

    protected :
    bool fParano;
    
};


#endif
