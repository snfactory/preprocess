/* === Doxygen Comment ======================================= */
/*! 
 * \file          image.hxx
 * \copyright     (c) 2004 SNIFS Collaboration
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

*/

/* ----- library include ----- */
#include <vector>
using namespace std;

/* ----- gsl include ----- */
#include <gsl/gsl_vector.h>

/* ----- local include ----- */
class Section;
#include "IFU_io.h"
#include "iomethod.hxx"

/* ===== IMAGE SIMPLE ============================== */

class ImageSimple {
  public :
  

    ImageSimple(IoMethod_t Method=kIoPlain,int MParam=0);
    ImageSimple(const ImageSimple &image,char* newname,short newtype = 0,int copydata=0,IoMethod_t Method=kIoPlain,int MParam=0);
    ImageSimple(char* name, char* mode="Input",IoMethod_t Method=kIoPlain,int MParam=0);
  
    virtual ~ImageSimple();

    // Wrappers to IMAGE2D content
    int Nx() const { return fIo->Nx();}
    int Ny() const { return fIo->Ny();}
    char* Name() const { return Io()->Name(); }

    // Wrappers to IMAGE2D methods (and very simple methods)
    void WrFrame(int col, int line, double value) {fIo->WrFrame(col,line,value);};
    double RdFrame(int col, int line) const {return fIo->RdFrame(col,line);}
    int WrDesc(char* Descr, short Type, int NbElements, const void* Values);
    int RdDesc(char* Descr, short Type, int NbElements, void* Values) const;
    int RdIfDesc(char* Descr, short Type, int NbElements, void* Values) const;
    int DeleteDesc(char* Descr);
    int OpenFrame(char *name, char *mode="Input") {return Io()->OpenFrame(name,mode);}
    int CloseFrame() {return Io()->CloseFrame();}
    int DeleteFrame() {return Io()->DeleteFrame();}
    // for create_frame from an existing image, consider copy contructor
    int CreateFrame(char *name,int nx, int ny, short Type=FLOAT ) {return Io()->CreateFrame(name,nx,ny,Type);}
    int CreateFrameFull(char *name,int *Npix ,double*Start, double*Step, short Type, char* a, char* b )
        {return Io()->CreateFrameFull(name,Npix,Start,Step,Type,a,b);}

    // internal setters and getters
    // Io is not user settable. It is defined at the construction time
    IoMethod* Io() const {return fIo;}
    IoMethod_t IoType() const {return Io()->IoType();}
    void SetNLines(int NLines);

    // Utility Methods
    int Inside(int Xc,int Yc) const;
    void MinMax(Section* Sec, double * min, double * max);
  
    void ImportHeader(ImageSimple * From){ CP_non_std_desc(From->Frame(),Frame());}
    // internal routine which does not import the variance.
    void ImportSectionFrame(ImageSimple * From, Section* Sec, int X1Start, int Y1Start,int XDir=1,int YDir=1,double ZScale=1 );
    void ImportSection(ImageSimple * From, Section* Sec, int X1Start, int Y1Start,int XDir=1,int YDir=1,double ZScale=1 );
    void ImportFlip(ImageSimple * From, int XDir=1,int YDir=1,double ZScale=1 );
    void Add(ImageSimple* ToAdd, double Scale=1);
    void Add(double Constant);
    void Scale(double Scale);
    void SetTo(double Value);
    void CutLow(double LowBound);
    int AbsThreshold(double Threshold);
    void Divide(ImageSimple* Denom);
    double MeanValue(Section* Sec,int step=1);
    double GetSignificance(int nstep=1);
    void Mask(ImageSimple* in, double cut, int margin, double replace);

    // variance settings
    void SetVarianceFrame(ImageSimple* Var){fVariance = Var;}
    ImageSimple* Variance() const {return fVariance;}  

    // Algorithms

    // Overscan
    double SectionRms(Section* Sec,double sigcut=0);
    gsl_vector* FitBy(vector<ImageSimple*> Refs,double cut=3.0);


    void AddPoissonNoise();
    int HandleSaturation(double Level);
    int CleanWith(ImageSimple * Ref, double SigCut);

  // because we need to hack the keywords
    IMAGE2D * Frame() const {return Io()->Frame();}
  protected :
  
    ImageSimple* fVariance;
  //    static const int nLinesDefault; // parameter for the overscan
    IoMethod *fIo;
  
  private :
    void SetIoMethod(IoMethod_t Method, int MParam);

};


#endif
