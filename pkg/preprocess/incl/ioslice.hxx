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

#ifndef IOSLICE_H
#define IOSLICE_H
/*
 IoSlice : the image is only loaded for at most nlines
    - wrappers to IFU methods
    - wrappers for methods without i/o
    Only the virtual methods, i.e. the ones that may be overloaded 
      are included here 
*/   

#include "ioplain.hxx"

/* ===== IO SLICE ============================== */

class IoSlice : public IoPlain{
  public :

    IoSlice(int MaxLines);
    virtual ~IoSlice();

    // overloaded functions
    virtual void WrFrame(int line, int col, double value);
    virtual double RdFrame(int line, int col); 
    virtual int OpenFrame(char *name, char *mode="Input");
    virtual int CloseFrame();
    virtual int CreateFrame(char *name,int nx, int ny, short Type=FLOAT );
    virtual int CreateFrameFull(char *name,int *Npix,double* Start, double* Step, short Type,char* Ident, char* U);
    virtual IoMethod_t IoType() const {return kIoSlice;}


    // Specific functions (not in the interface now)
    void SetMaxLines(int Max);
    void LoadFromLine(int Line, int Read=1);
    void WriteCurrent();
  
    int FirstLine() const {return fFirstLine;}

    int MaxLine() const {return fMaxLines;}
    int LastLine() const {return fLastLine;}
  

  protected :
      // user not allowed to call directly these
    int MemAlloc();
    int MemFree();
    int MemIsAlloc();
    void LoadLineRange(int First, int Last);
  

    int fFirstLine; // first line in the data region
    int fLastLine; // first line in the data region
    int fMaxLines; // maximum length (in lines) allowed to be stored
    int fAllocLines; // number of lines really allocated
    int* fExists; // marks the possibility to load line from fitsfile
                  // (means LoadFromLine will accept them)
  
};


#endif
