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

#ifndef IOMETHOD_H
#define IOMETHOD_H
/*
 IoMethod : the virtual interface to ios
    - wrappers to IFU methods
    - wrappers for methods without i/o
    Only the virtual methods, i.e. the ones that may be overloaded 
      are included here 
*/   

#include "IFU_io.h"
class ImageSimple;


/* ===== IO METHOD T ============================== */
enum IoMethod_t {kIoPlain, kIoSlice, kIoPlainDouble};
const int kIoAll = -1; // flag to open all lines in Slice mode

/* ===== IO METHOD ============================== */

class IoMethod {
  public :
    int Nx() const { return fFrame->nx;}
    int Ny() const { return fFrame->ny;}
    char* Name() const { return Frame()->name; }

    // Wrappers to IMAGE2D methods (and very simple methods)
    virtual void WrFrame(int line, int col, double value) =0;
    // not const because of derived ImageSlice class
    virtual double RdFrame(int line, int col) =0;
    virtual int OpenFrame(char *name, char *mode="Input") =0;
    virtual int CloseFrame() =0;
    virtual int DeleteFrame() =0;
    // for create_frame from an existing image, consider copy contructor
    virtual int CreateFrame(char *name,int nx, int ny, short Type ) =0;
    virtual int CreateFrameFull(char *name,int *Npix,double*start, double*step, short Type, char* a, char* b ) =0;
    virtual void SetMaxLines(int NLines) {};
  
    virtual IoMethod_t IoType() const =0;
  

protected :
  friend class ImageSimple;
  IMAGE2D * Frame() const {return fFrame;}
  IMAGE2D * fFrame;

};

#endif
