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

#ifndef IOPLAIN_H
#define IOPLAIN_H
/*
 IoPlain : the image is entirely loaded
    - wrappers to IFU methods
    - wrappers for methods without i/o
    Only the virtual methods, i.e. the ones that may be overloaded 
      are included here 
*/   

#include "iomethod.hxx"

/* ===== IO PLAIN ============================== */

class IoPlain : public IoMethod{
  public :

    IoPlain();
    virtual ~IoPlain();

    // Wrappers to IMAGE2D methods (and very simple methods)
    virtual void WrFrame(int line, int col, double value);
    // not const because of derived ImageSlice class
    virtual double RdFrame(int line, int col); 
    virtual int OpenFrame(char *name, char *mode="Input");
    virtual int CloseFrame();
    virtual int DeleteFrame();
    // for create_frame from an existing image, consider image copy contructor
    virtual int CreateFrame(char *name,int nx, int ny, short Type=FLOAT );
    virtual int CreateFrameFull(char *name,int *Npix,double* Start, double* Step, short Type,char* Ident, char* U);

    virtual IoMethod_t IoType() const {return kIoPlain;}


  protected :
  bool fLoaded;

};


#endif
