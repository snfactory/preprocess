/* === Doxygen Comment ======================================= */
/*! 
 * \file          image.cxx
 * \copyright     (c) 2003 CRAL-Observatoire de Lyon
 * \date          Wed Aug  6 17:44:35 2003
 * \author        Emmanuel Gangler <e.gangler@ipnl.in2p3.fr>
 * \version       0.0
 * \brief         
 *                
 * $Id$
 **/
/* =========================================================== */

#include "ioplain.hxx"

/* ##### IO PLAIN ################################################# */


/* ===== constructor/destructor ======================================= */

/* ----- IoPlain -------------------------------------------------- */
IoPlain::IoPlain() {
  fLoaded = false;
  fFrame = new IMAGE2D;
}

/* ----- ~IoPlain ------------------------------------------------- */
IoPlain::~IoPlain(){
  if (fLoaded)
    CloseFrame();
  delete fFrame;
}

/* ===== Wrappers to IMAGE2D methods ================================== */

/* ----- WrFrame ------------------------------------------------------ */
void IoPlain::WrFrame(int line, int col, double value){
  //Wrapper to WR_frame
  WR_frame(Frame(),line,col,value);
}

/* ----- RdFrame ------------------------------------------------------ */
double IoPlain::RdFrame(int line, int col) {
  //Wrapper to RD_frame
  return (double) RD_frame(Frame(),line,col);
}

/* ----- OpenFrame ---------------------------------------------------- */
int IoPlain::OpenFrame(char *name, char *mode){
  // Wrapper to open_frame
  fLoaded=true;
  return open_frame(Frame(),name,mode);
}

/* ----- CloseFrame --------------------------------------------------- */
int IoPlain::CloseFrame(){
  // Wrapper to close_frame
  if (fLoaded) {
    fLoaded=false;
    return close_frame(Frame());
  } else {
    return 0;
  }
}

/* ----- CloseFrame --------------------------------------------------- */
int IoPlain::DeleteFrame(){
  // Wrapper to delete_frame
  fLoaded=false;
  return delete_frame(Frame());
}


/* ----- CreateFrame -------------------------------------------------- */
int IoPlain::CreateFrame(char *Name,int Nx, int Ny, short Type ){
  // this CreateFrame is less powerfull than teh original, but
  // this is because I don't understand all parameters !
  int npix[2];
  double start[2];
  double step[2];

  npix[0]=Nx;
  npix[1]=Ny;
  start[0]=start[1]=step[0]=step[1]=1;
  fLoaded=true;

  return create_frame(Frame(),Name,npix,start,step,Type,"","");
}

/* ----- CreateFrameFull -------------------------------------------------- */
int IoPlain::CreateFrameFull(char *name,int *Npix,double* Start, double* Step, short Type,char* Ident, char* U) {

  fLoaded=true;
  return create_frame(Frame(),name,Npix,Start,Step,Type,Ident,U);
}
