/* === Doxygen Comment ======================================= */
/*! 
 * \file          image.cxx
 * \copyright     (c) 2004 SNIFS Collaboration
 * \date          Wed Aug  6 17:44:35 2003
 * \author        Emmanuel Gangler <e.gangler@ipnl.in2p3.fr>
 * \version       0.0
 * \brief         
 *                
 * $Id$
 **/
/* =========================================================== */

/*----- IFU include -----*/
#include "IFU_io.h"

/* ----- local includes  -----*/
#include "bichip.hxx"
#include "preprocessor.hxx"
#include "overscan.hxx"

int main(int argc, char **argv) {

  char **argval, **arglabel;
  
  set_arglist("-in none -out none");
  init_session(argv,argc,&arglabel,&argval);

  BiChipSnifs * in=new BiChipSnifs(argval[0]);
  BiChipSnifs * out= new BiChipSnifs(*in,argval[1],FLOAT,1);
  
  out->HackFitsKeywords();
  
  //out->OddEvenCorrect();
  out->CreateVarianceFrame();
  Preprocessor P;
  P.SetOverscanAuto(out->Chip(0));
  P.Overscan()->Correct(out);
  
  delete in;
  delete out;

  exit_session(0);
  
}
