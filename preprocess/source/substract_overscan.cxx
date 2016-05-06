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

#include "imagesnifs.hxx"
#include "preprocessor.hxx"
#include "overscan.hxx"

int main(int argc, char **argv) {

  char **argval, **arglabel;
  
  set_arglist("-in none -out none");
  init_session(argv,argc,&arglabel,&argval);

  ImageSnifs *in = new ImageSnifs(argval[0]);
  ImageSnifs *out = new ImageSnifs(*in,argval[1],FLOAT,1);
  
  Preprocessor P;
  out->HackFitsKeywords();
  P.Overscan()->Correct(out);

  delete in;
  delete out;
  

  exit_session(0);
  
}
