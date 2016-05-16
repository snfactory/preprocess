/* === Doxygen Comment ======================================= */
/*! 
 * \file          stack_dark.cxx
 * \copyright     (c) 2012 IPNL
 * \date          Wed Aug  6 17:44:35 2003
 * \author        Emmanuel Gangler <e.gangler@ipnl.in2p3.fr>
 * \version       0.0
 * \brief         
 *                
 * $Id$
 **/
/* =========================================================== */

/* this will stack together the outputs of dark_fit */

#include "imagestacksnifs.hxx"
#include "catorfile.hxx"

int main(int argc, char **argv) {

  char **argval, **arglabel;
  
  set_arglist("-in none -out none -nlines 4100");
  init_session(argv,argc,&arglabel,&argval);

  CatOrFile catIn(argval[0]);
  int nlines;
  get_argval(2,"%d",&nlines);
  
  ImageStackStack * in=new ImageStackStack(&catIn,"I",nlines);


  ImageStackSnifs* out = in->Kombine(argval[1]);
  delete out;
  delete in;

  exit_session(0);
  
}
