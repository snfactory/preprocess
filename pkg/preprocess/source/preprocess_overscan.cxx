/* === Doxygen Comment ======================================= */
/*! 
 * \file          preprocess.cxx
 * \copyright     (c) 2004 SNIFS Collaboration
 * \date          Wed Aug  6 17:44:35 2003
 * \author        Emmanuel Gangler <e.gangler@ipnl.in2p3.fr>
 * \version       0.0
 * \brief         
 *                
 * $Id$
 **/
/* =========================================================== */

/* ----- local include ----- */
#include "bichip.hxx"
#include "imagesnifs.hxx"
#include "catorfile.hxx"
#include "preprocessor.hxx"

/* ----- main ---------------------------------------- */

int main(int argc, char **argv) {

  char **argval, **arglabel;
  
  set_arglist("-in none -out none -fast -all");
  init_session(argv,argc,&arglabel,&argval);

  char inName[lg_name+1],outName[lg_name+1];

  CatOrFile inCat(argval[0]);
  CatOrFile outCat(argval[1]);
  
  Preprocessor P;
  if (is_true(argval[2]))
    P.SetFastMode(1);
  if (is_true(argval[3]))
    P.SetAllImage(1);

  while (inCat.NextFile(inName) && outCat.NextFile(outName)) {
    
    print_msg("Processing %s",inName);
    BiChipSnifs *out = P.PreprocessOverscan(inName,outName);
    delete out;
  }

  exit_session(0);
  
}
