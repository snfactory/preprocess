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
  
  set_arglist("-in none -out none -bias null -dark null -flat null -fast");
  init_session(argv,argc,&arglabel,&argval);

  char inName[lg_name+1],outName[lg_name+1];

  CatOrFile inCat(argval[0]);
  CatOrFile outCat(argval[1]);
  
  BiChipSnifs * bias=0;
  ImageSnifs *dark=0, *flat=0;

  // Load once auxilliary files
  if (is_set(argval[2]))
    bias = new BiChipSnifs(argval[2]);
  if (is_set(argval[3]))
    dark = new ImageSnifs(argval[3]);
  if (is_set(argval[4]))
    flat = new ImageSnifs(argval[4]);

  Preprocessor P;
  if (is_true(argval[5]))
    P.SetFastMode(1);

  while (inCat.NextFile(inName) && outCat.NextFile(outName)) {
    
    print_msg("Processing %s",inName);
    ImageSnifs *out = P.Preprocess(inName,outName,bias,dark,flat);
    delete out;
  }

  if (flat) delete flat;
  if (dark) delete dark;
  if (bias) delete bias;

  exit_session(0);
  
}
