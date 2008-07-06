/* === Doxygen Comment ======================================= */
/*! 
 * \file          preprocess_dark.cxx
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
#include "darkmodel.hxx"

/* ----- main ---------------------------------------- */

int main(int argc, char **argv) {

  char **argval, **arglabel;
  
  set_arglist("-in none -out none -bias null -fast -all -bm null -dm null");
  init_session(argv,argc,&arglabel,&argval);

  char inName[lg_name+1],outName[lg_name+1];

  CatOrFile inCat(argval[0]);
  CatOrFile outCat(argval[1]);
  
  ImageSnifs * bias=0;
  ImageSnifs *dark=0, *flat=0;

  // Load once auxilliary files
  if (is_set(argval[2]))
    bias = new ImageSnifs(argval[2]);

  Preprocessor P;
  if (is_true(argval[3]))
    P.SetFastMode(1);
  if (is_true(argval[4]))
    P.SetAllImage(1);

  DarkModel *biasModel=0;
  DarkModel *darkModel=0;
  if (is_set(argval[5]))
    biasModel = new DarkModel(argval[7]);
  if (is_set(argval[6]))
    darkModel = new DarkModel(argval[8]);

  while (inCat.NextFile(inName) && outCat.NextFile(outName)) {
    
    print_msg("Processing %s",inName);
    ImageSnifs *out = P.PreprocessDark(inName,outName,bias,biasModel,darkModel);
    delete out;
  }

  if (flat) delete flat;
  if (dark) delete dark;
  if (bias) delete bias;

  exit_session(0);
  
}
