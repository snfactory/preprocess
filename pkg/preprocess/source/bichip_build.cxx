/* === Doxygen Comment ======================================= */
/*! 
 * \file          bichip_build.cxx
 * \copyright     (c) 2003 SNIFS-Supernova Factory Experiment
 * \date          Wed Aug  6 17:44:35 2003
 * \author        Emmanuel Gangler <e.gangler@ipnl.in2p3.fr>
 * \version       0.0
 * \brief         
 *                
 * $Id$
 **/
/* =========================================================== */

/* ----- Preprocess includes ----- */
#include "preprocessor.hxx"
#include "bichip.hxx"

int main(int argc, char **argv) {

  char **argval, **arglabel;
  
  set_arglist("-in none -out none");
  init_session(argv,argc,&arglabel,&argval);
  Preprocessor P;
  
  BiChipSnifs * in = P.BuildRawBiChip(argval[0],argval[1]);
  delete in;
  exit_session(0);
  
}
