/* === Doxygen Comment ======================================= */
/*! 
 * \file          rootify_data.cxx
 * \copyright     (c) 2003 SNIFS-Supernova Factory Experiment
 * \date          Wed Aug  6 17:44:35 2003
 * \author        Emmanuel Gangler <e.gangler@ipnl.in2p3.fr>
 * \version       0.0
 * \brief         
 *                
 * $Id$
 **/
/* =========================================================== */

/* ----- ROOT includes ------------------------------ */

/* ----- image includes ------------------------------ */
#include "imagesnifs.hxx"
#include "catorfile.hxx"

/* ----- local includes ------------------------------ */
#include "analimage.hxx"


/*------------ main ----------------------*/
int main(int argc, char **argv) {

  char **argval, **arglabel;
  char inName[lg_name+1];
  ImageSignature sig;

  set_arglist("-in none -blurb");
  init_session(argv,argc,&arglabel,&argval);

  CatOrFile inCat(argval[0]);

  if (is_true(argval[1]))
    sig.PrintBlurb();
    sig.PrintHeader();
    
  /* loop on catalog */
  while(inCat.NextFile(inName)) {
    ImageSnifs *in = new ImageSnifs(inName);
    
    sig.Reset();
    sig.Fill(in);
    sig.PrintContent();
    
    delete in;
  }

  if (is_true(argval[1]))
    sig.PrintTrailer();
             

  exit_session(0);
  return(0);
}
