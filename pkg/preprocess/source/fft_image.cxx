/* === Doxygen Comment ======================================= */
/*! 
 * \file          rootify.cxx
 * \copyright     (c) 2003 SNIFS-Supernova Factory Experiment
 * \date          Wed Aug  6 17:44:35 2003
 * \author        Emmanuel Gangler <e.gangler@ipnl.in2p3.fr>
 * \version       0.0
 * \brief         
 *                
 * $Id$
 **/
/* =========================================================== */

/* ----- image includes ------------------------------ */
#include "image.hxx"
#include "catorfile.hxx"
#include "section.hxx"
#include "analyser.hxx"

/*------------ main ----------------------*/
int main(int argc, char **argv) {

  char **argval, **arglabel;
  ImageSnifs* in,*out;
  char inName[lg_name+1], outName[lg_name+1];

  set_arglist("-in none -out none -sec [4:1027,1:4102]");
  init_session(argv,argc,&arglabel,&argval);

  CatOrFile inCat(argval[0]);
  CatOrFile outCat(argval[1]);
  Section sec(argval[2]);

  /* loop on catalog */
  while(inCat.NextFile(inName) && outCat.NextFile(outName)) {
    printf("Opening now %s\n",inName);
    in = new ImageSnifs(inName);
    ImageAnalyser ana(in,&sec);

    out = ana.LineFft(outName);

    delete in;
    delete out;
  }
  
  exit_session(0);

  return(0);
}



