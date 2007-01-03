/* === Doxygen Comment ======================================= */
/*! 
 * \file          filter_imaged.cxx
 * \copyright     (c) 2004 SNIFS Collaboration
 * \date          Wed Aug  6 17:44:35 2003
 * \author        Emmanuel Gangler <e.gangler@ipnl.in2p3.fr>
 * \version       0.0
 * \brief         
 *                
 * $Id$
 **/
/* =========================================================== */

/* ----- std includes ----- */
#include <stdlib.h>

/* ----- local include ----- */
#include "imagesnifs.hxx"
#include "filter.hxx"
#include "catorfile.hxx"
#include "utils.h"

/* ----- main ---------------------------------------- */

int main(int argc, char **argv) {

  char **argval, **arglabel;
  double cut;
  ImageFilterRemCosmic *F;

  set_arglist("-in none -out none -cut 2.0");
  // cut=2.0 preserves all arcs
  // cut=1.5 preserves most of the arcs - especially the R channel is losing
  // cut=1.0 is recomended for removing most of the R cosmics
  init_session(argv,argc,&arglabel,&argval);

  char inName[lg_name+1],outName[lg_name+1];

  CatOrFile inCat(argval[0]);
  CatOrFile outCat(argval[1]);
  cut = atof(argval[2]);
  
  F = new ImageFilterRemCosmic();

  while (inCat.NextFile(inName) && outCat.NextFile(outName)) {
    print_msg("Filtering %s",inName);
    ImageSnifs *in = new ImageSnifs(inName,"I");
    ImageSnifs *out = new ImageSnifs(*in,outName,FLOAT,1);
    F->SetInputImage(in->Image());
    F->SetOutputImage(out->Image());
    F->SetRatio(cut);
    F->Filter();
    delete in;
    delete out;
  }


  exit_session(0);
  
}
