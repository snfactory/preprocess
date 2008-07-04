/* === Doxygen Comment ======================================= */
/*! 
 * \file          remove_dark.cxx
 * \copyright     (c) 2007 SNIFS Collaboration
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
#include <vector>
using namespace std;

/* ----- local include ----- */
//#include <gsl/gsl_vector.h>
//#include <gsl/gsl_matrix.h>
//#include <gsl/gsl_linalg.h>


/* ----- local include ----- */
#include "imagesnifs.hxx"
//#include "imagestacksnifs.hxx"
//#include "filter.hxx"
#include "catorfile.hxx"
#include "darkmodel.hxx"
//#include "kombinator.hxx"
//#include "analyser.hxx"
//#include "section.hxx"
//#include "utils.h"

/* ----- main ---------------------------------------- */

int main(int argc, char **argv) {

  char **argval, **arglabel;
  char inName[lg_name+1],outName[lg_name+1];

  set_arglist("-in none -out none -model none -timeon null");
  init_session(argv,argc,&arglabel,&argval);

  CatOrFile inCat(argval[0]);
  CatOrFile outCat(argval[1]);
  DarkModel Model(argval[2]);


  while (inCat.NextFile(inName) && outCat.NextFile(outName)) {
    ImageSnifs *in = new ImageSnifs(inName,"I");
    ImageSnifs *out = new ImageSnifs(*in,outName,0,1);

  if (is_set(argval[3])) 
    out->WrDesc("TIMEON",CHAR,lg_name+1,argval[3]);
  else {
    char* noneStr="None";
    out->WrDesc("TIMEON",CHAR,lg_name+1,noneStr);
  }
    out->SubstractBiasModel(&Model);

    delete in;
    delete out;


  }
  exit_session(0);
  
}
