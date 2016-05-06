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
//#include "filter.hxx"
#include "catorfile.hxx"
//#include "utils.h"

/* ----- main ---------------------------------------- */

int main(int argc, char **argv) {

  char **argval, **arglabel;
  double cut;

  set_arglist("-in none -ref none -cut 3.0");

  init_session(argv,argc,&arglabel,&argval);

  char inName[lg_name+1];

  CatOrFile inCat(argval[0]);
  ImageSnifs * ref=new ImageSnifs(argval[1]);

  cut = atof(argval[2]);
  
  while (inCat.NextFile(inName)) {
    print_msg("Filtering %s",inName);
    ImageSnifs *in = new ImageSnifs(inName,"IO");
    int nbad = in->CleanWith(ref,cut);
    print_msg("Removed %d pixels",nbad);
    delete in;

  }

  delete ref;
  exit_session(0);
  
}
