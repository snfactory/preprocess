/* === Doxygen Comment ======================================= */
/*! 
 * \file          hack_B_headers.cxx
 * \copyright     (c) 2004 SNIFS Collaboration
 * \date          Wed Aug  6 17:44:35 2003
 * \author        $Author$
 * \version       $Name$
 * \brief         
 *                
 * $Id$
 **/
/* =========================================================== */

/* ----- C includes ----- */
#include <stdlib.h>

/* ----- IFU include ----- */
#include "IFU_io.h"

/* ----- local include ----- */
#include "bichip.hxx"
#include "imagesnifs.hxx"
#include "catorfile.hxx"
#include "utils.h"
//#include "preprocessor.hxx"

/* ----- main ---------------------------------------- */

int main(int argc, char **argv) {

  char **argval, **arglabel;
  
  set_arglist("-in none -out none");
  init_session(argv,argc,&arglabel,&argval);

  char inName[lg_name+1],outName[lg_name+1];

  CatOrFile inCat(argval[0]);
  CatOrFile outCat(argval[1]);
  
  while (inCat.NextFile(inName) && outCat.NextFile(outName)) {
    
    print_msg("Processing %s",inName);
    char command[lg_name*2+1];
    sprintf(command,"cp %s %s",inName, outName);
    system(command);
    BiChipSnifs *out = new BiChipSnifs(outName,"IO");
    char primary_name[lg_name+1];
    ut_primary_header_name(out->Chip(0)->Name(),primary_name);
    Anyfile primary_header;
    //    open_primary_hd(&primary_header,primary_name,"IO");
    open_primary_hd(&primary_header,primary_name,"I");
    for (int i=0;i<out->NChips();i++)
      

      CP_non_std_desc(&primary_header,out->Chip(i)->Frame());
      

    delete out;
  }

  exit_session(0);
  
}
