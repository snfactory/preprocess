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

/* This applies teh bias and dark removal to an already preprocessed file
no checks are preformed if the preprocessed file is already bias- or dark-bubtracted*/

/* ----- std includes ----- */
#include <stdlib.h>
#include <vector>
using namespace std;

/* ----- local include ----- */
#include "imagesnifs.hxx"
#include "catorfile.hxx"
#include "darkmodel.hxx"

/* ----- main ---------------------------------------- */

int main(int argc, char **argv) {

  char **argval, **arglabel;
  char inName[lg_name+1],outName[lg_name+1];

  set_arglist("-in none -out none -bm null -dm null -bias null -timeon null");
  init_session(argv,argc,&arglabel,&argval);

  CatOrFile inCat(argval[0]);
  CatOrFile outCat(argval[1]);

  DarkModel *biasModel=0;
  DarkModel *darkModel=0;
  if (is_set(argval[2]))
    biasModel = new DarkModel(argval[2]);
  if (is_set(argval[3]))
    darkModel = new DarkModel(argval[3]);

  ImageSnifs *bias=0;
  ImageSnifs *dark=0;
  if (is_set(argval[4]))
    bias = new ImageSnifs(argval[4]);



  while (inCat.NextFile(inName) && outCat.NextFile(outName)) {
    ImageSnifs *in = new ImageSnifs(inName,"I");
    ImageSnifs *out = new ImageSnifs(*in,outName,0,1);

    if (is_set(argval[3])) 
      out->WrDesc("TIMEON",CHAR,lg_name+1,argval[3]);

    if (biasModel)
      out->SubstractBiasModel(biasModel);
    if (bias) 
      out->SubstractBias(bias);
    if (darkModel)
      out->SubstractDarkModel(darkModel);
    if (dark) 
      out->SubstractDark(dark);

    delete in;
    delete out;


  }
  exit_session(0);
  
}
