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

/* ----- local include ----- */
#include "image.hxx"
#include "filter.hxx"
#include "catorfile.hxx"
#include "utils.h"

/* ----- main ---------------------------------------- */

int main(int argc, char **argv) {

  char **argval, **arglabel;
  int wx,wy;
  double *args;
  ImageFilter *F;

  set_arglist("-in none -out none -halfwindow 2,2 -method hf -args 0");
  init_session(argv,argc,&arglabel,&argval);

  char inName[lg_name+1],outName[lg_name+1];

  CatOrFile inCat(argval[0]);
  CatOrFile outCat(argval[1]);

  sscanf(argval[2],"%d,%d",&wx,&wy);

  if (!strcmp(argval[3],"hf")) {
    args = new double[2];
    int nargs = sscanf(argval[4],"%lf,%lf",args,args+1);
    if (nargs<2) {
      print_error("filter_image : expects -args <double>,<double> option");
      exit_session(1);
    }
    ImageFilterHF * FHF= new ImageFilterHF(wx,wy,ImageFilter::kShrinks);
    FHF->SetNoDataAnswer(1,ut_big_value);
    FHF->SetThreshold(args[0]);
    FHF->SetSignificance(args[1]);
    F= FHF;
  } else {
    print_error("filter_image : unknown method %s",argval[3]);
  }
  

  while (inCat.NextFile(inName) && outCat.NextFile(outName)) {
    print_msg("Filtering %s",inName);
    ImageSimple *in = new ImageSimple(inName,"I");
    ImageSimple *out = new ImageSimple(*in,outName);
    F->SetInputImage(in);
    F->SetOutputImage(out);
    F->Filter();
    delete in;
    delete out;
  }

  delete args;

  exit_session(0);
  
}