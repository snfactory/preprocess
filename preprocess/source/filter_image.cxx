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
#include "imagesnifs.hxx"
#include "filter.hxx"
#include "catorfile.hxx"
#include "utils.h"

/* ----- main ---------------------------------------- */

int main(int argc, char **argv) {

  char **argval, **arglabel;
  int wx,wy;
  double *args;
  ImageFilter *F;

  set_arglist("-in none -out none -halfwindow 2,2 -med|hf|max|sclip|laplacian|pix|chi2 -args 0");
  // remember : in pix mode, halfwindow is the pixel size
  init_session(argv,argc,&arglabel,&argval);


  char inName[lg_name+1],outName[lg_name+1];

  CatOrFile inCat(argval[0]);
  CatOrFile outCat(argval[1]);

  sscanf(argval[2],"%d,%d",&wx,&wy);

  if (!strcmp(arglabel[3],"-hf")) {
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
  } else if (!strcmp(arglabel[3],"-med")) {
    F = new ImageFilterMedian(wx,wy,ImageFilter::kShrinks);
  } else if (!strcmp(arglabel[3],"-chi2")) {
    F = new ImageFilterChi2(wx,wy,ImageFilter::kShrinks);
  } else if (!strcmp(arglabel[3],"-max")) {
    F = new ImageFilterMax(wx,wy,ImageFilter::kShrinks);
  } else if (!strcmp(arglabel[3],"-sclip")) {
    double args;
    int nargs = sscanf(argval[4],"%lf",&args);
    if (nargs<1) {
      print_error("filter_image : expects -args <double> option");
      exit_session(1);
    }
    ImageFilterSigmaClip * FSclip = new ImageFilterSigmaClip(wx,wy,ImageFilter::kShrinks);
    FSclip->SetSigma(args);
    F=FSclip;
  } else if (!strcmp(arglabel[3],"-pix")) {
    double args;
    int nargs = sscanf(argval[4],"%lf",&args);
    if (nargs<1) {
      print_error("filter_image : expects -args <double> option");
      exit_session(1);
    }
    ImageFilterSigmaClip * FSclip = new ImageFilterSigmaClip(wx,wy,ImageFilter::kPixelize);
    FSclip->SetSigma(args);
    F=FSclip;
  } else if (!strcmp(arglabel[3],"-laplacian")){
     F = new ImageFilterLaplacian();
  } else {
    print_error("filter_image : unknown method %s",argval[3]);
  }
  

  while (inCat.NextFile(inName) && outCat.NextFile(outName)) {
    print_msg("Filtering %s",inName);
    ImageSnifs *in = new ImageSnifs(inName,"I");
    ImageSnifs *out = new ImageSnifs(*in,outName);
    F->SetInputImage(in->Image());
    F->SetOutputImage(out->Image());
    F->Filter();
    delete in;
    delete out;
  }

  //delete args;

  exit_session(0);
  
}
