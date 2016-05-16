/* === Doxygen Comment ======================================= */
/*! 
 * \file          preprocess_bias.cxx
 * \copyright     (c) 2004 SNIFS Collaboration
 * \date          Wed Aug  6 17:44:35 2003
 * \author        Emmanuel Gangler <e.gangler@ipnl.in2p3.fr>
 * \version       0.0
 * \brief         
 *                
 * $Id$
 **/
/* =========================================================== */

#include "imagestacksnifs.hxx"
#include "catorfile.hxx"
#include "kombinator.hxx"
#include "preprocessor.hxx"

int main(int argc, char **argv) {

  char **argval, **arglabel;
  
  set_arglist("-in none -out none -stackout null -sigcut 3.0 -nlines 5000 -all");
  init_session(argv,argc,&arglabel,&argval);

  CatOrFile catIn(argval[0]);
  CatOrFile catOut(argval[1]);
  int nlines;
  char inName[lg_name+1],outName[lg_name+1];
  Preprocessor P;
  if (is_true(argval[5]))
    P.SetAllImage(1);

  get_argval(4,"%d",&nlines);
  
  if (is_set(argval[2])) {
    
    ImageStackSnifs * outStack = new ImageStackSnifs(nlines);
  
    while (catIn.NextFile(inName) && catOut.NextFile(outName)) {
      print_msg("Opening %s",inName);
      P.SetIoMethod(kIoSlice);
      // Don't want to add Poisson Noise on it.
      ImageSnifs * tmpOut = P.PreprocessAssemble(inName,outName);
      tmpOut->HandleCosmetics();
      outStack->AddImage(tmpOut);
    }
    
    double sigma;
    get_argval(3,"%lf", &sigma);
    KGauss k(sigma);
    ImageSnifs* out = outStack->Kombine(argval[2],&k);

    delete outStack;
    delete out;
  
  } else { // no need to store temporary information
  
    while (catIn.NextFile(inName) && catOut.NextFile(outName)) {
      print_msg("Opening %s",inName);
      ImageSnifs * out = P.PreprocessAssemble(inName,outName);
      delete out;  
    }
  }

  exit_session(0);
}
