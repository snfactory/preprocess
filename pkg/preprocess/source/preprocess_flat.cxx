
#include "preprocessor.hxx"
#include "bichip.hxx"
#include "imagesnifs.hxx"
#include "catorfile.hxx"

int main(int argc, char **argv) {

  char **argval, **arglabel;
  
  set_arglist("-in none -out none -bias null -dark null");
  init_session(argv,argc,&arglabel,&argval);

  char tmp_name[lg_name+1],inName[lg_name+1],outName[lg_name+1];

  CatOrFile inCat(argval[0]);
  CatOrFile outCat(argval[1]);
  
  BiChipSnifs * bias=0;
  ImageSnifs *dark=0;
  
  // Load once auxilliary files
  if (is_set(argval[2]))
    bias = new BiChipSnifs(argval[2]);
  if (is_set(argval[3]))
    dark = new ImageSnifs(argval[3]);

  Preprocessor P;
  
  while (inCat.NextFile(inName) && outCat.NextFile(outName)) {

    ImageSnifs* out = P.PreprocessFlat(inName,outName,bias,dark);
    delete out;
  }

  if (dark) delete dark;
  if (bias) delete bias;

  exit_session(0);
  
}
