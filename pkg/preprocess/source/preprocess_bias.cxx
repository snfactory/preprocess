
#include "bichip.hxx"
#include "imagestacksnifs.hxx"
#include "catorfile.hxx"

int main(int argc, char **argv) {

  char **argval, **arglabel;
  
  set_arglist("-in none -out none -stackout null -sigcut 3.0 -nlines 1");
  init_session(argv,argc,&arglabel,&argval);

  CatOrFile catIn(argval[0]);
  CatOrFile catOut(argval[1]);
  int nlines;
  get_argval(4,"%d",&nlines);

  BiChipStackSnifs * in=new BiChipStackSnifs(&catIn,"I",nlines);
  BiChipStackSnifs * tmp= in->PreprocessBias(&catOut,nlines);
  if (is_set(argval[3])) {
    double sigma;
    get_argval(3,"%lf", &sigma);
    BiChipSnifs* out = tmp->KombineGauss(argval[2],sigma);
    delete out;
  }
  
  delete in;
  delete tmp;

  exit_session(0);
  
}
