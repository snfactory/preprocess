
#include "bichip.hxx"
#include "imagestacksnifs.hxx"
#include "catorfile.hxx"
#include "kombinator.hxx"

int main(int argc, char **argv) {

  char **argval, **arglabel;
  int nlines;
  
  set_arglist("-in none -out none -stackout null -bias null -sigma 4.0 -nlines 5000");
  init_session(argv,argc,&arglabel,&argval);

  CatOrFile catIn(argval[0]);
  CatOrFile catOut(argval[1]);

  BiChipSnifs * bias=0;  
  // Load once auxilliary files
  if (is_set(argval[3]))
    bias = new BiChipSnifs(argval[3]);

  get_argval(4,"%d",&nlines);

  BiChipStackSnifs * in=new BiChipStackSnifs(&catIn,"I",nlines);
  ImageStackSnifs * tmp= in->PreprocessDark(&catOut,bias,nlines);
  delete in;

  if (is_set(argval[2])) {
    double sigma;
    get_argval(4,"%lf", &sigma);
    Kombinator * k = new KGaussPoisson(sigma);
    ImageSnifs *out = tmp->Kombine(argval[2],k);
    delete k;
    delete out;
  }
  delete tmp;
  if (bias) 
    delete bias;

  exit_session(0);
}
