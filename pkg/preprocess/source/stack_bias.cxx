/* or more generically : bichip stack gauss */

#include "bichip.hxx"
#include "imagestacksnifs.hxx"
#include "catorfile.hxx"
#include "kombinator.hxx"

int main(int argc, char **argv) {

  char **argval, **arglabel;
  
  set_arglist("-in none -out none -sigcut 3.0 -nlines 1");
  init_session(argv,argc,&arglabel,&argval);

  CatOrFile catIn(argval[0]);
  int nlines;
  get_argval(3,"%d",&nlines);
  
  BiChipStackSnifs * in=new BiChipStackSnifs(&catIn,"I",nlines);
  double sigma;
  get_argval(2,"%lf", &sigma);
  KGauss k(sigma);
  BiChipSnifs* out = in->Kombine(argval[1],&k);
  delete out;
  
  delete in;

  exit_session(0);
  
}
