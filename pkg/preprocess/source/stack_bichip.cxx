/* or more generically : bichip stack gauss */

#include "bichip.hxx"
#include "imagestacksnifs.hxx"
#include "catorfile.hxx"
#include "kombinator.hxx"

int main(int argc, char **argv) {

  char **argval, **arglabel;
  
  set_arglist("-in none -out none -sigcut 3.0 -nlines 5000 -method Gauss|Spread" );
  init_session(argv,argc,&arglabel,&argval);

  CatOrFile catIn(argval[0]);
  int nlines;
  get_argval(3,"%d",&nlines);
  

  double sigma;
  get_argval(2,"%lf", &sigma);
  
  Kombinator * k;
  if (argval[4][0]=='S')
    k = new KSpread(sigma,1);
  else
    k = new KGauss(sigma);

  BiChipStackSnifs * in=new BiChipStackSnifs(&catIn,"I",nlines);
  BiChipSnifs* out = in->Kombine(argval[1],k);
  delete out;
  
  delete in;
  delete k;
  
  exit_session(0);
  
}
