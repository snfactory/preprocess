/* or more generically : bichip stack gauss */

#include "bichip.hxx"
#include "imagestacksnifs.hxx"
#include "catorfile.hxx"
#include "kombinator.hxx"

int main(int argc, char **argv) {

  char **argval, **arglabel;
  
  set_arglist("-in none -out none -sigcut 3.0 -nlines 4100 -method Average|Gauss|Spread|Median");
  init_session(argv,argc,&arglabel,&argval);

  CatOrFile catIn(argval[0]);
  int nlines;
  get_argval(3,"%d",&nlines);
  
  ImageStackSnifs * in=new ImageStackSnifs(&catIn,"I",nlines);
  double sigma;
  get_argval(2,"%lf", &sigma);

  Kombinator * k;
  if (argval[4][0]=='S')
    k = new KSpread(sigma,1);
  else if (argval[4][0]=='A')
    k = new KSpread(sigma,0);
  else if (argval[4][0]=='M')
    k= new KMedian();
  else
    k = new KGauss(sigma);

  ImageSnifs* out = in->Kombine(argval[1],k);
  delete out;
  
  delete in;

  exit_session(0);
  
}
