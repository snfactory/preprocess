/* or more generically : bichip stack gauss */

#include "bichip.hxx"
#include "imagestacksnifs.hxx"
#include "catorfile.hxx"
#include "kombinatorfit.hxx"
#include "valuegetter.hxx"
#include "section.hxx"

int main(int argc, char **argv) {

  char **argval, **arglabel;
  
  set_arglist("-in none -out none -sigcut 5.0 -nlines 5000 " );
  init_session(argv,argc,&arglabel,&argval);

  CatOrFile catIn(argval[0]);
  int nlines;
  get_argval(3,"%d",&nlines);

  double sigma;
  get_argval(2,"%lf", &sigma);
  
  KombinatorFit *k = new KFGain(sigma);
  Section S("[4:1027,1:4102]");
  ValueGetter *vg = new ValueAnalyserMean(&S);
  
  BiChipStackSnifs * in=new BiChipStackSnifs(&catIn,"I",nlines);
  BiChipSnifs* out[1];
  in->KombineFit(out,&argval[1],k,vg);
  delete out[0];
  
  delete vg;
  delete in;
  delete k;
  
  exit_session(0);
  
}
