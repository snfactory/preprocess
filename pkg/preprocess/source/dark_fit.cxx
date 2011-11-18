/* This will take an input catalog and perform a dark fit to the residuals */

#include "imagestacksnifs.hxx"
#include "catorfile.hxx"
#include "kombinatorfit.hxx"
#include "darkmodel.hxx"
#include "valuegetter.hxx"

int main(int argc, char **argv) {

  char **argval, **arglabel;
  
  set_arglist("-in none -out none -sigcut 3.0 -nlines 4100 -darkmodel none");
  init_session(argv,argc,&arglabel,&argval);


  CatOrFile catIn(argval[0]);

  double sigma;
  get_argval(2,"%lf", &sigma);

  int nlines;
  get_argval(3,"%d",&nlines);  
  ImageStackSnifs * in=new ImageStackSnifs(&catIn,"I",nlines);

  DarkModel darkModel(argval[4]);

  ValuesGetter * vg = new ValuesGetterDarkFitter(&darkModel);
  KombinatorFitND * k = new KFLinear(sigma,vg->NParams());

  ImageStackSnifs* out = in->KombineFitND(argval[1],k,vg);
  delete out;
  
  delete in;

  exit_session(0);
  
}
