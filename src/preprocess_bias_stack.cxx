#include "bichipstacksnifs.hxx"
#include "bichip.hxx"

int main(int argc, char **argv) {

  char **argval, **arglabel;
  
  set_arglist("-in none -out none -stackout null -sigcut 3.0");
  init_session(argv,argc,&arglabel,&argval);

  BiChipStackSnifs * in=new BiChipStackSnifs(argval[0]);
  BiChipStackSnifs * tmp= in->PreprocessBias(argval[1]);
  if (is_set(argval[2])) {
    double sigcut;
    get_argval(3,"%f",&sigcut)
    BiChipSnifs* out = tmp->KombineGauss(argval[2],sigcut);
  }
  
  delete in;
  delete tmp;
  delete out;

  exit_session(0);
  
}
