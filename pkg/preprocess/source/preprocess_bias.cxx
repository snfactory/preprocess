
#include "bichip.hxx"
#include "image.hxx"

int main(int argc, char **argv) {

  char **argval, **arglabel;
  
  set_arglist("-in none -out none");
  init_session(argv,argc,&arglabel,&argval);

  BiChipSnifs * in=new BiChipSnifs(argval[0]);
  BiChipSnifs * out= new BiChipSnifs(*in,argval[1],FLOAT,1);

  out->PreprocessBias();

  delete in;
  delete out;

  exit_session(0);
  
}
