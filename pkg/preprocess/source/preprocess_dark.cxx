
#include "bichip.hxx"
#include "image.hxx"

int main(int argc, char **argv) {

  char **argval, **arglabel;
  
  set_arglist("-in none -out none -bias null");
  init_session(argv,argc,&arglabel,&argval);

  char tmp_name[lg_name+1];
  sprintf(tmp_name,"mem://%s",argval[1]);

  BiChipSnifs * in=new BiChipSnifs(argval[0]);
  BiChipSnifs * tmp_out= new BiChipSnifs(*in,tmp_name,FLOAT,1);
  BiChipSnifs * bias=0;
  
  if (is_set(argval[2]))
    bias = new BiChipSnifs(argval[2]);

  ImageSnifs *out = tmp_out->PreprocessDark(argval[1],bias);

  delete in;
  delete tmp_out;
  delete out;
  if (bias) delete bias;

  exit_session(0);
  
}
