
#include "bichip.hxx"
#include "image.hxx"

int main(int argc, char **argv) {

  char **argval, **arglabel;
  
  set_arglist("-in none -out none -bias null -dark null");
  init_session(argv,argc,&arglabel,&argval);

  char tmp_name[lg_name+1];
  sprintf(tmp_name,"mem://%s",argval[1]);

  BiChipSnifs * in=new BiChipSnifs(argval[0]);
  BiChipSnifs * tmp_out= new BiChipSnifs(*in,tmp_name,FLOAT,1);
  BiChipSnifs * bias=0;
  ImageSnifs *dark=0;
  
  if (is_set(argval[2]))
    bias = new BiChipSnifs(argval[2]);
  if (is_set(argval[3]))
    dark = new ImageSnifs(argval[3]);

  ImageSnifs *out = tmp_out->PreprocessFlat(argval[1],bias,dark);

  delete in;
  delete tmp_out;
  delete out;
  if (dark) delete dark;
  if (bias) delete bias;
  

  exit_session(0);
  
}
