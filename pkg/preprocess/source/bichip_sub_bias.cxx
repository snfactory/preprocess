
#include "IFU_io.h"
#include "bichip.hxx"

int main(int argc, char **argv) {

  char **argval, **arglabel;
  
  set_arglist("-in none -out none -bias null");
  init_session(argv,argc,&arglabel,&argval);

  BiChipSnifs * in=new BiChipSnifs(argval[0]);
  BiChipSnifs * out= new BiChipSnifs(*in,argval[1],FLOAT,1);
  
  if is_set(argval[2]) 
    {
      BiChipSnifs * bias=new BiChipSnifs(argval[2]);
      out->SubstractBias(bias);
      delete bias;
    }
  
  delete in;
  delete out;

  exit_session(0);
  
}
