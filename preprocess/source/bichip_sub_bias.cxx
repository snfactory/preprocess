
#include "IFU_io.h"
#include "bichip.hxx"
#include "catorfile.hxx"

int main(int argc, char **argv) {

  char **argval, **arglabel;
  char inName[lg_name+1],outName[lg_name+1];

  set_arglist("-in none -out none -bias null");
  init_session(argv,argc,&arglabel,&argval);

  CatOrFile inCat(argval[0]);
  CatOrFile outCat(argval[1]);

  BiChipSnifs * bias=0;
  if is_set(argval[2])
      bias=new BiChipSnifs(argval[2]);
    
  while (inCat.NextFile(inName) && outCat.NextFile(outName)) {
  
    BiChipSnifs * in=new BiChipSnifs(inName);
    BiChipSnifs * out =  new BiChipSnifs(*in,outName,FLOAT,1);
    if (bias)
      out->SubstractBias(bias);    

    delete out;
    delete in;
  }

  delete bias;

  exit_session(0);
  
}
