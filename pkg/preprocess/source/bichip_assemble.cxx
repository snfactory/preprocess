
//#include "IFU_io.h"
#include "bichip.hxx"
#include "imagesnifs.hxx"

int main(int argc, char **argv) {

  char **argval, **arglabel;
  
  set_arglist("-in none -out none");
  init_session(argv,argc,&arglabel,&argval);

  // IO because HackGainRatio puts a key in there
  BiChipSnifs * in=new BiChipSnifs(argval[0],"IO");
  in->HackGainRatio();
  
  ImageSnifs *out = in->Assemble(argval[1]);

  delete in;
  delete out;

  exit_session(0);
  
}
