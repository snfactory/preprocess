#include "image.hxx"

int main(int argc, char **argv) {

  char **argval, **arglabel;
  
  set_arglist("-in none -out none -dark none");
  init_session(argv,argc,&arglabel,&argval);

  ImageSnifs *in = new ImageSnifs(argval[0]);
  ImageSnifs *out = new ImageSnifs(*in,argval[1],FLOAT,1);
  ImageSnifs *dark = new ImageSnifs(argval[2]);
  
  out->SubstractDark(dark);

  delete in;
  delete out;
  delete dark;

  exit_session(0);
  
}
