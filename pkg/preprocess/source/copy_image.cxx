
#include "image.hxx"

int main(int argc, char **argv) {

  char **argval, **arglabel;

  set_arglist("-in none -out none");
  init_session(argv,argc,&arglabel,&argval);

  ImageSnifs *in = new ImageSnifs(argval[0]);
  ImageSnifs *out = new ImageSnifs(*in,argval[1]);
  
  delete in;
  delete out;

  exit_session(0);
  
}
