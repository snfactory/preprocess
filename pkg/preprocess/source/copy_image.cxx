
#include "image.hxx"

int main(int argc, char **argv) {

  char **argval, **arglabel;

  set_arglist("-in none -out none");
  init_session(argv,argc,&arglabel,&argval);

  ImageSimple *in = new ImageSimple(argval[0]);
  ImageSimple *out = new ImageSimple(*in,argval[1],0,1);
  
  delete in;
  delete out;

  exit_session(0);
  
}
