#include "image.hxx"

int main(int argc, char **argv) {

  char **argval, **arglabel;
  
  set_arglist("-in none -out none -flat none");
  init_session(argv,argc,&arglabel,&argval);

  ImageSnifs *in = new ImageSnifs(argval[0]);
  ImageSnifs *out = new ImageSnifs(*in,argval[1],FLOAT,1);
  ImageSnifs *flat = new ImageSnifs(argval[2]);
  
  out->ApplyFlat(flat);

  delete in;
  delete out;
  delete flat;

  exit_session(0);
  
}
