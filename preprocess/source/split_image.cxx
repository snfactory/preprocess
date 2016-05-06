
#include "bichip.hxx"
#include "imagesnifs.hxx"
#include "section.hxx"

int main(int argc, char **argv) {

  char **argval, **arglabel;

  set_arglist("-in none -out none");
  init_session(argv,argc,&arglabel,&argval);

  // in shall be a preprocessed snifs image
  // out shall be an arg name recipee with a %d
  // which will be translated to 'a', 'b', ...
  ImageSnifs *in = new ImageSnifs(argval[0]);
  ImageSnifs *out[5];
  
  in->SplitMultiFilter(argval[1],out);
  
  delete in;
  for (int i=0;i<5;i++)
    delete out[i];

  exit_session(0);
  
}
