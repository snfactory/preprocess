#include "imagesnifs.hxx"

int main(int argc, char **argv) {

  char **argval, **arglabel;
  
  set_arglist("-in1 none -in2 none -out none -sc1 1.0 -sc2 1.0 -add -divide");
  init_session(argv,argc,&arglabel,&argval);

  double sc1,sc2;
  get_argval(3,"%lf",&sc1);
  get_argval(4,"%lf",&sc2);

  if (is_true(argval[5]) && is_true(argval[6])){
    print_error("Must choose between add and divide options");
    exit_session(1);
  }
  
  ImageSnifs *in1 = new ImageSnifs(argval[0],"I");
  ImageSnifs *in2 = new ImageSnifs(argval[1],"I");
  ImageSnifs *out=new ImageSnifs(*in1,argval[2],0,1);
  
  if (is_true(argval[6]))
    out->Divide(in2);
  else { // add 
    out->Scale(sc1);
    out->Add(in2,sc2);
  }

  delete out;
  delete in1;
  delete in2;
  
  exit_session(0);
}
