#include "image.hxx"

int main(int argc, char **argv) {

  char **argval, **arglabel;
  
  set_arglist("-in1 none -in2 none -out none -sc1 1.0 -sc2 1.0");
  init_session(argv,argc,&arglabel,&argval);

  double sc1,sc2;
  get_argval(3,"%lf",&sc1);
  get_argval(4,"%lf",&sc2);

  ImageSimple *in1 = new ImageSimple(argval[0],"I");
  ImageSimple *in2 = new ImageSimple(argval[1],"I");
  ImageSimple *out=new ImageSimple(*in1,argval[2],0,1);
  out->Scale(sc1);
  out->Add(in2,sc2);
  
  delete out;
  delete in1;
  delete in2;
  
  exit_session(0);
}
