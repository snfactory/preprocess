
#include "image.hxx"

int main(int argc, char **argv) {

  char **argval, **arglabel;

  set_arglist("-in none -out none -nlines 100");
  
  init_session(argv,argc,&arglabel,&argval);

  int nlines;
  get_argval(2,"%d",&nlines);
  ImageSimple *in = new ImageSimple(argval[0],"I",kIoSlice,nlines);
  ImageSimple *out;
  out = new ImageSimple(*in,argval[1],0,1,kIoSlice,nlines);
  /*
out = new ImageSimple(kIoSlice,nlines);
  int npix[2];
  double start[2];
  double step[2];

  npix[0] = 1056;
  npix[1] = 4128;
  start[0] = 0;
  start[1] = 0;
  step[0] = 1;
  step[1] = 1;
  out->CreateFrameFull(argval[1],npix,start,step,FLOAT,"id","un");
  */
  delete out;
  delete in;

  exit_session(0);
  
}
