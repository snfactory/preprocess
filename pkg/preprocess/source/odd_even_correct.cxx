#include "image.hxx"

int main(int argc, char **argv) {

  char **argval, **arglabel;
  
  set_arglist("-in none -out none");
  init_session(argv,argc,&arglabel,&argval);

  ImageSnifs *in = new ImageSnifs(argval[0]);
  ImageSnifs *out = new ImageSnifs(*in,argval[1],FLOAT,1);
  
  // hack because value stored is not correct
  char bias[lg_name+1];
  int x1,x2,y1,y2;
  out->RdDesc("BIASSEC",CHAR,lg_name+1,bias);
  sscanf(bias,"[%d:%d,%d:%d]",&x1,&x2,&y1,&y2);
  sprintf(bias,"[%d:%d,%d:%d]",x1+2,x2,1,out->Ny());
  out->WrDesc("BIASSEC",CHAR,lg_name+1,bias);
  
  out->OddEvenCorrect();
  

  delete in;
  delete out;

  exit_session(0);
  
}
