
#include "bichip.hxx"
#include "image.hxx"
#include "section.hxx"

int main(int argc, char **argv) {

  char **argval, **arglabel;
  
  set_arglist("-in none -out none -model none");
  init_session(argv,argc,&arglabel,&argval);

  // in shall be a preprocessed snifs image
  ImageSnifs *in = new ImageSnifs(argval[0]);
  ImageSnifs *model = new ImageSnifs(argval[2]);
  ImageSnifs *out = new ImageSnifs(*model,argval[1],FLOAT,0);

  int naxis1, naxis2;
  in->RdDesc("NAXIS1",INT,1,&naxis1);
  in->RdDesc("NAXIS2",INT,1,&naxis2);
  Section sec(1,naxis1,1,naxis2);
  int out_naxis1, out_naxis2;
  in->RdDesc("NAXIS1",INT,1,&out_naxis1);
  in->RdDesc("NAXIS2",INT,1,&out_naxis2);
  out->ImportSectionFrame(in,&sec,1,1+out_naxis2-naxis2);
  out->Variance()->ImportSectionFrame(in,&sec,1+out_naxis1-naxis1,1);
  // the Readout noise
  out->Variance()->Add(5);

  delete model;
  delete in;
  delete out;
  

  exit_session(0);
  
}
