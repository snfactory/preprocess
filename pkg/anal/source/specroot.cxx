#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*--- iolib include ---*/
//#include "gendef.h"
//#include "items.h"
//#include "iofuncdecl.h"

/*--- mathlib include */
#include "gsl/gsl_statistics.h"

/*--- Ccd/preprocess includes ---*/
#include "image.hxx"

/*--- Specific ROOT includes */
#include "TString.h"
#include "TFile.h"
#include "TMatrix.h"
#include "TProfile.h"
#include "TH1F.h"

#define MAX_LINE 8192

/* --------------- FillHistoLine ------------------- */
void FillHistoLine(ImageSnifs* image, int npix) {

  TH1F * histo = new TH1F("line","Pixel values",npix,1,npix+1);
  
  for (int ipix=0;ipix<npix;ipix++) {
    int nx =  image->Nx();
    int i = ipix % nx;
    int j = ipix / nx;
    histo->SetBinContent(ipix+1,image->RdFrame(i,j));
  }
  histo->Write();
  delete histo;
}




/*------------ main ----------------------*/
int main(int argc, char **argv) {

  char **argval, **arglabel;
  char im_name[lg_name+1],out_name[lg_name+1];
  int npix;
  

  set_arglist("-in none -out none -number 21120");
  init_session(argv,argc,&arglabel,&argval);

  strcpy(im_name,argval[0]);
  strcpy(out_name,argval[1]);
  get_argval(2,"%d", &npix);

  ImageSnifs * image = new ImageSnifs(im_name);
  TFile * file = new TFile(out_name,"RECREATE");

  FillHistoLine(image,npix);

  file->Close();
  delete file;
  delete image;

  exit_session(0);
  return 0 ;
}



