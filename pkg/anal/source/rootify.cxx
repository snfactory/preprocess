#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "gsl/gsl_statistics.h"

# ifdef __cplusplus
extern "C" {
#endif

#include "gendef.h"
#include "items.h"
#include "iofuncdecl.h"

#include "utils.h"

# ifdef __cplusplus
}
#endif



/*--- Specific ROOT includes */
#include "TString.h"
#include "TFile.h"
#include "TMatrix.h"
#include "TProfile.h"
#include "TH1F.h"

#define MAX_LINE 8192

/*-------------- builds a few root usefull histos ------------- */
/* this section depends on ROOT package */
void HorizontalProfile(IMAGE2D* image,int yin, int yout,float sigma){
  
  TProfile * prof = new TProfile("horiprfl","Horizontal Profile",image->nx,0,image->nx);
  
  int i,j,remain;
  double *over, *overptr;
  int nover = yout-yin;

  over = new double[nover];

  for (i=0;i<image->nx;i++) {
    for (j=yin;j<yout;j++) 
      over[j-yin] = RD_frame(image,i,j);

    /* purify the overscan from spurious */
    overptr=over;
    remain = nover;
    if (sigma>0) 
      ut_trunc_sigma_unknown(&overptr,&remain,sigma);
    
    for (j=0;j<remain;j++) {
      prof->Fill(i,overptr[j]);
    }
  }
  prof->Write();
  delete prof;
  delete[]over;
  
}

/*-------------- builds a few root usefull histos ------------- */
/* this section depends on ROOT package */
void VerticalProfile(IMAGE2D* image,int xin, int xout,float sigma){
  
  TProfile * prof = new TProfile("vertprfl","Vertical Profile",image->ny,0,image->ny);
  
  int i,j,remain;
  double over[MAX_LINE],*overptr;
  int nover = xout-xin;

  for (j=0;j<image->ny;j++) {
    for (i=xin;i<xout;i++) 
      over[i-xin] = RD_frame(image,i,j);

    /* purify the overscan from spurious */
    remain=nover;
    overptr=over;
    if (sigma>0) 
     ut_trunc_sigma_unknown(&overptr,&remain,sigma);
    
    for (i=0;i<remain;i++) {
      prof->Fill(j,overptr[i]);
    }
  }
  prof->Write();
  delete prof;
  
}

/*-------------- builds a few root usefull histos ------------- */
/* this section depends on ROOT package */
void OverscanProfile(IMAGE2D* image,int xin, int xout,float sigma){
  
  TProfile * prof = new TProfile("ovscprfl","Overscan Strip",image->ny,0,image->ny);
  
  int i,j,remain;
  double over[MAX_LINE],*overptr;
  int nover = xout-xin;

  for (j=0;j<image->ny;j++) {
    for (i=xin;i<xout;i++) 
      over[i-xin] = RD_frame(image,i,j);

    /* purify the overscan from spurious */
    remain=nover;
    overptr=over;
    if (sigma>0) 
      ut_trunc_sigma_unknown(&overptr,&remain,sigma);

    for (i=0;i<remain;i++) {
      prof->Fill(j,overptr[i]);
    }
  }
  prof->Write();
  delete prof;
  
}

/*-------------- builds a few root usefull histos ------------- */
/* this section depends on ROOT package */
void OddEvenVerticalProfile(IMAGE2D* image,int xin, int xout,float sigma){
  
  TProfile * prof = new TProfile("oevertprfl","Odd-Even Vertical profile Strip",image->ny,0,image->ny);
  TProfile * prof2= new TProfile("oevp2","Odd-Even Vertical profile Strip 2",image->ny,0,image->ny);
  
  int i,j,remainOdd,remainEven,rOE;
  double odd[MAX_LINE],*oddptr;
  double even[MAX_LINE],*evenptr;
  double oddeven[MAX_LINE],*oddevenptr;
  double wodd,weven;

  for (j=0;j<image->ny;j++) {
    for (i=xin;i<xout;i+=2) { /* 'odd' */ 
      odd[(i-xin)/2] = RD_frame(image,i,j);
      if (i+1<xout) {
        even[(i-xin)/2] = RD_frame(image,i+1,j);
        oddeven[(i-xin)/2] = odd[(i-xin)/2]-even[(i-xin)/2];
      }
    }
  
    /* purify from spurious */
    remainOdd=(xout+1-xin)/2;
    rOE = remainEven=(xout-xin)/2;
    oddptr=odd;
    evenptr=even;
    oddevenptr=oddeven;

    if (sigma>0) {
      ut_trunc_sigma_unknown(&oddptr,&remainOdd,sigma);
      ut_trunc_sigma_unknown(&evenptr,&remainEven,sigma);
      ut_trunc_sigma_unknown(&oddevenptr,&rOE,sigma);
    }

    /* wheighting condition : sum w = Ntot, sum wodd=sum weven */
    wodd = (remainOdd+remainEven)/(2.0*remainOdd);
    weven = (remainOdd+remainEven)/(2.0*remainEven);
    
    /* fill the profile */
    for (i=0;i<remainOdd;i++)
      prof->Fill(j,oddptr[i],wodd);
    for (i=0;i<remainEven;i++)
      prof->Fill(j,-evenptr[i],weven);
    for (i=0;i<rOE;i++)
      prof2->Fill(j,oddevenptr[i]);

  }
  prof2->Write();
  prof->Write();
  delete prof;
  delete prof2;
  
}

/*-----------------------------------------------------------------------*/
/* histo to make difference between left and right strip */
void LRDiff(IMAGE2D* image,int x1in, int x1out, int x2in, int x2out, float sigma){

  /* x1 is 'left' and x2 'right' */
  /* diff is R - L */
  int i,j,remain;
  double *over,*overptr;
  int nover1 = x1out-x1in;
  int nover2 = x2out-x2in;
  float offset1, offset2 ;

  if (nover1 > nover2)
    over = new double[nover1];
  else 
    over = new double[nover2];
  
  /* 2nd : do the job */
  TH1F * histo = new TH1F("lrdiff","Diff between L and R Strip",100,0,0);
  TProfile * prof = new TProfile("lrprof","Diff between L and R Strip",image->ny,0,image->ny);

  for (j=0;j<image->ny;j++) {
    for (i=x1in;i<x1out;i++) 
      over[i-x1in] = RD_frame(image,i,j);

    /* purify the overscan from spurious */
    remain = nover1;
    overptr=over;
    if (sigma>0) 
      ut_trunc_sigma_unknown(&overptr,&remain,sigma);

    offset1 = gsl_stats_mean(overptr,1,remain);
    
    /* once again for strip 2*/
    for (i=x2in;i<x2out;i++) 
      over[i-x2in] = RD_frame(image,i,j);

    /* purify the overscan from spurious */
    overptr=over;
    remain=nover2;
    if (sigma>0) 
      ut_trunc_sigma_unknown(&overptr,&remain,sigma);

    offset2 = gsl_stats_mean(overptr,1,remain);
    
    histo->Fill(offset2-offset1);
    prof->Fill(j,offset2-offset1);
    
  }
 
  histo->Write();
  prof->Write();

  delete histo;
  delete prof;
  delete[]over;

}


/* ----------- error from overscan strip ---------------- */
void OverscanError(IMAGE2D* image,int xin, int xout,float sigma){
  
  int i,j,remain;
  double over[MAX_LINE],*overptr; 
  /* m : mean substracted ; a : all data*/
  int nover = xout-xin;
  float offset,var;

  TH1F * histo = new TH1F("overerr","Error from Oversvcan",100,0,0);
  TH1F * histom = new TH1F("overerrm","Error from Oversvcan mean substracted",100,0,0);
  

  /* then fill the histos*/
  for (j=0;j<image->ny;j++) {
    for (i=xin;i<xout;i++) 
      over[i-xin] = RD_frame(image,i,j);

    
    /* purify the overscan from spurious */
    remain=nover;
    overptr=over;
    if (sigma>0) 
      ut_trunc_sigma_unknown(&overptr,&remain,sigma);

    offset = gsl_stats_mean(overptr,1,remain);
    var=gsl_stats_variance_m(overptr,1,remain,offset);

    for (i=0;i<nover;i++) {
      float val = (overptr[i]-offset)*sqrt((remain+1.0)/remain);
      histom->Fill(val);
      histo->Fill(overptr[i]);
    }
  }
  
  histom->Write();
  histo->Write();
  
  delete histom;
  delete histo;
  
  
}



/*-------------- main ROOT routine ---------------------------- */
void rootify(IMAGE2D * image,int xin, int xout,int x2in, int x2out,float sigma) 
{
  /* note : x2 is for left overscan : should be less than xout*/

  /* makes the ROOT file */
  TString name(image->name);
  char match[] = ".fits";  

  if (name.Contains(match))
    name.ReplaceAll(".fits",".root");
  else
    name.Append(".root");
  //name.Prepend(image->path);
  printf("creating %s\n",name.Data());
  

  TFile * file = new TFile(name.Data(),"RECREATE");
  TH1F * histo = new TH1F("data","Histo of data values",1000,0,0);
  TH1F * histoo = new TH1F("ovsc","Histo of overscan values",1000,0,0);

  histo->SetBuffer(10000);
  

  /* table of values */
  TMatrix * mat = new TMatrix(image->nx,image->ny);
  
  for (int j=0;j<image->ny;j++) {
    for (int i=0;i<image->nx;i++)
      (*mat)(i,j) = RD_frame(image,i,j);
    for (int i=xin;i<xout;i++)
      histo->Fill(RD_frame(image,i,j));
    for (int i=x2in;i<x2out;i++)
      histoo->Fill(RD_frame(image,i,j));
  }
  mat->Write("DataMat");
  histo->Write();
  histoo->Write();
  
  delete mat;
  delete histo;
  delete histoo;
  

  /* others */
  VerticalProfile(image,xin,xout,sigma);
  OverscanProfile(image,x2in,x2out,sigma);
  OddEvenVerticalProfile(image,xin,xout,sigma);
  OverscanError(image,xin,xout,sigma);
  HorizontalProfile(image,0,image->ny,sigma);
  /* note the inversion */
  //LRDiff(image,x2in,x2out,xin,xout,sigma);
  
  
  /* close output*/
  file->Close();
  delete file;
}




/*------------ main ----------------------*/
int main(int argc, char **argv) {

  char **argval, **arglabel;
  IMAGE2D in;
  int xin,xout,x2in,x2out,cat;
  float sigma;
  char nom_cat[lg_name+1],text[81], im_name[lg_name+1];

  set_arglist("-in null -xin 6 -xout 1026 -x2in 1028 -x2out 1056 -sigma 5.0");
  init_session(argv,argc,&arglabel,&argval);

  /* xin and xout follows the convention :
     first number is 1 
     xin and xout are inside selection
  */

  get_argval(1,"%d", &xin);
  get_argval(2,"%d", &xout);
  get_argval(3,"%d", &x2in);
  get_argval(4,"%d", &x2out);
  get_argval(5,"%f", &sigma);

  /* check catalog */
  strcpy(nom_cat, argval[0]);
  cat = (strstr(nom_cat,".cat") != NULL);
  if (cat) {
    sprintf(text,"Image Catalog: %s", nom_cat);
    print_msg(text);
    RD_catalog(nom_cat, im_name);
  } else {
    strcpy(im_name, nom_cat);
  }

  /* loop on catalog */
  do {
    printf("Opening now %s\n",im_name);
    open_frame(&in,im_name,"I");
    /* parameters in rootify are for i=xin, i<xout */
    rootify(&in,xin-1,xout,x2in,x2out,sigma); 
    close_frame(&in);
    

  } while (cat && RD_catalog(nom_cat, im_name));


  exit_session(0);
  return(0);
}



