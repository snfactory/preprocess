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

#include "section.hxx"


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

/* ----- MeanValue ----------------------------------------------------- */
double MeanValue(IMAGE2D* image,Section* Sec,int step=1){
  //Computes the mean value of the section data
  // the step factor is here to speed-up the computation
  
  double sum;
  int i,j,n=0;

  i = Sec->XFirst();
  for ( j=Sec->YFirst() ; j<Sec->YLast() ; ) {
    for ( ; i<Sec->XLast() ; i+=step ) {
      sum += RD_frame(image,i,j);
      n++;
    }
    while(i>=Sec->XLast()) {
      i-=Sec->XLength();
      j++;
    }
  }
  return sum/n;
}

/*-------------- élairement moyen ----------------------------- */
double Meanlight(IMAGE2D* image, int xin, int xout)
{
  Section sec(xin,xout,1,image->ny);
  return MeanValue(image,&sec,1);
}

/*-------------- hf ------------------------------------------- */
void HF(IMAGE2D* image, int xin, int xout)
{
  double rapporth,rapportv,minh,minv,mintot,maxh,maxv,maxtot;

  minh = 1;
  minv = 1;
  maxh = 0;
  maxv = 0;

  for (int j=0; j<(image->ny)-1;j++){
    for (int i=xin; i<xout-1;i++){
      rapporth = (RD_frame(image,i,j))/(RD_frame(image,i,j+1));
      rapportv = (RD_frame(image,i,j))/(RD_frame(image,i+1,j));
      if (rapporth<minh) 
        minh = rapporth;
      else if (rapporth>maxh) 
        maxh = rapporth; 
      if (rapportv<minv) 
        minv = rapportv;
      else if (rapportv>maxv) 
        maxv = rapportv;
    }
  }
  
  if (minh<minv){
    mintot = minh;
  }
  else mintot = minv;
  
  if (maxh<maxv){
    maxtot = maxv;
  }
  else maxtot = maxh;
  
  
  TH1F* histoh = new TH1F("rapporth","Rapport horizontal",1000,minh-1,maxh+1);
  TH1F* histov = new TH1F("rapportv","Rapport vertical",1000,minv-1,maxv+1);
  TH1F* histotot = new TH1F("rapportot","Rapport aux voisins",1000,mintot-1,maxtot+1);
    

  for (int j=0; j<(image->ny)-1;j++){
    for (int i=xin; i<xout-1;i++){
      rapporth = (RD_frame(image,i,j))/(RD_frame(image,i,j+1));
      rapportv = (RD_frame(image,i,j))/(RD_frame(image,i+1,j));
      histoh->Fill(rapporth);
      histov->Fill(rapportv);
      histotot->Fill(rapporth);
      histotot->Fill(rapportv);
    }
  }
  
  histoh->Write();
  histov->Write();
  histotot->Write();
  
  delete histoh;
  delete histov;
  delete histotot;
  
}


/*-------------- main ROOT routine ---------------------------- */
void rootify(IMAGE2D * image, char* out_name, Section sec, int x2in, int x2out,float sigma) 
{
  /* note : x2 is for left overscan : should be less than xout*/


  double min, max;
  min = RD_frame(image,0,0);
  max = RD_frame(image,0,0);
  
  for (int j=sec.YFirst();j<sec.YLast();j++){
    for (int i=sec.XFirst();i<sec.XLast();i++){
      if (RD_frame(image,i,j)<min) 
        min = RD_frame(image,i,j);
      else if (RD_frame(image,i,j)>max) 
        max = RD_frame(image,i,j);
    }
  }
  

  printf("creating %s\n",out_name);
  TFile * file = new TFile(out_name,"RECREATE");
  TH1F * histo = new TH1F("data","Histo of data values",1000,min-1,max+1);
  TH1F * histoo = new TH1F("ovsc","Histo of overscan values",1000,0,0);

  histo->SetBuffer(10000);
  

  /* table of values */
  TMatrix * mat = new TMatrix(sec.XLast(),sec.YLast());
  
  for (int j=sec.YFirst();j<sec.YLast();j++) {
    for (int i=sec.XFirst();i<sec.XLast();i++){
      (*mat)(i,j) = RD_frame(image,i,j);
      histo->Fill(RD_frame(image,i,j));
    }
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
  VerticalProfile(image,sec.XFirst(),sec.XLast(),sigma);
  OverscanProfile(image,x2in,x2out,sigma);
  OddEvenVerticalProfile(image,sec.XFirst(),sec.XLast(),sigma);
  OverscanError(image,sec.XFirst(),sec.XLast(),sigma);
  HorizontalProfile(image,sec.YFirst(),sec.YLast(),sigma);
  HF(image,sec.XFirst(),sec.XLast());
  printf("éclairement moyen = %f\n", Meanlight(image, sec.XFirst(), sec.XLast()));
  /* note the inversion */
  //LRDiff(image,x2in,x2out,sec.XFirst(),sec.XLast(),sigma);
  
  
  /* close output*/
  file->Close();
  delete file;
}




/*------------ main ----------------------*/
int main(int argc, char **argv) {

  char **argval, **arglabel;
  IMAGE2D in;
  TFile out;
  int x2in,x2out,cat;
  float sigma;
  char nom_cat[lg_name+1],text[81], im_name[lg_name+1], nom_cat_out[lg_name+1], out_name[lg_name+1];

  set_arglist("-in null -out null -sec [6:1026,1:4102] -x2in 1028 -x2out 1056 -sigma 5.0");
  init_session(argv,argc,&arglabel,&argval);

  /* xin and xout follows the convention :
     first number is 1 
     xin and xout are inside selection
  */

  Section sec(argval[2]);
  
  get_argval(3,"%d", &x2in);
  get_argval(4,"%d", &x2out);
  get_argval(5,"%f", &sigma);

  /* check catalog */
  strcpy(nom_cat, argval[0]);
  strcpy(nom_cat_out, argval[1]);
  cat = (strstr(nom_cat,".cat") != NULL);
  if (cat) {
    if (!strstr(nom_cat_out,".cat")) {
      print_error ("-out shall be a catalog");
      exit_session(1);
    }
    
    sprintf(text,"Image Catalog: %s", nom_cat);
    print_msg(text);
    RD_catalog(nom_cat, im_name);
    RD_catalog(nom_cat_out, out_name);
  } else {
    strcpy(im_name, nom_cat);
    strcpy(out_name, nom_cat_out);
  }

  /* loop on catalog */
  do {
    printf("Opening now %s\n",im_name);
    open_frame(&in,im_name,"I");
    /* parameters in rootify are for i=xin, i<xout */
    rootify(&in,out_name,sec,x2in,x2out,sigma);  
  close_frame(&in);
    

  } while (cat && RD_catalog(nom_cat, im_name) &&  RD_catalog(nom_cat_out, out_name));


  exit_session(0);
  return(0);
}



