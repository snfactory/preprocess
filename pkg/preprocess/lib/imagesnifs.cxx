/* === Doxygen Comment ======================================= */
/*! 
 * \file          image.cxx
 * \copyright     (c) 2003 CRAL-Observatoire de Lyon
 * \date          Wed Aug  6 17:44:35 2003
 * \author        Emmanuel Gangler <e.gangler@ipnl.in2p3.fr>
 * \version       0.0
 * \brief         
 *                
 * $Id$
 **/
/* =========================================================== */

#include "utils.h"
#include "imagesnifs.hxx"
#include "section.hxx"


/* ##### IMAGE SNIFS ################################################## */

/* ===== constructor/Destructor ======================================= */
ImageSnifs::ImageSnifs (IoMethod_t Method, int MParam) : ImageSimple(Method, MParam){
  SetParanoMode(true);
}

/* ----- ImageSnifs copy ------------------------------ */
ImageSnifs::ImageSnifs(const ImageSnifs &image,char* newname,short newtype,int copydata,IoMethod_t Method, int MParam) 
  : ImageSimple(image, ut_create_check_name(newname), newtype, copydata,Method,MParam) 
{
  
    
  SetParanoMode(true);

  // if variance, try to copy if a meaningful naem can be done
  if (image.Variance()) {
    char varname[lg_name+1];
    ut_varname_from_imname(newname,varname);
    if (varname[0]) {
      ImageSimple* var = new ImageSimple(*image.Variance(),varname,newtype,copydata,Method,MParam);
      SetVarianceFrame(var);
    }
  }
}

/* ----- ImageSnifs open ------------------------------ */
ImageSnifs::ImageSnifs(char* name, char* mode,IoMethod_t Method, int MParam) 
  : ImageSimple( ut_open_check_name(name), mode,Method, MParam)
{

  SetParanoMode(true);

  // load varianece if available
  char varname[lg_name+1] ;
  ut_varname_from_imname(Name(),varname);
  if (varname[0] && exist(varname)) {
    ImageSimple* var = new ImageSimple(varname,mode,Method,MParam);
    SetVarianceFrame(var);
  }
}


/* ----- ~ImageSnifs ------------------------------ */
ImageSnifs::~ImageSnifs(){  
}

/* ===== Methods ======================================= */

/* ----- BuildSubImage  ------------------------------ */
ImageSnifs* ImageSnifs::BuildSubImage(Section* Sec,char *Name){
  // Build a raster image

  char varname[lg_name+1], secString[lg_name+1] ;
  ut_varname_from_imname(Name,varname);
  
  // raw contructor (create_frame has to be thought about)
  ImageSnifs* sub = new ImageSnifs();
  if (Variance() && varname[0]!='\0') {
    ImageSimple* variance = new ImageSimple();
    sub->SetVarianceFrame(variance);
  }

  // create the frames
  sub->CreateFrame(Name, Sec->XLength() ,Sec->YLength());
  sub->ImportHeader(this);
  if (sub->Variance()) {
    sub->Variance()->CreateFrame(varname, Sec->XLength() ,Sec->YLength());
    sub->Variance()->ImportHeader(Variance());
  }
  sub->ImportSection(this,Sec,1,1);

  // update the keys : RASTER
  Sec->GetString(secString);
  sub->WrDesc("RASTER",CHAR,lg_name+1,secString);
  if (sub->Variance())
    sub->Variance()->WrDesc("RASTER",CHAR,lg_name+1,secString);
  // DATASEC
  sprintf(secString,"[%d:%d,%d:%d]",1,Sec->XLength() ,1,Sec->YLength());
  sub->WrDesc("DATASEC",CHAR,lg_name+1,secString);
  if (sub->Variance())
    sub->Variance()->WrDesc("DATASEC",CHAR,lg_name+1,secString);

  return sub;
}



/* ----- SubstractOverscan ------------------------------ */
void ImageSnifs::SubstractOverscan() {

  int overscanDone;

  if (ParanoMode()) {
    // check fcalsses !

    // check overscan was not already substracted
    if ( RdIfDesc("OVSCDONE",INT, 1, &overscanDone) >= 0 
         && overscanDone ){
      print_error("Overscan already substracted for %s\n",Name());
      print_error("Nothing Done%s\n",Name());
      return;
    }
  }
  
  // Get the overscan bounds from file
  char bias[lg_name+1];
  RdDesc("BIASSEC",CHAR,lg_name+1,bias);
  Section * Sec = new Section(bias);
  Image()->SubstractOverscan(Sec);  
  delete Sec;

  overscanDone=1;
  WrDesc("OVSCDONE",INT,1,&overscanDone);

  // set the biasframe if allowed
  char obstype[lg_name+1];
  RdDesc("OBSTYPE",CHAR,lg_name+1,obstype);
  if (! strcmp(obstype,"BIAS")) {
    int biasFrame=1;
    WrDesc("BIASFRAM",INT,1,&biasFrame);
  }

}

/* ----- odd-Even ----------------------------------------------------- */
void ImageSnifs::OddEvenCorrect() {
  // odd-even substraction
  //

  double param[2];
  
  if (ParanoMode()) {
    // check fcalsses !
    
    // check type of CCD (is operation allowed ?)

    // check overscan was not already substracted
    // returns 0 if not done !
    if ( RdIfDesc("OEPARAM",DOUBLE, 2, param) > 0 ){
      print_error("Odd-Even already substracted for %s\n",Name());
      print_error("Nothing Done%s\n",Name());
      return;
    }
  }

  // Get the overscan bounds from file
  char bias[lg_name+1];
  RdDesc("BIASSEC",CHAR,lg_name+1,bias);
  Section* Sec=new Section(bias);

  // the magic 2.0 parameter is because of CTE on CCD EEV#1 (blue)
  Image()->OddEvenCorrect(Sec,param,2.0);
  delete Sec;

  WrDesc("OEPARAM",DOUBLE, 2, param);
}

/* -----  SubstractBias ------------------------------------------------- */
void ImageSnifs::SubstractBias(ImageSnifs* Bias) {
  // substracts a bias image
  int biasDone,overscanDone,biasFrame;

  if (ParanoMode()) {
    // check fclasses
    Bias->RdDesc("BIASFRAM",INT,1,&biasFrame);
    if (!biasFrame) {
      print_error("ImageSnifs:SubstractBias %s is not a bias frame",Bias->Name());
      return;
    }

    // check bias was not already substracted
    if ( RdIfDesc("BIASDONE",INT, 1, &biasDone) >= 0 
         && biasDone ){
      print_error("ImageSnifs::SubstractBias already substracted for %s",Name());
      return;
    }

    // check overscan was already substracted
    RdDesc("OVSCDONE",INT, 1, &overscanDone);
    if (!overscanDone ){
      print_error("Overscan was not substracted yet for %s\n",Name());
      return;
    }

    // check images are of the same size
    if (Nx()!=Bias->Nx() || Ny()!=Bias->Ny()) {
      print_error("ImageSnifs::SubstractBias : images are not of same size");
      return;
    }    
  }
  
  // do the real job
  Image()->Add(Bias->Image(),-1);

  biasDone=1;
  WrDesc("BIASDONE",INT,1,&biasDone);

  // set the darkframe if allowed
  char obstype[lg_name+1];
  RdDesc("OBSTYPE",CHAR,lg_name+1,obstype);
  if (! strcmp(obstype,"DARK")) {
    int darkFrame=1;
    WrDesc("DARKFRAM",INT,1,&darkFrame);
  }

}

/* -----  SubstractDark ------------------------------------------------- */
void ImageSnifs::SubstractDark(ImageSnifs* Dark) {
  // substracts a dark image
  int darkDone,overscanDone,darkFrame;
  double exptime,darktime;

  if (ParanoMode()) {
    // check fclasses
    Dark->RdDesc("DARKFRAM",INT,1,&darkFrame);
    if (!darkFrame) {
      print_error("ImageSnifs::SubstractDark %s is not a dark frame",Dark->Name());
      return;
    }

    // check dark was not already substracted
    if ( RdIfDesc("DARKDONE",INT, 1, &darkDone) >= 0 
         && darkDone ){
      print_error("ImageSnifs::SubstractDark already substracted for %s",Name());
      return;
    }

    // check overscan was already substracted
    RdDesc("OVSCDONE",INT, 1, &overscanDone);
    if (!overscanDone ){
      print_error("ImageSnifs::SubstractDark Overscan was not substracted yet for %s\n",Name());
      return;
    }

    // check images are of the same size
    if (Nx()!=Dark->Nx() || Ny()!=Dark->Ny()) {
      print_error("ImageSnifs::SubstractDark : images are not of same size");
      return;
    }
  }

  // do the real job
  RdDesc("EXPTIME",DOUBLE,1,&exptime);
  Dark->RdDesc("EXPTIME",DOUBLE,1,&darktime);
  Image()->Add(Dark->Image(),-exptime/darktime);

  darkDone=1;
  WrDesc("DARKDONE",INT,1,&darkDone);


}


/* -----  ApplyFlat     ------------------------------------------------- */
void ImageSnifs::ApplyFlat(ImageSnifs* Flat) {
  // apply a flat frame
  int flatDone,biasDone,flatFrame;

  if (ParanoMode()) {
    // check fclasses
    RdDesc("BIASDONE",INT, 1, &biasDone);
    if ( !biasDone ){
      print_error("ImageSnifs::ApplyFlat : Bias was not substracted yet for %s\n",Name());
      return;
    }

    Flat->RdDesc("FLATFRAM",INT, 1, &flatFrame);
    if ( !flatFrame ){
      print_error("ImageSnifs::ApplyFlat : Bias was not substracted yet for %s\n",Name());
      return;
    }

    // check bias (so also overscan) was already substracted
    RdDesc("BIASDONE",INT, 1, &biasDone);
    if ( !biasDone ){
      print_error("ImageSnifs::ApplyFlat : Bias was not substracted yet for %s\n",Name());
      return;
    }

    // check dark was not already substracted
    if ( RdIfDesc("FLATDONE",INT, 1, &flatDone) >= 0 
         && flatDone ){
      print_warning("ImageSnifs::ApplyFlat already done for %s",Name());
    }

    // check images are of the same size
    if (Nx()!=Flat->Nx() || Ny()!=Flat->Ny()) {
      print_error("ImageSnifs::ApplyFlat : images are not of same size");
      return;
    }
  }

  // do the real job
  Image()->Divide(Flat->Image());

  flatDone=1;
  WrDesc("FLATDONE",INT,1,&flatDone);
  int eleCalib=0;
  WrDesc("ELECALIB",INT,1,&eleCalib);
}

/* -----  BuildFlat     ------------------------------------------------- */

void ImageSnifs::BuildFlat() {
  // transforms the image into a flat frame
  int biasDone,flatFrame;

  if (ParanoMode()) {
    // check fclasses
    char obstype[lg_name+1];
    RdDesc("OBSTYPE",CHAR,lg_name+1,obstype);
    if ( strcmp(obstype,"OBJECT")) {
       print_error("ImageSnifs::BuildFlat : %s is not a flat frame",Name());
      return;
    }
    // check it is not already a flat frame
    if ( RdIfDesc("FLATFRAM",INT, 1, &flatFrame) >= 0 
         && flatFrame ){
      print_warning("ImageSnifs::buildFlat %s is already a flat frame %s",Name());
    }

    // check bias was substracted
    RdDesc("BIASDONE",INT, 1, &biasDone);
    if ( !biasDone ){
      print_error("ImageSnifs::BuildFlat : Bias was not substracted yet for %s",Name());
      return;
    }
  }

  // Get the data bounds from file
  char data[lg_name+1];
  RdDesc("DATASEC",CHAR,lg_name+1,data);
  Section* Sec=new Section(data);
  // do the real job
  double norm = Image()->MeanValue(Sec,100);
  Image()->Scale(1/norm);
  delete Sec;
  
  flatFrame=1;
  WrDesc("FLATFRAM",INT,1,&flatFrame);
  int eleCalib=0;
  WrDesc("ELECALIB",INT,1,&eleCalib);
}

/* -----  AddPoissonNoise ----------------------------------------------- */

void ImageSnifs::AddPoissonNoise() {
  // Adds the poisson noise to the image.
  //
  // Side effect : image is bounded to reject negative values.

  int poisNoise=0;

  // ParanoChecks
  if (ParanoMode()) {
    // check if in electrons
    // and if noise was not already added
    int eleCalib=0;
    if ( RdIfDesc("ELECALIB",INT, 1, &eleCalib) <= 0 
         || !eleCalib) {
      print_error("ImageSnifs::AddPoissonNoise %s is not in electrons",Name());
      return;
    }
    if ( Variance()->RdIfDesc("POISNOIS",INT, 1, &poisNoise) > 0 
         && poisNoise) {
      print_error("ImageSnifs::AddPoissonNoise already performed for %s",Name());
      return;
    }
  }
  if (!Variance()) {
    print_error("ImageSnifs::AddPoissonNoise but %s has no variance",Name());
    return;  
  }
  // End of checks

  // Go
  Image()->AddPoissonNoise();

  // Set variables
  poisNoise=1;
  Variance()->WrDesc("POISNOIS",INT, 1, &poisNoise);
}

/* -----  AddPoissonNoise ----------------------------------------------- */
void ImageSnifs::HandleSaturation() {
  // puts variance to infinity for saturated pixels
  int saturate;
  RdDesc("SATURATE",INT,1,&saturate);
  Image()->HandleSaturation(saturate-0.5);
}

/* ===== Hacks ======================================= */

/* ----- HackFitsKeywords ---------------------------------------------- */
void ImageSnifs::HackFitsKeywords() {

  // The following is not valid for any raster yet !
  char key[lg_name+1];
  RdDesc("RASTER",CHAR,lg_name+1,key);
  if (strcmp(key,"FULL")){
    print_error("ImageSnifs::HackFitsKeywords : only works for FULL raster");
    return;
  }

  // The BiasSec is not correct
  int x1,x2,y1,y2;
  RdDesc("BIASSEC",CHAR,lg_name+1,key);
  sscanf(key,"[%d:%d,%d:%d]",&x1,&x2,&y1,&y2);
  sprintf(key,"[%d:%d,%d:%d]",x1+2,x2,1,Ny());
  WrDesc("BIASSEC",CHAR,lg_name+1,key);
  
  // the DATASEC needs help too
  RdDesc("DATASEC",CHAR,lg_name+1,key);
  sscanf(key,"[%d:%d,%d:%d]",&x1,&x2,&y1,&y2);
  sprintf(key,"[%d:%d,%d:%d]",4,x2+2,1,y2+5);
  WrDesc("DATASEC",CHAR,lg_name+1,key);

}


/* ===== Processing Managment ======================================= */

/* ----- CreateVarianceFrame ------------------------------ */
void ImageSnifs::CreateVarianceFrame(char* name) {

  if (ParanoMode()) {
    if (Variance()) {
      print_error(" ImageSnifs::CreateVarianceFrame : already a variance in %s",Name());
      return;
    }
  
    if (!name[0]) {
      ut_varname_from_imname(Name(),name);
      if (!name[0]) {
         print_error(" ImageSnifs::CreateVarianceFrame : no name provided",Name());
      return;
      }
    }

  } // ParanoMode

  // Builds the variance frame
  ImageSimple * var = new ImageSimple(*Image() , name,FLOAT,0,IoType(),Ny());
  var->SetTo(0);
  SetVarianceFrame(var);
  HandleSaturation();
  
}


/* ----- AddOverscanVariance ------------------------------ */
void ImageSnifs::AddOverscanVariance() {

  int ovscNoise;
  if (ParanoMode()) {
    if (Variance()->RdIfDesc("OVSCNOIS",INT, 1, &ovscNoise) > 0 
     && ovscNoise ) {
       print_error(" ImageSnifs::AddOverscanVariance already done in %s",Name());
      return;
      }
  }

  // Fills with analysis of RMS strip
  char bias[lg_name+1];
  RdDesc("BIASSEC",CHAR,lg_name+1,bias);
  Section * Sec = new Section(bias);
  // we allow 1 outlier...
  //double rms = Image()->OverscanRms(Sec,ut_fraction_sigcut(1.0/(Sec->XLength()-0.5)));
  // no outlier
  double rms = Image()->OverscanRms(Sec,0);
  WrDesc("RDNOISE",DOUBLE,1,&rms);
  Variance()->Add(rms*rms);
  delete Sec;

}

/* ----- Assembled2Dark ------------------------------ */
void ImageSnifs::Assembled2Dark() {
  AddPoissonNoise();
}

/* ----- Assembled2Flat ------------------------------ */
void ImageSnifs::Assembled2Flat(ImageSnifs * Dark) {
  Assembled2Dark();
  if (Dark) 
    SubstractDark( Dark );
  BuildFlat();
}

/* ----- Assembled2Preprocessed ------------------------------ */
void ImageSnifs::Assembled2Preprocessed(ImageSnifs * Dark,ImageSnifs* Flat) {
  Assembled2Dark();
  if (Dark) 
    SubstractDark(Dark);
  if (Flat)
    ApplyFlat(Flat);
}

/* ----- CompareIntDesc ------------------------------ */
int ImageSnifs::IdenticalPreprocDesc(ImageSnifs * ToCheck, char* Desc, int* val1, int* val2) {
  // a very specialized routine to help check if files are alike
  if ( RdIfDesc(Desc,INT, 1, val1) <0 )
    *val1 =0;
  if ( ToCheck->RdIfDesc(Desc,INT, 1, val2) <0 )
    *val2 =0;
  if (*val1 != *val2) {
    print_msg("ImageSnifs::IdenticalPreprocDesc %s, %d is not %d",Desc,*val1, *val2);
    return 0;
  }
  else 
    return 1;
}

/* ----- CanBeStackedWith ------------------------------ */
int ImageSnifs::CanBeStackedWith(ImageSnifs * ToCheck) {
  if (ToCheck->Nx() != Nx() || ToCheck->Ny() != Ny()) {
    print_msg("ImageSnifs::CanBeStackedWith Size differs");
    return 0;
  }
  int val1,val2;
  if (!IdenticalPreprocDesc(ToCheck,"OVSCDONE",&val1,&val2))
    return 0;
  if (!IdenticalPreprocDesc(ToCheck,"BIASFRAM",&val1,&val2))
    return 0;
  if (!IdenticalPreprocDesc(ToCheck,"BIASDONE",&val1,&val2))
    return 0;
  if (!IdenticalPreprocDesc(ToCheck,"DARKFRAM",&val1,&val2))
    return 0;
  if (!IdenticalPreprocDesc(ToCheck,"DARKDONE",&val1,&val2))
    return 0;
  if (!IdenticalPreprocDesc(ToCheck,"FLATFRAM",&val1,&val2))
    return 0;
   if (!IdenticalPreprocDesc(ToCheck,"FLATDONE",&val1,&val2))
    return 0;
   if (!IdenticalPreprocDesc(ToCheck,"ELECALIB",&val1,&val2))
    return 0;
   if ((Variance() && !ToCheck->Variance()) || (!Variance() && ToCheck->Variance())) {
     print_msg("ImageSnifs::CanBeStackedWith Variance frame differs");
     return 0;
   }
   if (Variance() ) {
     if ( RdIfDesc("POISNOIS",INT, 1, &val1) <0 )
       val1 =0;
     if ( ToCheck->RdIfDesc("POISNOIS",INT, 1, &val2))
       val2 =0;
     if (val1 != val2) {
       print_msg("ImageSnifs::CanBeStackedWith POISNOIS, %d is not %d",val1, val2);
       return 0;
     }
   }
   return 1;
}



/*====== Specific photometry tools ================== */

/* ----- BuildFilterSection ------------------------------ */
void ImageSnifs::FilterSectionFill(Section * Sec,int Filter) {

  // Sec will be overwritten in the process

  switch (Filter) {
  case 0 : 
    Sec->SetString("[1:2048,1:686]");
    break;
  case 1 :
    Sec->SetString("[1:772,1192:2393]");
    break;
  case 2:
    Sec->SetString("[1277:2048,1192:2393]");
    break;
  case 3:
    Sec->SetString("[1:772,2899:4096]");
    break;
  case 4:
    Sec->SetString("[1:772,2899:4096]");
    break;
  default:
    print_error("ImageSnifs::BuildFilterSection : paramter out of 0-4 range");
  }
}

/* ----- SplitMultifilter ------------------------------ */
void ImageSnifs::SplitMultiFilter(char* NameRecipee,ImageSnifs** Subs){
  // the NameRecipee shall contain one %s which will
  // be transformated with 'a' ... 'e'
  // The Subs should contain 5 allocated ImageSnifs*

  // as usual, paranochecks
  if (ParanoMode()) {
  }
  // and checks
  if (!strstr(NameRecipee,"%c")){
    print_error("ImageSnifs::SplitMultifilter Name shall contain %c");
    return;
  }

  // Do it :
  char imageName[lg_name+1];
  Section sec;

  for (int i=0;i<5;i++) {
    // names for the image
    sprintf(imageName,NameRecipee,'a'+i);
    FilterSectionFill(&sec,i);
    Subs[i] = BuildSubImage(&sec,imageName);
  }
}


