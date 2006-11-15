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

/* ----- global includes ----- */
#include "fclass_snifs.h"
#include "snifs_defs.h"

/* ----- local includes ----- */
#include "utils.h"
#include "imagesnifs.hxx"
#include "section.hxx"
#include "algocams.hxx"
#include "snifs_const.h"
#include "analyser.hxx"

//static const char kChannelName[4][lg_name+1]={"Blue channel","Red channel","Photometry","Guiding"} ;


/* ##### IMAGE SNIFS ################################################## */

/* ===== constructor/Destructor ======================================= */
ImageSnifs::ImageSnifs (IoMethod_t Method, int MParam) : ImageSimple(Method, MParam){
  SetParanoMode(true);
  fAlgo=0;
}

/* ----- ImageSnifs copy ------------------------------ */
ImageSnifs::ImageSnifs(const ImageSnifs &image,char* newname,short newtype,int copydata,IoMethod_t Method, int MParam) 
  : ImageSimple(image, ut_create_check_name(newname), newtype, copydata,Method,MParam) 
{
  
    
  SetParanoMode(image.ParanoMode());

  // if variance, try to copy if a meaningful naem can be done
  if (image.Variance()) {
    char varname[lg_name+1];
    ut_varname_from_imname(newname,varname);
    if (varname[0]) {
      ImageSimple* var = new ImageSimple(*image.Variance(),varname,newtype,copydata,Method,MParam);
      SetVarianceFrame(var);
    }
  }
  fAlgo=0;
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
  fAlgo=0;
}


/* ----- ~ImageSnifs ------------------------------ */
ImageSnifs::~ImageSnifs(){  
  if (fAlgo)
    delete fAlgo;
}

/* ===== Utilities ======================================= */

/* ----- FClass  ------------------------------ */
int ImageSnifs::GetFClass(){
  int fclass=DONT_KNOW;
  RdIfDesc("FCLASS",INT,1,&fclass);
  return fclass;
}

void ImageSnifs::SetFClass(int Fclass){
  WrDesc("FCLASS",INT,1,&Fclass);
}

/* ----- Algo  ------------------------------ */
void ImageSnifs::SetAlgo(char* soft){
  WrDesc("IMAGESW",CHAR,lg_name+1,soft);
  if (fAlgo!=0)
    delete fAlgo;
  fAlgo=0;
}

AlgoCams* ImageSnifs::Algo(){
  char soft[lg_name+1];
  if (fAlgo)
    return fAlgo;
  if (RdIfDesc("IMAGESW",CHAR,lg_name+1,soft) <=0 )
    return 0;
  if (!strcmp(soft,"DETCOM"))
    fAlgo = new AlgoDetcom;
  if (!strcmp(soft,"SNFDETCOM"))
    fAlgo = new AlgoSnfDetcom;
  if (!strcmp(soft,"OTCOM"))
    fAlgo = new AlgoOtcom;
  return fAlgo;
}

/* ----- Channel  ------------------------------ */
void ImageSnifs::SetChannel(int channel){
  char name[lg_name+1];
  name[0]= '\0';
  
  for (int i=0;i<kNChannel;i++) {
    if (channel & (1<<i)) {
      if ( name[0] )
        sprintf(name,"%s + %s",name,kChannelName[i]);
      else
        strcpy(name,kChannelName[i]);
    }
  }
  WrDesc("CHANNEL",CHAR,lg_name+1,name);
}

/* ----- Channel  ------------------------------ */
int ImageSnifs::GetChannel(){
  char name[lg_name+1];
  int channel = kUnknown;

  if (RdIfDesc("CHANNEL",CHAR,lg_name+1,name) <= 0)
    return kUnknown;
  
  for (int i=0;i<kNChannel;i++)
    if (strchr(name,kChannelName[i][0]))
      channel += (1<<i);

  return channel;
}

/* ----- DataSec  ------------------------------ */
Section* ImageSnifs::DataSec(){
  char data[lg_name+1];
  RdDesc("DATASEC",CHAR,lg_name+1,data);
  Section* sec=new Section(data);
  return sec;
}

/* ----- BiasSec  ------------------------------ */
Section* ImageSnifs::BiasSec(){
  char data[lg_name+1];
  RdDesc("BIASSEC",CHAR,lg_name+1,data);
  Section* sec=new Section(data);
  return sec;
}


/* ----- HasOverscan  ------------------------------ */
int ImageSnifs::HasOverscan(){
  // return 1 if there is a valid overscan region in the fits keywords

  char biasSecString[lg_name+1];
  Image()->RdDesc("BIASSEC",CHAR,lg_name+1,biasSecString);
  int x1,x2,y1,y2,nmatch;
  nmatch = sscanf(biasSecString,"[%d:%d,%d:%d]",&x1,&x2,&y1,&y2);
  if (nmatch<4) {
    return 0;
  }
  return 1;
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
#ifdef OLD
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
  Section * Sec = Algo()->SafeOverscanStrip(this);
  Image()->SubstractOverscan(Sec);  
  delete Sec;

  overscanDone=1;
  WrDesc("OVSCDONE",INT,1,&overscanDone);

  // set the biasframe if allowed
  char obstype[lg_name+1];
  RdDesc("OBSTYPE",CHAR,lg_name+1,obstype);
  if (! strcmp(obstype,"BIAS") ||GetFClass() == RAW_BIAS ) {
    int biasFrame=1;
    WrDesc("BIASFRAM",INT,1,&biasFrame);
  }
  if (GetFClass() == RAW_BIAS)
    SetFClass(BIAS_FRAME);
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

  // Get the overscan bounds
  Section* Sec=Algo()->SafeOverscanStrip(this);

  // the magic 2.0 parameter is because of CTE on CCD EEV#1 (blue)
  Image()->OddEvenCorrect(Sec,param,2.0);
  delete Sec;

  WrDesc("OEPARAM",DOUBLE, 2, param);
}
#endif

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
  if (! strcmp(obstype,"DARK") || GetFClass()==RAW_DARK_FRAME ) {
    int darkFrame=1;
    WrDesc("DARKFRAM",INT,1,&darkFrame);
  }
}


/* -----  SubstractBias ------------------------------------------------- */
void ImageSnifs::UpdateFClass() {
  if (GetFClass()==RAW_DARK_FRAME)
    SetFClass(PRE_DARK_FRAME);
  if (GetFClass()==RAW_CAL_FRAME)
    SetFClass(PRE_CAL_FRAME);
  if (GetFClass()==RAW_CON_FRAME)
    SetFClass(PRE_CON_FRAME);
  if (GetFClass()==RAW_SKY_FRAME)
    SetFClass(PRE_SKY_FRAME);
  if (GetFClass()==RAW_OBJ_FRAME)
    SetFClass(PRE_OBJ_FRAME);
  if (GetFClass()==RAW_DOM_FRAME)
    SetFClass(PRE_DOM_FRAME);
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
  Section* Sec= DataSec();
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

  // by-pass check
  if (!Variance()) {
    return;  
  }

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
  // End of checks

  // Go
  Image()->AddPoissonNoise();

  // Set variables
  poisNoise=1;
  Variance()->WrDesc("POISNOIS",INT, 1, &poisNoise);
}

/* -----  HandleSaturation ----------------------------------------------- */
void ImageSnifs::HandleSaturation() {
  // puts variance to infinity for saturated pixels
  int saturate, nsat;
  RdDesc("SATURATE",INT,1,&saturate);
  nsat = Image()->HandleSaturation(saturate-0.5);
  // not so good idea ...
  //  WrDesc("NSATU",INT,1,&nsat);
}

/* -----  HandleCosmetics ----------------------------------------------- */
void ImageSnifs::HandleCosmetics() {
  // puts variance to infinity for some definite pixels
  if (Image()->Variance()) {

    // Exclude non-standard binning
    int nbin[2];
    Image()->RdDesc("CCDBIN",INT,2,&nbin);
    if (nbin[0]!=1 || nbin[1]!=1) {
      print_msg("ImageSnifs::HandleCosmetics refuses non-standard binning");
      return ;
    }

    int nmax=0;
    const int * data;
    if (GetChannel()==kRedChannel) {
      if (Variance())
        nmax = kBadSectionsRed;
      else
        nmax=kBadSectionsRedNoVar;
      data=kBadSectionsDataRed;
    }
    if (GetChannel()&kPhotometric) {
      nmax = kBadSectionsPhot;
      data=kBadSectionsDataPhot;
    }

    // Non-standard raster
    char ccdSecString[lg_name+1];
    Image()->RdDesc("CCDSEC",CHAR,lg_name+1,ccdSecString);
    if (!strcmp(ccdSecString,"Undefined")) {
      print_warning("ImageSnifs::HandleCosmetics Undefined CCDSEC");
      return;
    }
    Section ccdSec(ccdSecString);

    for (int n=0;n<nmax;n++) {
      Section Sec(data[n*4],data[n*4+1],data[n*4+2],data[n*4+3]);
      for (int iy = Sec.YFirst(); iy<Sec.YLast();iy++)
        if (  iy-ccdSec.YFirst() < Ny() && 
              iy-ccdSec.YFirst() >=0)
          for (int ix = Sec.XFirst(); ix<Sec.XLast();ix++) {
            if (ix-ccdSec.XFirst() < Nx() && 
                ix-ccdSec.XFirst() >=0 )
              Image()->Variance()->WrFrame(ix-ccdSec.XFirst(),iy-ccdSec.YFirst(),ut_big_value);
        }
    }
  } // if variance
}

/* -----  HandleCosmetics ----------------------------------------------- */
void ImageSnifs::CheatCosmetics() {
  // puts some number where there is a cosmetic in order to have
  // a fast processing working
  // it builds a linear interpolation to the data

    // Exclude non-standard binning
    int nbin[2];
    Image()->RdDesc("CCDBIN",INT,2,&nbin);
    if (nbin[0]!=1 || nbin[1]!=1) {
      print_msg("ImageSnifs::CheatCosmetics refuses non-standard binning");
      return ;
    }

    int nmax=0;
    const int * data;
    if (GetChannel()==kRedChannel) {
      if (Variance())
        nmax = kBadSectionsRed;
      else
        nmax=kBadSectionsRedNoVar;
      data=kBadSectionsDataRed;
    }
    if (GetChannel()&kPhotometric) {
      nmax = kBadSectionsPhot;
      data=kBadSectionsDataPhot;
    }

    // Non-standard raster
    char ccdSecString[lg_name+1];
    Image()->RdDesc("CCDSEC",CHAR,lg_name+1,ccdSecString);
    if (!strcmp(ccdSecString,"Undefined")) {
      print_warning("ImageSnifs::CheatCosmetics Undefined CCDSEC");
      return;
    }
    Section ccdSec(ccdSecString);

    for (int n=0;n<nmax;n++) {
      Section Sec(data[n*4],data[n*4+1],data[n*4+2],data[n*4+3]);
      for (int iy = Sec.YFirst(); iy<Sec.YLast();iy++) {
        if (  iy-ccdSec.YFirst() < Ny() && 
              iy-ccdSec.YFirst() >=0) {
          
          double p1[2];
          if (Sec.XFirst()-ccdSec.XFirst()<=0) {
            if (Sec.XLast()-ccdSec.XFirst() < Nx() && Sec.XLast()-ccdSec.XFirst() >=0)
              p1[0]=Image()->RdFrame(Sec.XLast()-ccdSec.XFirst(),iy-ccdSec.YFirst());
            // nothing can be done
            else p1[0]=0;
            p1[1]=0;
          } else if (Sec.XLast()-ccdSec.XFirst()>=Nx()) {
            if (Sec.XFirst()-1-ccdSec.XFirst() < Nx() && 
                Sec.XFirst()-1-ccdSec.XFirst() >=0)
              p1[0]=Image()->RdFrame(Sec.XFirst()-1-ccdSec.XFirst(),iy-ccdSec.YFirst());
            else p1[0]=0;
            p1[1]=0;
          } else {
            p1[0] = Image()->RdFrame(Sec.XFirst()-1-ccdSec.XFirst(),iy-ccdSec.YFirst());
            p1[1] = (Image()->RdFrame(Sec.XLast()-ccdSec.XFirst(),iy-ccdSec.YFirst()) - p1[0])/(Sec.XLength()+1);
          }

          for (int ix = Sec.XFirst(); ix<Sec.XLast();ix++) {
            if (ix-ccdSec.XFirst() < Nx() && 
                ix-ccdSec.XFirst() >=0 )
              Image()->WrFrame(ix-ccdSec.XFirst(),iy-ccdSec.YFirst(),p1[0]+p1[1]*(ix-Sec.XFirst()+1));
          }
        }
      }
    }
}

/* -----  SpecialRedCosmetics -------------------------------------------- */
void ImageSnifs::SpecialRedCosmetics() {
  // special treatment for R channel hot lines.

  if (GetChannel()!=kRedChannel) {
    return;
  }  
  
  //
  // Handling non-standard raster
  //
  char ccdSecString[lg_name+1];
  Image()->RdDesc("CCDSEC",CHAR,lg_name+1,ccdSecString);
  if (!strcmp(ccdSecString,"Undefined")) {
    print_warning("ImageSnifs::SpecialRedCosmetics Undefined CCDSEC");
    return;
  }
  Section ccdSec(ccdSecString);
  
  // saturate value : 1 as there is a flip !!!
  double saturate;
  Image()->RdDesc("SATURAT1",DOUBLE,1,&saturate);
  
  const int kBadCols[2]={1575,1580};
  
  for (int icol=0;icol<2;icol++) {

    int ix=kBadCols[icol]-ccdSec.XFirst();
    if (ix<0) 
      continue;
    
    //
    // Get the saturating range
    //
    int ybeg=-1,yend=-1;
    ybeg = 2938-ccdSec.YFirst(); // bad pixel which sometimes does not saturate
    
    // try to find the saturating range.
    // ybeg : the first saturated point
    // yend-1 : the last saturated point
    for(int iy=0;iy<Ny();iy++) {
      if ( RdFrame(ix,iy)>saturate ) {
        if (iy<ybeg) 
          ybeg=iy;
        yend = iy+1;
      }
    } // 2943 for the 2nd case !
    if (yend<2942+icol-ccdSec.YFirst())
      yend=2942+icol-ccdSec.YFirst();
    
    // the profile approach is not robust in case of a signal : over 10 pixels,
    // the computed step can be 25% off of the real data in the worst case,
    // and an error of 4% in the best case is usual
    
    // make the correction from data computed once
    // we have to treat the first points outside the bounds also, but treatment will be different
    if (ybeg>=Ny()) 
      ybeg=Ny()-1;

    if (ybeg>0 && RdFrame(ix,ybeg)>saturate ) 
      ybeg--;
    for(int iy=ybeg;iy<=yend && iy<Ny();iy++){
      
      double prop=1;
      int bound=0;
      
      // check for boundaries
      if ( RdFrame(ix,iy)<saturate ) {
        bound=1;
        if (ix>0 && ix+1<Nx()) {
          double guess=(RdFrame(ix-1,iy)+RdFrame(ix+1,iy))/2;
          prop= (RdFrame(ix,iy)-guess)/(saturate-guess);
        }
        else // very rare case : we are at the edge in x
          prop= RdFrame(ix,iy)/saturate;
      }
      if (  RdFrame(ix,iy)>saturate && (iy-1>=0 &&  RdFrame(ix,iy-1)<saturate)
            || (iy+1<Ny() &&  RdFrame(ix,iy+1)<saturate) )
        bound=1;
      
      for (int dx=0;dx<kSpecialRed[1];dx++){
        double val=RdFrame(ix-dx-1,iy)-kSpecialCorr[1][dx]*prop;
        WrFrame(ix-1-dx,iy,val);
        if (Variance()) {
          double var;
          if (bound)
            var=Variance()->RdFrame(ix-1-dx,iy)+SQ(kSpecialCorr[1][dx]*prop*kSpecialConservative);
          else
            var=Variance()->RdFrame(ix-1-dx,iy)+SQ(kSpecialVar[1][dx]*prop*kSpecialConservative);
          
          Variance()->WrFrame(ix-dx-1,iy,var);
        }
      }      
    }
    
    //
    // now, the line itself.
    //
    
    //
    // 1st, the quick-clocked part
    //
    // that is, the beginning of the line = high y as it is flipped
    
    ImageAnalyser a;
    a.SetImage(this);
    // 4087 as the 8 last lines are spurious sometimes
    // 4047 to avoid being contaminated by 1st order light
    const int kHeight=40;
    const int kLast=4087;
    if (Ny()>kLast){    
      Section lastGround(ix+1-4,ix+1+4,kLast-kHeight,kLast);
      a.SetSection(&lastGround);
      double ground=a.Quantile(0.5);
      Section lastLine(ix+1,ix+1,kLast-kHeight,kLast);
      a.SetSection(&lastLine);
      double correction=ground-a.Quantile(0.5);
      // apply the correction
      for (int iy=yend+1;iy<Ny();iy++) {
        double val=RdFrame(ix,iy)+correction;
        WrFrame(ix,iy,val);
        // we neglect the scaling of +0.04% effect on the variance.
        // but there is a systematic uncertainty we take into account.
        // called kSpecialErrorFast
        // [the computation is however NOT accurate for rasters]
        if (Variance()){
          double var=Variance()->RdFrame(ix,iy)
            + SQ((Ny()-iy)*1.0/(Ny()-yend)*correction*kSpecialErrorFast*
                 kSpecialConservative);
          Variance()->WrFrame(ix,iy,var);
        }
      }
      WrDesc("CORRLOW",DOUBLE,1,&correction);
    }
    
    //
    // The medium part : assume we lost all the information
    //
    for (int iy=ybeg;iy<=yend&&iy<Ny();iy++) {
      if (ix>0 && ix+1<Nx()) {
        double guess=(RdFrame(ix-1,iy)+RdFrame(ix+1,iy))/2;
        WrFrame(ix,iy,guess);
      }
      if (Variance())
        Variance()->WrFrame(ix,iy,ut_big_value);
    }
    
    //
    // the hot part of the line is difficult to handle.
    //
    // the function to fit is p0+p1/(x-p2), where p2 is very close to ybeg.
    // but the best estimate we can get for p0 comes from the "end" of the line
    // that is, the beginning once it is flipped.
    // we know that an estimate using the mean of left and right makes a 
    // systematic error, and we apply it :
    // 1.0*the mean (measurment on the data from continnum lamp)
    // 
    // the varying part is not handled yet.
    
    if (ix>0 && ix+1<Nx()) {
      const int kLength=50;
      double guesses[kLength],diffs[kLength];
      for(int iy=0;iy<kLength;iy++){
        guesses[iy]=(RdFrame(ix-1,iy)+RdFrame(ix+1,iy))/2;
        diffs[iy]=RdFrame(ix,iy)-guesses[iy];
      }
      double correction=ut_median(diffs,kLength);
      //    double uncertainty=ut_median(guesses,kLength)*1.0*kSpecialConservative;

      // the correction is applied, but the variance is set to something big
      // a special control is added in order to get the numbers not too bad
      if (Variance()) {
        for (int iy=0;iy<ybeg;iy++) {
          double val=RdFrame(ix,iy)-correction;
          // check it is acceptable : shall not be outside the 
          // bounds by more than 1 sigma...
          if ( fabs(val - (RdFrame(ix+1,iy)+RdFrame(ix-1,iy))/2 )
               >  sqrt( SQ(1.0*(RdFrame(ix+1,iy)+RdFrame(ix-1,iy))/2)
                        + Variance()->RdFrame(ix,iy)))
            WrFrame(ix,iy,(RdFrame(ix+1,iy)+RdFrame(ix-1,iy))/2);
          else
            WrFrame(ix,iy,val);
          // double var= Variance()->RdFrame(ix,iy)+uncertainty*uncertainty;
          // better loose the line than doing something complex and maby wrong
          Variance()->WrFrame(ix,iy,ut_big_value);
        }
      } else {
        for (int iy=0;iy<ybeg;iy++) {
          double val=RdFrame(ix,iy)-correction;
          // check it is acceptable : shall not be outside the 
          // bounds by more than 1 sigma...
          if ( fabs(val - (RdFrame(ix+1,iy)+RdFrame(ix-1,iy))/2) 
               >  fabs(1.0*(RdFrame(ix+1,iy)+RdFrame(ix-1,iy))/2 ))
            WrFrame(ix,iy,(RdFrame(ix+1,iy)+RdFrame(ix-1,iy))/2);
          else
            WrFrame(ix,iy,val);
        }
      }
    }
  }
}

/* -----  CustomFlat -------------------------------------------- */
void ImageSnifs::CustomFlat() {
  // special treatment for R channel hot lines.

  if (GetChannel()!=kRedChannel) {
    return;
  }  

  // Exclude non-standard binning
  int nbin[2];
  Image()->RdDesc("CCDBIN",INT,2,&nbin);
  if (nbin[0]!=1 || nbin[1]!=1) {
    print_msg("ImageSnifs::CustomFlat refuses non-standard binning");
      return ;
  }

  // Non-standard raster
  char ccdSecString[lg_name+1];
  Image()->RdDesc("CCDSEC",CHAR,lg_name+1,ccdSecString);
  if (!strcmp(ccdSecString,"Undefined")) {
    print_warning("ImageSnifs::CustomFlat Undefined CCDSEC");
    return;
  }
  Section ccdSec(ccdSecString);

  for (int amp=0;amp<2;amp++){
    int xstart = amp*1024-ccdSec.XFirst();
    int xend=(amp+1)*1024-ccdSec.XFirst();
    if (xstart >= Nx() || xend < 0)
      continue;
    if (xstart<0) 
      xstart=0;
    if (xend>Nx())
      xend=Nx();
    
    for (int i=0;i<kRedNHfff[amp];i++) {
      int iy=kRedHfffLine[amp][i]-1-ccdSec.YFirst();
      if (iy<0 || iy>=Ny())
        continue;
      for (int ix=xstart;ix<xend;ix++) { 
        if (Variance()) {
          double var = Variance()->RdFrame(ix,iy)*kRedHfffVal[amp][i]*kRedHfffVal[amp][i]
            + RdFrame(ix,iy)* RdFrame(ix,iy)*kRedHfffSigma*kRedHfffSigma;
          Variance()->WrFrame(ix,iy,var);
        }
        double val = RdFrame(ix,iy)*kRedHfffVal[amp][i];
        WrFrame(ix,iy,val);
      }
    }
  }
      
}



/* ===== Hacks ======================================= */

/* ----- HackFitsKeywords ---------------------------------------------- */
void ImageSnifs::HackFitsKeywords() {
  Algo()->HackFitsKeywords(this);
}

/* ===== Processing Managment ======================================= */

/* ----- CreateVarianceFrame ------------------------------ */
void ImageSnifs::CreateVarianceFrame(char* name) {

    char varName[lg_name+1];
    strcpy(varName,name);

  if (ParanoMode()) {
    if (Variance()) {
      print_error(" ImageSnifs::CreateVarianceFrame : already a variance in %s",Name());
      return;
    }  
  } // ParanoMode
  
  if (!varName[0]) {
    ut_varname_from_imname(Name(),varName);
    if (!varName[0]) {
      print_error(" ImageSnifs::CreateVarianceFrame : no name provided",Name());
      return;
    }
  }


  // Builds the variance frame
  ImageSimple * var = new ImageSimple(*Image() , varName,FLOAT,0,IoType(),Ny());
  var->SetTo(0);
  SetVarianceFrame(var);
}

#ifdef OLD
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
  Section * Sec = Algo()->SafeOverscanStrip(this);
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
#endif


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


