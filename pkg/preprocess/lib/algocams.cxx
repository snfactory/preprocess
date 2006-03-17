/* === Doxygen Comment ======================================= */
/*! 
 * \file          algocams.cxx
 * \copyright     (c) 2003 SNIFS-Supernova Factory Experiment
 * \date          Wed Aug  6 17:44:35 2003
 * \author        Emmanuel Gangler <e.gangler@ipnl.in2p3.fr>
 * \version       0.0
 * \brief         
 *                
 * $Id$
 **/
/* =========================================================== */

/* ----- IFU lib includes ----- */
#include "IFU_io.h"

/* ----- local includes ----- */
#include "algocams.hxx"
#include "imagesnifs.hxx"
#include "section.hxx"
#include "utils.h"

/* ##### AlgoDetcom ################################################# */

/* ===== Methods ======================================= */

/* ----- HackFitsKeywords -------------------------------------------------- */
void AlgoDetcom::HackFitsKeywords(ImageSnifs* I) {

  // get all keywords from primary header in the image
  //char primary_name[lg_name+1];
  //ut_primary_header_name(I->Name(),primary_name);
  //Anyfile * header;
  //open_primary_hd(header,primary_name,"I");
  //CP_non_std_desc(header,I->Frame());
  //close_primary_hd(header);

  // First : hack for normal raster
  char key[lg_name+1];
  I->RdDesc("RASTER",CHAR,lg_name+1,key);
  int nbin[2];
  I->RdDesc("CCDBIN",INT,2,nbin);
  if (!strcmp(key,"FULL") ){
    // The BiasSec is not correct
    int x1,x2,y1,y2;
    I->RdDesc("BIASSEC",CHAR,lg_name+1,key);
    sscanf(key,"[%d:%d,%d:%d]",&x1,&x2,&y1,&y2);
    sprintf(key,"[%d:%d,%d:%d]",x1+2,x2,1,I->Ny());
    I->WrDesc("BIASSEC",CHAR,lg_name+1,key);
    
    // the DATASEC needs help too
    I->RdDesc("DATASEC",CHAR,lg_name+1,key);
    sscanf(key,"[%d:%d,%d:%d]",&x1,&x2,&y1,&y2);
    // the truth is this one, but it is not compatible with Snf Detcom
    // so we truncate at 4096...
    //    sprintf(key,"[%d:%d,%d:%d]",4,x2+2,1,y2+5/nbin[1]);    
    sprintf(key,"[%d:%d,%d:%d]",4,x2+2,1,y2-1);
    I->WrDesc("DATASEC",CHAR,lg_name+1,key);

  } else { // rasters 
    //if (nbin[0] != 1 || nbin[1] !=1 )
    //  print_warning("ImageSnifs::HackFitsKeywords : tested only for BIN 1 1 raster");
    
    int x1,x2,y1,y2;
    I->RdDesc("BIASSEC",CHAR,lg_name+1,key);
    if (sscanf(key,"[%d:%d,%d:%d]",&x1,&x2,&y1,&y2)==4) {
       print_error("ImageSnifs::HackFitsKeywords : RASTER with Overscan not handled");
      return;
    }
    
    // OK, correct the DATASEC -- I don't know what to do with CCDSEC
    //I->RdDesc("DATASEC",CHAR,lg_name+1,key);
    //sscanf(key,"[%d:%d,%d:%d]",&x1,&x2,&y1,&y2);
    sprintf(key,"[%d:%d,%d:%d]",4,I->Nx(),1,I->Ny());
    I->WrDesc("DATASEC",CHAR,lg_name+1,key);
  }

  
  int fclass;
  if (I->RdIfDesc("FCLASS",INT,1,&fclass)<=0) {
    // The logic is that if no fclass, don't propagate.
    // no fclass -> put one
    // print_warning("AlgoDetcom::HackFitsKeywords : no FCLASS found");
    //    I->SetFClass(DONT_KNOW);
  }

  int channel;
  if (I->RdIfDesc("CHANNEL",INT,1,&channel) <=0) {
    // assume DETCOM = BLUE
    I->SetChannel(kBlueChannel);
  }

  // Then add the julian date, as it is cool to have
  if (I->RdIfDesc("JD",INT,1,&channel) <=0) {
    // set it !
    char keyVal[lg_name+1];
    int year,month,day, hour, minute;
    double second,jd;
    I->RdDesc("DATE-OBS",CHAR,lg_name+1,keyVal);
    sscanf(keyVal,"%d-%d-%d",&year,&month,&day);
    I->RdDesc("TIME-OBS",CHAR,lg_name+1,keyVal);
    sscanf(keyVal,"%d:%d:%lf",&hour,&minute,&second);
    jd = juldat(year,month,day,hour + minute/60.0 + second/3600.0);
    I->WrDesc("JD",DOUBLE,1,&jd);
  }
}

/* ----- SafeOverscanStrip -------------------------------------------------- */
Section* AlgoDetcom::SafeOverscanStrip(ImageSnifs* I) {
  // reduces the overscan region in order to have safe measurments for the 
  // algorithms

  // Get the overscan bounds from file
  char bias[lg_name+1];
  I->RdDesc("BIASSEC",CHAR,lg_name+1,bias);
  Section* Sec=new Section(bias);
  // remove the first pixel, as there is a bias pattern on it.
  // remove a few pixels (10) then, to avoid serial CTE contaminating the 
  // overscan measurments.
  Sec->SetX1(Sec->X1()+11);
  return Sec;
}

/* ##### AlgoSnfDetcom ################################################# */

/* ===== Methods ======================================= */

/* ----- HackFitsKeywords -------------------------------------------------- */
void AlgoSnfDetcom::HackFitsKeywords(ImageSnifs* I) {

  // First : hack for normal raster

  // THIS HAS TO BE CORRECTED !!!! 
  char key[lg_name+1];

  int channel;
  if (I->RdIfDesc("CHANNEL",INT,1,&channel) <=0) {
    // assume DETCOM = BLUE
    I->SetChannel(kBlueChannel);
  }

  // No hack for JD as there is MJD

  // SATURATE missing
  int saturate = 65535;
  if (I->RdIfDesc("SATURATE",INT,1,&saturate) <=0) {
    I->WrDesc("SATURATE",INT,1,&saturate);
  }
  
  // OBSTYPE missing
  char type[lg_name+1];
  if (I->RdIfDesc("OBSTYPE",CHAR,lg_name+1,type)<=0) {
    if (I->RdIfDesc("IMTYPE",CHAR,lg_name+1,type) <=0){
      print_warning("AlgoSnfDetcom::HackFitsKeywords no IMTYPE keyword found");
      type[0]=0;
    }
    I->WrDesc("OBSTYPE",CHAR,lg_name+1,type);
  }
  
  int nbin[2];
  if (I->RdIfDesc("CCDBIN",INT,2,&nbin)<=1) {
    char sum[lg_name+1];
    I->RdDesc("CCDSUM",CHAR,lg_name+1,sum);
    sscanf(sum,"%d %d",nbin,nbin+1);
    I->WrDesc("CCDBIN",INT,2,&nbin);
  }
}

/* ##### AlgoOtcom ################################################# */

/* ===== Methods ======================================= */

/* ----- HackFitsKeywords -------------------------------------------------- */
void AlgoOtcom::HackFitsKeywords(ImageSnifs* I) {

  int saturate = 65535;
  if (I->RdIfDesc("SATURATE",INT,1,&saturate) <=0) {
    I->WrDesc("SATURATE",INT,1,&saturate);
  }
  
  char type[lg_name+1];
  if (I->RdIfDesc("OBSTYPE",CHAR,lg_name+1,type)<=0) {
    if (I->RdIfDesc("IMAGETYP",CHAR,lg_name+1,type) <=0){
      print_warning("AlgoOtcom::HackFitsKeywords no IMAGETYP keyword found");
      type[0]=0;
    }
    I->WrDesc("OBSTYPE",CHAR,lg_name+1,type);
  }

  // The BiasSec contains a bad line (the last one)
  int x1,x2,y1,y2;
  char key[lg_name+1];
  I->RdDesc("BIASSEC",CHAR,lg_name+1,key);
  sscanf(key,"[%d:%d,%d:%d]",&x1,&x2,&y1,&y2);
  sprintf(key,"[%d:%d,%d:%d]",x1,x2-1,1,I->Ny());
  I->WrDesc("BIASSEC",CHAR,lg_name+1,key);
  
  int fclass;
  if (I->RdIfDesc("FCLASS",INT,1,&fclass)<=0) {
    // no fclass -> put one
    // print_warning("AlgoOtcom::HackFitsKeywords : no FCLASS found");
    //    I->SetFClass(DONT_KNOW);
  }

}

/* ----- SafeOverscanStrip -------------------------------------------------- */
Section* AlgoOtcom::SafeOverscanStrip(ImageSnifs* I) {
  // reduces the overscan region in order to have safe measurments for the 
  // algorithms

  // Get the overscan bounds from file
  char bias[lg_name+1];
  I->RdDesc("BIASSEC",CHAR,lg_name+1,bias);
  Section* Sec=new Section(bias);
  // remove the first pixel, as there is a bias pattern on it.
  // remove a few pixels (10) then, to avoid serial CTE contaminating the 
  // overscan measurments.
  // We do not know yet the bias pattern on this.
  Sec->SetX1(Sec->X1()+10);
  return Sec;
}

