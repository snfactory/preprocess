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

/* ##### AlgoDetcom ################################################# */

/* ===== Methods ======================================= */

/* ----- HackFitsKeywords -------------------------------------------------- */
void AlgoDetcom::HackFitsKeywords(ImageSnifs* I) {

  // First : hack for normal raster
  char key[lg_name+1];
  I->RdDesc("RASTER",CHAR,lg_name+1,key);
  if (!strcmp(key,"FULL")){
    // The BiasSec is not correct
    int x1,x2,y1,y2;
    I->RdDesc("BIASSEC",CHAR,lg_name+1,key);
    sscanf(key,"[%d:%d,%d:%d]",&x1,&x2,&y1,&y2);
    sprintf(key,"[%d:%d,%d:%d]",x1+2,x2,1,I->Ny());
    I->WrDesc("BIASSEC",CHAR,lg_name+1,key);
    
    // the DATASEC needs help too
    I->RdDesc("DATASEC",CHAR,lg_name+1,key);
    sscanf(key,"[%d:%d,%d:%d]",&x1,&x2,&y1,&y2);
    sprintf(key,"[%d:%d,%d:%d]",4,x2+2,1,y2+5);
    I->WrDesc("DATASEC",CHAR,lg_name+1,key);

  } else { // rasters 
    int nbin1,nbin2;
    I->RdDesc("CCDBIN1",INT,1,&nbin1);
    I->RdDesc("CCDBIN2",INT,1,&nbin2);
    if (nbin1 != 1 || nbin2 !=1 ) {
      print_error("ImageSnifs::HackFitsKeywords : only works for BIN 1 1 raster");
      return;  }
    int x1,x2,y1,y2;
    I->RdDesc("BIASSEC",CHAR,lg_name+1,key);
    if (sscanf(key,"[%d:%d,%d:%d]",&x1,&x2,&y1,&y2)==4) {
       print_error("ImageSnifs::HackFitsKeywords : RASTER with Overscan not handled");
      return;
    }
    
    // OK, correct the DATASEC -- I don't know what to do with CCDSEC
    I->RdDesc("DATASEC",CHAR,lg_name+1,key);
    sscanf(key,"[%d:%d,%d:%d]",&x1,&x2,&y1,&y2);
    sprintf(key,"[%d:%d,%d:%d]",4,x2-1,1,y2-1);
    I->WrDesc("DATASEC",CHAR,lg_name+1,key);
  }

  
  int fclass;
  if (I->RdIfDesc("FCLASS",INT,1,&fclass)<=0) {
    // no fclass -> put one
    I->SetFClass(DONT_KNOW);
  }
}

/* ===== ifdef ================================================== */
#ifdef FOR_FUTURE_WHY_NOT
/* ----- OddEvenCorrect -------------------------------------------------- */
void AlgoDetcom::OddEvenCorrect(ImageSnifs* I) {
  // odd-even substraction
  //

  double param[2];
  
  if (I->ParanoMode()) {
    // check overscan was not already substracted
    // returns 0 if not done !
    if ( I->RdIfDesc("OEPARAM",DOUBLE, 2, param) > 0 ){
      print_error("Odd-Even already substracted for %s\n",Name());
      print_error("Nothing Done%s\n",Name());
      return;
    }
  }

  // Get the overscan bounds from file
  Section* Sec=SafeoverscanStrip(I);

  // the magic 2.0 parameter is because of CTE on CCD EEV#1 (blue)
  // this could also duplicate the magic effect of the reduced overscan section
  I->Image()->OddEvenCorrect(Sec,param,2.0);
  delete Sec;

  I->WrDesc("OEPARAM",DOUBLE, 2, param);
}


/* ----- AddOverscanVariance ------------------------------ */
void AlgoDetcom::AddOverscanVariance(ImageSnifs* I) {

  int ovscNoise;
  if (I->ParanoMode()) {
    if (I->Variance()->RdIfDesc("OVSCNOIS",INT, 1, &ovscNoise) > 0 
     && ovscNoise ) {
       print_error(" ImageSnifs::AddOverscanVariance already done in %s",Name());
      return;
      }
  }

  // Fills with analysis of RMS strip
  Section * Sec = SafeOverscanStrip(I);
  double rms = I->Image()->OverscanRms(Sec,0);
  WrDesc("RDNOISE",DOUBLE,1,&rms);
  I->Variance()->Add(rms*rms);
  delete Sec;

}

#endif

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

/* ##### AlgoOtcom ################################################# */

/* ===== Methods ======================================= */

/* ----- HackFitsKeywords -------------------------------------------------- */
void AlgoOtcom::HackFitsKeywords(ImageSnifs* I) {

  int saturate = 65535;
  if (I->RdIfDesc("SATURATE",INT,1,&saturate) <=0) {
    print_warning("AlgoOtcom::HackFitsKeywords no SATURATE keyword found");
    I->WrDesc("SATURATE",INT,1,&saturate);
  }
  
  char type[lg_name+1];
  if (I->RdIfDesc("IMAGETYP",CHAR,lg_name+1,type) <=0){
    print_warning("AlgoOtcom::HackFitsKeywords no IMAGETYP keyword found");
    type[0]=0;
  }
  I->WrDesc("OBSTYPE",CHAR,lg_name+1,type);

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
    I->SetFClass(DONT_KNOW);
  }
}

/* ===== ifdef ================================================== */
#ifdef FOR_FUTURE_WHY_NOT
/* ----- OddEvenCorrect -------------------------------------------------- */
void AlgoOtcom::OddEvenCorrect(ImageSnifs* I) {
  // odd-even substraction
  //

  double param[2];
  
  if (I->ParanoMode()) {
    // check overscan was not already substracted
    // returns 0 if not done !
    if ( I->RdIfDesc("OEPARAM",DOUBLE, 2, param) > 0 ){
      print_error("Odd-Even already substracted for %s\n",Name());
      print_error("Nothing Done%s\n",Name());
      return;
    }
  }

  // Get the overscan bounds from file
  Section* Sec=SafeoverscanStrip(I);

  // the magic 2.0 parameter is because of CTE on CCD EEV#1 (blue)
  // this could also duplicate the magic effect of the reduced overscan section
  I->Image()->OddEvenCorrect(Sec,param,2.0);
  delete Sec;

  I->WrDesc("OEPARAM",DOUBLE, 2, param);
}


/* ----- AddOverscanVariance ------------------------------ */
void AlgoOtcom::AddOverscanVariance(ImageSnifs* I) {

  int ovscNoise;
  if (I->ParanoMode()) {
    if (I->Variance()->RdIfDesc("OVSCNOIS",INT, 1, &ovscNoise) > 0 
     && ovscNoise ) {
       print_error(" ImageSnifs::AddOverscanVariance already done in %s",Name());
      return;
      }
  }

  // Fills with analysis of RMS strip
  Section * Sec = SafeOverscanStrip(I);
  double rms = I->Image()->OverscanRms(Sec,0);
  WrDesc("RDNOISE",DOUBLE,1,&rms);
  I->Variance()->Add(rms*rms);
  delete Sec;

}

#endif

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

