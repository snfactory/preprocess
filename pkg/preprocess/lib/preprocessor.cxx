/* === Doxygen Comment ======================================= */
/*! 
 * \file          preprocessor.cxx
 * \copyright     (c) 2003 SNIFS-Supernova Factory Experiment
 * \date          Wed Aug  6 17:44:35 2003
 * \author        Emmanuel Gangler <e.gangler@ipnl.in2p3.fr>
 * \version       0.0
 * \brief         
 *                
 * $Id$
 **/
/* =========================================================== */

/* ----- local includes ----- */
#include "preprocessor.hxx"
#include "bichip.hxx"
#include "imagesnifs.hxx"
#include "section.hxx"
#include "utils.h"
#include "algocams.hxx"
#include "overscan.hxx"

/* Note about the fast mode : it suppresses OddEven correction and Variance computation */

/* ##### Preprocessor ################################################# */

/* ===== constructor/destructor ======================================= */

/* ----- constructor ---------------------------------------- */
Preprocessor::Preprocessor () {
  fMode = kIoPlain;
  fOverscanSnifs = new OverscanSnifs();
  fOverscanRescue = new OverscanFromData();
  fOverscan = fOverscanSnifs; // default current overscan
  fFast=0;
  fAllImage=0;
}

/* ----- destructor ---------------------------------------- */
Preprocessor::~Preprocessor () {
  if (fOverscanSnifs)
    delete fOverscanSnifs;
  if (fOverscanRescue)
    delete fOverscanRescue;
}

/* ===== Setters ======================================= */

/* ----- SetOverscanAuto -------------------------------------------------- */
void Preprocessor::SetOverscanAuto(ImageSnifs* Image){
  /* auto - toggling of the overscan between regular / rescue */
  fOverscan = fOverscanSnifs;

  if (Image->HasOverscan()) {
    fOverscan = fOverscanRescue;
  }
  
}


/* ===== Methods ======================================= */

/* ----- BuildRawBiChip -------------------------------------------------- */
BiChipSnifs* Preprocessor::BuildRawBiChip(char* name, char* outName){
  // returns the original bichip for detcom image
  // reconstructs a bichip from raw data for otcom image.

  // detcom bichips are with extensions [chip00] or must contain a %d in name
  // so we check the %
  // by default, the bichip is loaded read only in totality

  // build the names
  char imname[lg_name+1];
  if (!outName[0]) {
    strcpy(imname,name);
    ut_build_tmp_name(imname,"bichip");
  }
  else
    strcpy(imname,outName);

  BiChipSnifs* bichip;

  // First case : DETCOM

  // We could have made a copy of the bichip at this stage,
  // but we want to allow for an utilitary routine converting
  // otcom->detcom without time overhead in case of an original
  // detcom image
  if (ut_is_bichip_detcom(name)) {
    bichip = new BiChipSnifs(name,"I",fMode,kIoAll);
    // not allowed because I mode
    //    bichip->SetAlgo("DETCOM");
    // no new image asked for
    if (!outName[0] || !strcmp(name,outName))
      return bichip;
    else {
      BiChipSnifs* out = new BiChipSnifs(*bichip,outName,FLOAT,1,fMode,kIoAll);
      out->SetAlgo("DETCOM");
      delete bichip;
      return out;
    }
  }
  
  // Default : OTCOM
  //
  // otcom image is already assembled, with bias put somewhere at the end
  // so I assume :
  // chip0 1->1024 chip1 1024->1 bias0 1025->1057 bias 1 1057->1025

  ImageSnifs* image = new ImageSnifs(name);

  // read necessary values for the image
  char dataSecString[lg_name+1],biasSecString[lg_name+1];
  image->RdDesc("DATASEC",CHAR,lg_name+1,dataSecString);
  image->RdDesc("BIASSEC",CHAR,lg_name+1,biasSecString);
  Section dataSec(dataSecString);
  Section biasSec(biasSecString);

  // detect if it is a raster : 
  // in case of a raster, substract 1 column of bad data. 
  // (last one 'before'(=in 1 channel readout sequence) the overscan )
  int isRaster=0;
  int nAmp;
  char ccdSecString[lg_name+1];
  image->RdDesc("CCDSEC",CHAR,lg_name+1,ccdSecString);
  image->RdDesc("CCDNAMP",INT,1,&nAmp);
  Section ccdSec(ccdSecString);
  if (ccdSec.XLength()!=1024*nAmp) {
    print_msg("Preprocessor::BuildRawBiChip detected a Raster image");
    isRaster=1;
  }

  // build the bichip
  int newNamp=nAmp;
  if (nAmp != 2 && GetAllImage())
    bichip = new BiChipSnifs(nAmp);
  else {
    bichip = new BiChipSnifs(2);
    newNamp=2;
  }
  
  ImageSnifs* im[bichip->NChips()];

  // build now the images and the headers
  for (int chip=0;chip<bichip->NChips();chip++) {
    char extName[lg_name+1];
    int newDataXLength=dataSec.XLength()/nAmp - isRaster;
    int newBiasXLength=biasSec.XLength()/nAmp;

    im[chip] = new ImageSnifs(fMode,kIoAll);
    sprintf(extName,"%s[chip0%d]",imname,chip);
    im[chip]->CreateFrame(extName,newDataXLength + newBiasXLength, dataSec.YLength());
    im[chip]->ImportHeader(image);
    double gain;
    char gainKey[lg_name+1];
    sprintf(gainKey,"CCD%dGAIN",chip);
    image->RdDesc(gainKey,DOUBLE,1,&gain);
    im[chip]->WrDesc("GAIN",DOUBLE,1,&gain);
    sprintf(dataSecString,"[%d:%d,%d:%d]",1,newDataXLength,1,dataSec.YLength());
    if (newNamp != nAmp)
      im[chip]->WrDesc("CCDNAMP",INT,1,&newNamp);
    
    im[chip]->WrDesc("DATASEC",CHAR,lg_name+1,dataSecString);
    sprintf(biasSecString,"[%d:%d,%d:%d]",newDataXLength+1,newDataXLength+newBiasXLength,1,dataSec.YLength());
    
    im[chip]->WrDesc("BIASSEC",CHAR,lg_name+1,biasSecString);
    bichip->SetChip(chip,im[chip]);

    // fill the image
    Section sec;
    if (chip%2==0) {
      sec.SetX1(dataSec.X1()+(chip*dataSec.XLength())/nAmp);
      sec.SetX2(sec.X1() + newDataXLength - 1 );
      sec.SetY1(dataSec.Y1());
      sec.SetY2(dataSec.Y2());
      im[chip]->ImportSection(image,&sec,1,1,1,1);
      sec.SetX1(biasSec.X1()+(chip*biasSec.XLength())/nAmp);
      sec.SetX2(sec.X1() + newBiasXLength - 1);
      sec.SetY1(biasSec.Y1());
      sec.SetY2(biasSec.Y2());
      im[chip]->ImportSection(image,&sec,newDataXLength+1,1,1,1);
    } else {
      sec.SetX1(dataSec.X1()+(chip*dataSec.XLength())/nAmp + isRaster);
      sec.SetX2(sec.X1() + newDataXLength - 1);
      sec.SetY1(dataSec.Y1());
      sec.SetY2(dataSec.Y2());
      im[chip]->ImportSection(image,&sec,newDataXLength,1,-1,1);
      sec.SetX1(biasSec.X1()+(chip*biasSec.XLength())/nAmp);
      sec.SetX2(sec.X1() + newBiasXLength - 1);
      sec.SetY1(biasSec.Y1());
      sec.SetY2(biasSec.Y2());
      im[chip]->ImportSection(image,&sec,newDataXLength+newBiasXLength,1,-1,1);
    }
    
  }

  // needs the images 
  bichip->SetAlgo("OTCOM");
  delete image;
  
  return bichip;
}

/* ----- PreprocessOverscan------------------------------------------------- */
BiChipSnifs * Preprocessor::PreprocessOverscan(char* name, char* outName){
  /* handles all the information relevant to the overscan */

  // Preliminary : getting a bichip in IO mode

  // build a temporary name if needed
  char imName[lg_name+1];
  if (!outName[0]) {
    strcpy(imName,name);
    ut_build_tmp_name(imName,"bias");
  } else
    strcpy(imName,outName);

  BiChipSnifs * out = BuildRawBiChip(name,imName);

  // keyword hacking

  // Detcom image -> make a special header hack
  if (!out->Chip(0)->Algo()->GetId() == kDetcom) {
    char primary_name[lg_name+1];
    ut_primary_header_name(out->Chip(0)->Name(),primary_name);
    out->HackFitsKeywords(primary_name);
  } else { 
    // keywords hacking
    out->HackFitsKeywords();
  }
    
  //
  // OK, we have now a working copy of the image !
  //


  // variance creation
  if (!FastMode()) {
    out->CreateVarianceFrame();
    out->HandleSaturation();
  }
  
  // overscan substraction
  if (out->Chip(0)->HasOverscan()) {
    // normal exposure with an overscan
    if (FastMode())
      fOverscanSnifs->SetOddEven(0);
    else
      fOverscanSnifs->SetOddEven(1);
    fOverscanSnifs->Correct(out);
  }
  // rescue procedure
  else {
    fOverscanRescue->Correct(out);
  }

  out->UpdateFClass();
  return out;
}

/* ----- PreprocessAssemble ------------------------------------------------ */
ImageSnifs * Preprocessor::PreprocessAssemble(char* name, char* outName, ImageSnifs* bias){

  // simply returns the debiased bichip
  // first builds teh debiassed bichip
  // The opening mode shall be the standard one (temporary creation)
  IoMethod_t mode = fMode;
  SetIoMethod(kIoPlain);
  
  BiChipSnifs * bichip = PreprocessOverscan(name);
  // no hack for the variating gain - no good hack found anyway!
  //bichip->HackGainRatio();

  SetIoMethod(mode);
  ImageSnifs *out = bichip->Assemble(outName,fMode,kIoAll);
  delete bichip;
  // remove bad lines
  out->HandleCosmetics();
  // Special request from Yannick
  if (FastMode())
    out->CheatCosmetics();
  if (bias) out->SubstractBias(bias);
  return out;
  
}


/* ----- PreprocessDark ------------------------------------------------ */
ImageSnifs * Preprocessor::PreprocessDark(char* name, char* outName,ImageSnifs* bias){
  // simply returns the debiased bichip
  
  ImageSnifs *out = PreprocessAssemble(name,outName,bias);
  out->AddPoissonNoise();
  return out;
}

/* ----- PreprocessFlat ------------------------------------------------ */
ImageSnifs * Preprocessor::PreprocessFlat(char* name, char* outName,ImageSnifs* bias, ImageSnifs* dark){

  ImageSnifs* out = PreprocessDark(name,outName,bias);
  if (dark) 
    out->SubstractDark( dark );
  out->BuildFlat();
  return out;
}

/* ----- Preprocess ------------------------------------------------ */
ImageSnifs* Preprocessor::Preprocess(char* name, char* outName,ImageSnifs *bias,ImageSnifs *dark,ImageSnifs* flat) {
  ImageSnifs* out = PreprocessDark(name, outName,bias);
  if (dark) 
    out->SubstractDark( dark);
  if (flat) 
    out->ApplyFlat(flat);
  return out;
}

