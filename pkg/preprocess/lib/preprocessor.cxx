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
    return bichip;
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



#ifdef OLD
  // Fill the images with the data content
  Section sec;
  sec.SetX1(dataSec.X1());
  sec.SetX2(dataSec.X2()/2-isRaster);
  sec.SetY1(dataSec.Y1());
  sec.SetY2(dataSec.Y2());
  im[0]->ImportSection(image,&sec,1,1,1,1);
  sec.SetX1(dataSec.X2()/2+1+isRaster);
  sec.SetX2(dataSec.X2());
  im[1]->ImportSection(image,&sec,dataSec.X2()/2,1,-1,1);
  sec.SetX1(biasSec.X1());
  sec.SetX2(biasSec.X1()-1+biasSec.XLength()/2);
  im[0]->ImportSection(image,&sec,dataSec.X2()/2+1,1,1,1);
  sec.SetX1(biasSec.X1()+biasSec.XLength()/2);
  sec.SetX2(biasSec.X2());
  im[1]->ImportSection(image,&sec,im[1]->Nx(),1,-1,1);
#endif
  // needs the images 
  bichip->SetAlgo("OTCOM");
  return bichip;
}

/* ----- PreprocessBias -------------------------------------------------- */
BiChipSnifs * Preprocessor::PreprocessBias(char* name, char* outName){
  // simply returns the debiased bichip

  BiChipSnifs * bichip = BuildRawBiChip(name);
  BiChipSnifs * out;

  // Preliminary : getting an IO copy

  // Detcom image -> copy to out to be able to set keywords, etc...
  // if algo not set, assume it is detcom
  if (!bichip->Chip(0)->Algo()) {
    // build the names
    char imName[lg_name+1];
    if (!outName[0]) {
      strcpy(imName,name);
      ut_build_tmp_name(imName,"bias");
    }
    else
      strcpy(imName,outName);

    out = new BiChipSnifs(*bichip,imName,FLOAT,1,fMode,kIoAll);
    out->SetAlgo("DETCOM");
    delete bichip;
  } else { 
    // otcom image : no need for an additional copy (was done in BuildRawBiChip )
    out = bichip;
  }
    
  //
  // OK, we have now a working copy of the image !
  //

  // keywords hacking
  out->HackFitsKeywords();

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
ImageSnifs * Preprocessor::PreprocessAssemble(char* name, char* outName, BiChipSnifs* bias){

  // simply returns the debiased bichip
  // first builds teh debiassed bichip
  // The opening mode shall be the standard one (temporary creation)
  IoMethod_t mode = fMode;
  SetIoMethod(kIoPlain);
  
  BiChipSnifs * bichip = PreprocessBias(name);
  if (bias) bichip->SubstractBias(bias);
  // no hack for the variating gain - no good hack found anyway!
  //bichip->HackGainRatio();

  SetIoMethod(mode);
  ImageSnifs *out = bichip->Assemble(outName,fMode,kIoAll);
  delete bichip;
  return out;
  
}

/* ----- PreprocessDark ------------------------------------------------ */
ImageSnifs * Preprocessor::PreprocessDark(char* name, char* outName,BiChipSnifs* bias){
  // simply returns the debiased bichip
  
  ImageSnifs *out = PreprocessAssemble(name,outName,bias);
  out->AddPoissonNoise();
  return out;
}

/* ----- PreprocessFlat ------------------------------------------------ */
ImageSnifs * Preprocessor::PreprocessFlat(char* name, char* outName,BiChipSnifs* bias, ImageSnifs* dark){

  ImageSnifs* out = PreprocessDark(name,outName,bias);
  if (dark) 
    out->SubstractDark( dark );
  out->BuildFlat();
  return out;
}

/* ----- Preprocess ------------------------------------------------ */
ImageSnifs* Preprocessor::Preprocess(char* name, char* outName,BiChipSnifs *bias,ImageSnifs *dark,ImageSnifs* flat) {
  ImageSnifs* out = PreprocessDark(name, outName,bias);
  if (dark) 
    out->SubstractDark( dark);
  if (flat) 
    out->ApplyFlat(flat);
  return out;
}

