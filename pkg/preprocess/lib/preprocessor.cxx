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

/* ##### Preprocessor ################################################# */

/* ===== constructor ======================================= */
Preprocessor::Preprocessor () {
  fMode = kIoPlain;
}

/* ===== Methods ======================================= */

/* ----- BuildRawBiChip -------------------------------------------------- */
BiChipSnifs* Preprocessor::BuildRawBiChip(char* name, char* outName){
  // detcom bichips are with extensions [chip00] or must contain a %d in name
  // so we check the %
  // by default, the bichip is loaded read only in totality
  BiChipSnifs* bichip;

  if (ut_is_bichip_detcom(name)) {
    bichip = new BiChipSnifs(name,"I",fMode,kIoAll);
    // not allowed because I mode
    //    bichip->SetAlgo("DETCOM");
    return bichip;
  }
  
  // we may have an otcom bichip
  ImageSnifs* image = new ImageSnifs(name);

  // otcom image is already assembled, with bias put somewhere at the end
  // a part of the bias is not OK, but I do not know yet how the bias is 
  // ordered with respect to readout time. 
  // so I assume :
  // chip0 1->1024 chip1 1024->1 bias0 1025->1057 bias 1 1057->1025

  // build the bichip and read necessary values from the image
  bichip = new BiChipSnifs();
  ImageSnifs* im[2];

  // build the names
  char imname[lg_name+1];
  strcpy(imname,outName);
  if (!outName[0])
    ut_build_tmp_name(imname,"bichip");
  
  // store some data
  char dataSecString[lg_name+1],biasSecString[lg_name+1];
  image->RdDesc("DATASEC",CHAR,lg_name+1,dataSecString);
  image->RdDesc("BIASSEC",CHAR,lg_name+1,biasSecString);
  Section dataSec(dataSecString);
  Section biasSec(biasSecString);

  // build now the images
  for (int chip=0;chip<2;chip++) {
    char extName[lg_name+1];
    
    im[chip] = new ImageSnifs(fMode,kIoAll);
    sprintf(extName,"%s[chip0%d]",imname,chip);
    im[chip]->CreateFrame(extName,(dataSec.XLength()+biasSec.XLength())/2 , dataSec.YLength());
    im[chip]->ImportHeader(image);
    double gain;
    char gainKey[lg_name+1];
    sprintf(gainKey,"CCD%dGAIN",chip);
    image->RdDesc(gainKey,DOUBLE,1,&gain);
    im[chip]->WrDesc("GAIN",DOUBLE,1,&gain);
    sprintf(dataSecString,"[%d:%d,%d:%d]",1,dataSec.XLength()/2,1,dataSec.YLength());
    
    im[chip]->WrDesc("DATASEC",CHAR,lg_name+1,dataSecString);
    sprintf(biasSecString,"[%d:%d,%d:%d]",dataSec.XLength()/2+1,dataSec.XLength()/2+biasSec.XLength()/2,1,dataSec.YLength());
    
    im[chip]->WrDesc("BIASSEC",CHAR,lg_name+1,biasSecString);
    bichip->SetChip(chip,im[chip]);
  }
  // build the data section for frame 1
  Section sec;
  sec.SetX1(dataSec.X1());
  sec.SetX2(dataSec.X2()/2);
  sec.SetY1(dataSec.Y1());
  sec.SetY2(dataSec.Y2());
  im[0]->ImportSection(image,&sec,1,1,1,1);
  sec.SetX1(dataSec.X2()/2+1);
  sec.SetX2(dataSec.X2());
  im[1]->ImportSection(image,&sec,dataSec.X2()/2,1,-1,1);
  sec.SetX1(biasSec.X1());
  sec.SetX2(biasSec.X1()-1+biasSec.XLength()/2);
  im[0]->ImportSection(image,&sec,dataSec.X2()/2+1,1,1,1);
  sec.SetX1(biasSec.X1()+biasSec.XLength()/2);
  sec.SetX2(biasSec.X2());
  im[1]->ImportSection(image,&sec,im[1]->Nx(),1,-1,1);
  
  bichip->SetAlgo("OTCOM");
  return bichip;
}

/* ----- PreprocessBias -------------------------------------------------- */
BiChipSnifs * Preprocessor::PreprocessBias(char* name, char* outName){
  // simply returns the debiased bichip
  IoMethod_t mode = fMode;
  SetIoMethod(kIoPlain);
  BiChipSnifs * bichip = BuildRawBiChip(name);
  SetIoMethod(mode);

  // build the names
  char imName[lg_name+1];
  strcpy(imName,outName);
  if (!outName[0])
    ut_build_tmp_name(imName,"bias");

  BiChipSnifs * out = new BiChipSnifs(*bichip,imName,FLOAT,1,fMode,kIoAll);
  delete bichip;
  // if algo not set, assume it is detcom
  if (!out->Chip(0)->Algo())
    out->SetAlgo("DETCOM");
  //  out->PreprocessBias();
  out->HackFitsKeywords();
  out->CreateVarianceFrame();
  out->HandleSaturation();
  out->OddEvenCorrect();
  out->AddOverscanVariance();
  out->SubstractOverscan();
  return out;
}

/* ----- PreprocessAssemble ------------------------------------------------ */
ImageSnifs * Preprocessor::PreprocessAssemble(char* name, char* outName, BiChipSnifs* bias){
  // simply returns the debiased bichip
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

