/* === Doxygen Comment ======================================= */
/*! 
 * \file          bichip.cxx
 * \copyright     (c) 2003 SNIFS-Supernova Factory Experiment
 * \date          Wed Aug  6 17:44:35 2003
 * \author        Emmanuel Gangler <e.gangler@ipnl.in2p3.fr>
 * \version       0.0
 * \brief         
 *                
 * $Id$
 **/
/* =========================================================== */

#include "bichip.hxx"
#include "imagesnifs.hxx"
#include "section.hxx"
#include "utils.h"
#include "algocams.hxx"

/* ##### Data ################################################# */

static double kGainBlue[2]={0.773,0.744};
static double kGainRed[2]={0.736,0.792};
static double kGainPhot[4]={1.618,1.576,1.51,1.52};

/* ##### BiChipSnifs ################################################# */

/* ===== constructor/destructor ======================================= */

/* ----- BiChipSnifs -------------------------------------------------- */
BiChipSnifs::BiChipSnifs(int Nchips) {
  // Default constructor :
  fNChips = Nchips;
  fChip = new ImageSnifs* [fNChips];
  for (int i=0;i<NChips();i++){
    fChip[i]=0;
  }
}


/* ----- BiChipSnifs -------------------------------------------------- */
BiChipSnifs::BiChipSnifs(char* ChipNameRecipee,char* mode,IoMethod_t Method, int MParam) {
  // Default constructor :
  // The recipee shall contain a %d, which will be either replaced by a 0 or 1
  char modName[lg_name+1];
  strcpy(modName,ChipNameRecipee);
  CheckNameRecipee(modName);
  char chipname[lg_name+1];
  // first count the number of chips
  int Nchip=0;
  sprintf(chipname, modName, Nchip);
  while (exist(chipname)) {
    Nchip++;
    sprintf(chipname, modName, Nchip);
  } 
  fNChips = Nchip;
  fChip = new ImageSnifs* [fNChips];  
  for (int chip=0;chip<NChips();chip++) {
    sprintf(chipname, modName, chip);
    fChip[chip] = new ImageSnifs(chipname,mode,Method,MParam);
  }
}

/* ----- BiChipSnifs -------------------------------------------------- */
BiChipSnifs::BiChipSnifs(const BiChipSnifs &Father,char* NewNameRecipee,short newtype,int copydata,IoMethod_t Method, int MParam) {
  // The copy constructor maps the ImageSnifs copy constructor
  // The recipee shall contain a %d, which will be either replaced by a 0 or 1
  char modName[lg_name+1];
  strcpy(modName,NewNameRecipee);
  CheckNameRecipee(modName);
  fNChips = Father.NChips();
  fChip = new ImageSnifs* [fNChips];  
  for (int chip=0;chip<NChips();chip++) {
    char chipName[lg_name+1];
    sprintf(chipName,modName,chip);
    fChip[chip]=new ImageSnifs((*Father.fChip[chip]),chipName,newtype,copydata,Method,MParam);
  }
}

/* ----- ~BiChipSnifs -------------------------------------------------- */
BiChipSnifs::~BiChipSnifs() {
  for (int chip=0 ; chip<NChips() ; chip++) {
    if ( fChip[chip] )
      delete fChip[chip];
  }
}

/* ===== Utilities ====================================================== */

/* ----- CheckNameRecipee ----------------------------------------------- */
void BiChipSnifs::CheckNameRecipee(char* ChipNameRecipee) {
  // checks wether the recepee contains a %d
  // else appends [chip0%d] to the recipee

  char* test;
  // everything OK
  if ((test =  strstr(ChipNameRecipee,"%d"))) {
    return;
  }
  // Shall not contain a []
  if (( test = strrchr(ChipNameRecipee,'[') )) {
    print_error("BiChipSnifs::CheckNameRecipee : %s contains an extension but no %d");
    return;
  }
  sprintf(ChipNameRecipee,"%s[chip0%%d]",ChipNameRecipee);
}

/* ----- SetNLines ----------------------------------------------- */
void BiChipSnifs::SetNLines(int NLines) {
  for (int iChip=0;iChip<NChips();iChip++) {
    Chip(iChip)->SetNLines(NLines);
  }
}

/* ----- SetParanoMode ----------------------------------------------- */
void BiChipSnifs::SetParanoMode(bool Mode) {
  for (int iChip=0;iChip<NChips();iChip++) {
    Chip(iChip)->SetParanoMode(Mode);
  }
  
}

/* ----- SetParanoMode ----------------------------------------------- */
void BiChipSnifs::SetAlgo(char* Soft) {
  for (int iChip=0;iChip<NChips();iChip++) {
    Chip(iChip)->SetAlgo(Soft);
  }
  
}


/* ===== Methods ====================================================== */

void BiChipSnifs::UpdateFClass() {
  for (int chip=0;chip<NChips();chip++) {
    fChip[chip]->UpdateFClass();
  }
  
}


/* ----- Overscan ----------------------------------------------------- */
#ifdef OLD
void BiChipSnifs::SubstractOverscan() {
  // Substracts the overscan from NChips() chips
  // We want a crash if chips not already set !
  for (int chip=0;chip<NChips();chip++) {
    fChip[chip]->SubstractOverscan();
  }
}

/* ----- odd-Even ----------------------------------------------------- */
void BiChipSnifs::OddEvenCorrect() {
  // substracts the odd-even from NChips() chips
  for (int chip=0;chip<NChips();chip++) {
    fChip[chip]->OddEvenCorrect();
  }
}
#endif 

/* ----- substract bias ------------------------------------------------ */
void BiChipSnifs::SubstractBias( BiChipSnifs* Bias) {
  // substracts the bias from NChips() chips
  for (int chip=0;chip<NChips();chip++) {
    fChip[chip]->SubstractBias(Bias->fChip[chip]);
  }
}

/* ----- substract dark ------------------------------------------------ */
void BiChipSnifs::SubstractDark( BiChipSnifs* Dark) {
  // substracts the dark from NChips() chips
  for (int chip=0;chip<NChips();chip++) {
    fChip[chip]->SubstractDark(Dark->fChip[chip]);
  }
}

/* ----- Assemble ----------------------------------------------------- */
ImageSnifs* BiChipSnifs::Assemble(char* ImageName,IoMethod_t Io, int Nlines) {

  // Assembles the NChips() chips in 1 image
  // if the image has no extension indication, creates a 
  // [image] and a [variance] frame if the bichip is compatible
  // with variance.

  // parano checks
  if (fChip[0]->ParanoMode() || fChip[1]->ParanoMode()) {
    // images shall be bias-substracted for gain correction guess
    int elecalib0, elecalib1;
    if ( (fChip[0]->RdIfDesc("ELECALIB",INT, 1, &elecalib0)>0 && elecalib0 )
         || (fChip[1]->RdIfDesc("ELECALIB",INT, 1, &elecalib1)>0 && elecalib1 ) ) {
      print_error("BiChipSnifs::Assemble : Electron calib already done\n");
      return 0;
    }
    // many more parano checks possible : images shall be identical up to chip name !
  }
  
  // synthesis of true image name
  char imageN[lg_name+1];
  if(strchr(ImageName,'['))
    strcpy(imageN,ImageName);
  else
    sprintf(imageN,"%s[image]",ImageName);

  // synthesis of variance name -> usefull to know if it exists
  char varname[lg_name+1];
  int hasVar=1;
  for (int chip=0;chip<NChips();chip++){
    hasVar = hasVar & (fChip[chip]->Variance()!=0);
  }
  if (hasVar )
    ut_varname_from_imname(imageN,varname);
  else
    varname[0]='\0';


  // Raw constructor as we will allocate the rest on the spot
  ImageSnifs* compound = new ImageSnifs(Io,Nlines);
  ImageSimple* variance=0;
  if (varname[0]) {
    variance =  new ImageSimple();
    compound->SetVarianceFrame(variance);
  }

  //
  // Create the frame :
  //
  char key[lg_name+1];
  fChip[0]->RdDesc("DATASEC",CHAR,lg_name+1,key);
  Section * Sec = new Section(key);
  // Some parano checks 
  // We could also check the entire header !!!
  char key2[lg_name+1];
  for (int chip=1;chip<NChips();chip++) {
    fChip[chip]->RdDesc("DATASEC",CHAR,lg_name+1,key2);
    if (strcmp(key,key2)) {
      print_error("BiChipSnifs::Assemble Chips of different data size");
      return 0;
    }
  }

  compound->CreateFrame(imageN, Sec->XLength()*NChips() ,Sec->YLength() );
  if (variance)
    variance->CreateFrame(varname, Sec->XLength()*NChips() ,Sec->YLength() );
  
  // Puts header
  compound->ImportHeader(fChip[0]);
  if (variance)
    variance->ImportHeader(fChip[0]->Variance() );

  //
  // Puts the image data from chips
  //
  char keyName[lg_name+1];
  double gain,rdnoise[NChips()];
  for (int chip=0;chip<NChips();chip++) {

  // Note that R channel is flipped (180deg rotation) ... it is economic to do it here
    int xDir=1;
    int yDir=1;
    int xStart=1;
    int yStart=1;
    if (fChip[0]->GetChannel() == kRedChannel){
      xDir=-1;
      yDir=-1;
      xStart = Sec->XLength()*NChips();
      yStart = Sec->YLength();
    }
    
    fChip[chip]->RdDesc("GAIN",DOUBLE,1,&gain);


    // ImportSection takes care of variance !
    if (!(chip%2))
      compound->ImportSection(fChip[chip],Sec,xStart+(chip*Sec->XLength())*xDir,yStart,xDir,yDir,gain);
    else
      compound->ImportSection(fChip[chip],Sec,xStart-xDir + (1+chip)*Sec->XLength()*xDir,yStart,-xDir,yDir,gain);
    sprintf(keyName,"CCD%dGAIN",chip);
    compound->WrDesc(keyName,DOUBLE,1,&gain);
    //   if (fChip[chip]->RdIfDesc("NSATU",INT,1,&chipSatu) > 0)
    //  nSatu += chipSatu;

    // propagates the RdNoises
    fChip[chip]->RdDesc("RDNOISE",DOUBLE,1,rdnoise+chip);
    rdnoise[chip]*=gain;
  }
  
  // ... remains to update DATASEC and remove BIASSEC
  // We also set the ELECALIB flag
  int  eleCalib=1;
  sprintf(key,"[%d:%d,%d:%d]",1,Sec->XLength()*NChips(),1,Sec->YLength());
  compound->WrDesc("DATASEC",CHAR,lg_name+1,key);
  compound->DeleteDesc("BIASSEC");
  compound->WrDesc("ELECALIB",INT,1,&eleCalib);
  if (variance) {
    variance->WrDesc("DATASEC",CHAR,lg_name+1,key);
    variance->DeleteDesc("BIASSEC");
  }
  
  compound->WrDesc("RDNOISE",DOUBLE,NChips(),rdnoise);

  // and update the number of saturating pixels
  // compound->WrDesc("NSATU",INT,1,&nSatu);
  // well ... bad idea

  delete Sec;
  return compound;
}

/* ----- GuessGainRatio ------------------------------------------------ */
double BiChipSnifs::GuessGainRatio(Section* S) {
  // the ratio is gain(0)

  double val[2];
  double line[S->XLength()];
  double column[S->YLength()];

  if (NChips()!=2) {
    print_error("BiChipSnifs::GuessGainRatio is not designed for more than 2 chips");
    return 1;
  }

  for (int iy=S->YFirst();iy<S->YLast();iy++){
    // note the order : we do the median of the ratios
    for (int chip=0;chip<2;chip++) {
    
      for (int ix=S->XFirst();ix<S->XLast();ix++) {
        line[ix-S->XFirst()] = fChip[chip]->RdFrame(ix,iy);
      }
      val[chip]= ut_mean(line, S->XLength());
    } //chip
    column[iy - S->YFirst()] = val[0]/val[1];
  }
  return ut_median(column,S->YLength());
}

/* ----- CreateVarianceFrame ------------------------------------------------ */
void BiChipSnifs::CreateVarianceFrame(char* VarianceNameRecipee) 
{
  // create the variance frame
  // if the recipee contains a %d -> names of teh variance
  // else if the name is a plain .fits -> crates a err00 and a err01 extension
  // else tries the classical varname from imname (i.e. : var00 and var01
  // on the same file than the image, if its extension is chip00 / chip01)

  // Check the recipee (see CheckNameRecipee)
  if (VarianceNameRecipee[0] && !(strstr(VarianceNameRecipee,"%d"))) {
    if (( strrchr(VarianceNameRecipee,'[') )) {
      print_error("BiChipSnifs::CheckNameRecipee : %s contains an extension but no %d");
      return;
    }
    sprintf(VarianceNameRecipee,"%s[var0%%d]",VarianceNameRecipee);
  }

  for (int chip=0;chip<NChips();chip++) {
    // builds actual name
    char varName[lg_name+1];
    if (VarianceNameRecipee[0])
      sprintf(varName,VarianceNameRecipee,chip);
    else
      varName[0]='\0';
    // do it
    fChip[chip]->CreateVarianceFrame(varName);
  }
}

#ifdef OLD
/* ----- AddOverscanVariance ---------------------------------------- */
void BiChipSnifs::AddOverscanVariance() {
  // substracts the odd-even from 2 chips
  for (int chip=0;chip<NChips();chip++) {
    fChip[chip]->AddOverscanVariance();
  }
  
}
#endif

/* ----- HandleSaturation ---------------------------------------- */
void BiChipSnifs::HandleSaturation() {
  // substracts the odd-even from 2 chips
  for (int chip=0;chip<NChips();chip++) {
    fChip[chip]->HandleSaturation();
  } 
}


 
/* ----- PreprocessBias ------------------------------------------------ */
#ifdef OLD
void BiChipSnifs::PreprocessBias() {
  CreateVarianceFrame();
  HandleSaturation();
  HackFitsKeywords();
  OddEvenCorrect();
  AddOverscanVariance();
  SubstractOverscan();
}

/* ----- PreprocessDark ------------------------------------------------ */
ImageSnifs* BiChipSnifs::PreprocessAssemble(char* OutName,BiChipSnifs *bias) {
  PreprocessBias();
  if (bias) SubstractBias(bias);
  HackGainRatio();
  ImageSnifs *out = Assemble(OutName);
  return out;
}

/* ----- PreprocessDark ------------------------------------------------ */
ImageSnifs* BiChipSnifs::PreprocessDark(char* OutName,BiChipSnifs *bias) {
  ImageSnifs *out = PreprocessAssemble(OutName,bias);
  out->AddPoissonNoise();
  return out;
}

/* ----- PreprocessFlat ------------------------------------------------ */
ImageSnifs* BiChipSnifs::PreprocessFlat(char* OutName,BiChipSnifs *bias,ImageSnifs *dark) {
  ImageSnifs* out = PreprocessDark(OutName,bias);
  if (dark) 
    out->SubstractDark( dark );
  out->BuildFlat();
  return out;
}

/* ----- Preprocess ------------------------------------------------ */
ImageSnifs* BiChipSnifs::Preprocess(char* OutName,BiChipSnifs *bias,ImageSnifs *dark,ImageSnifs* flat) {
  ImageSnifs* out = PreprocessDark(OutName,bias);
  if (dark) 
    out->SubstractDark( dark );
  if (flat) 
    out->ApplyFlat(flat);
  return out;
}
#endif
  

/* ===== Ugly part ====================================================== */

/* ----- HackFitsSecKeywords()  ------------------------------------------- */
void  BiChipSnifs:: HackFitsKeywords()  {
  for (int chip=0;chip<NChips();chip++) {
    fChip[chip]->HackFitsKeywords();
  }
  
  // Easier place for the channel hack
  int channel=kBlueChannel;
  int nAmp;

  if (fChip[0]->Algo()->GetId() == kOtcom) {
    fChip[0]->RdDesc("CCDNAMP",INT,1,&nAmp);
    // assume OTCOM 2amps = RED
    // OTCOM 4 amps 2 chips = Photomatric
    // OTCOM 4 amps 4 chan
    if (nAmp==2)
      channel = kRedChannel;
    if (nAmp==4 && NChips()==2)
      channel = kPhotometric;
    if (nAmp==4 && NChips()==4)
      channel = kPhotometric + kGuiding;
    for (int chip=0;chip<NChips();chip++) {
      if (fChip[chip]->GetChannel() == kUnknown) {
        fChip[chip]->SetChannel(channel);
      }
    }
  }

  // hack the gains
 // gain hack
  if (channel == kBlueChannel) {
    for (int i=0;i<NChips();i++) {
      fChip[i]->WrDesc("GAIN",DOUBLE,1,&kGainBlue[i]);
    }
  }
  if (channel == kRedChannel) {
    for (int i=0;i<NChips();i++) {
      fChip[i]->WrDesc("GAIN",DOUBLE,1,&kGainRed[i]);
    }
  }
  if (channel & kPhotometric) {
    for (int i=0;i<NChips();i++) {
      fChip[i]->WrDesc("GAIN",DOUBLE,1,&kGainPhot[i]);
    }
  }
}

/* ----- HackGainRatio()  ------------------------------------------- */

void  BiChipSnifs:: HackGainRatio()  {
  
  if (fChip[0]->ParanoMode() || fChip[1]->ParanoMode()) {
    // images shall be bias-substracted for gain correction
    int biasDone0, biasDone1;
    fChip[0]->RdDesc("BIASDONE",INT, 1, &biasDone0);
    fChip[1]->RdDesc("BIASDONE",INT, 1, &biasDone1);
    if (!biasDone0 || !biasDone1 ){
      print_error("BiChipSnifs::HackGainRatio : Bias was not yet substracted\n");
      return;
    }
  }


  float gain;
  fChip[0]->RdDesc("GAIN",FLOAT,1,&gain);
  char key[lg_name+1];
  fChip[0]->RdDesc("DATASEC",CHAR,lg_name+1,key);
  Section* sec = new Section(key);
  sec->SetXFirst(sec->XLast()-30);
  double ratio = GuessGainRatio(sec);

  gain = gain * ratio;
  print_msg("GainRatio is %f\n",ratio);

  fChip[1]->WrDesc("GAIN",FLOAT,1,&gain);

  delete sec;

}
