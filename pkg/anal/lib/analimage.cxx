/* === Doxygen Comment ======================================= */
/*! 
 * \file          analimage.cxx
 * \copyright     (c)  2003 SNIFS-Supernova Factory Experiment
 * \date          Wed Aug  6 18:32:01 2003
 * \author        Emmanuel Gangler <e.gangler@ipnl.in2p3.fr>
 * \version       0.0
 * \brief         
 *                
 * $Id$
 **/
/* =========================================================== */

/* ----- ROOT includes ------------------------------ */
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TH1F.h"

/* ----- image includes ------------------------------ */
#include "imagesnifs.hxx"
#include "analyser.hxx"
#include "section.hxx"

/* ----- local includes ------------------------------ */
#include "analimage.hxx"
#include "roothistos.hxx"
#include "rootutils.hxx"

/* ===== PrintInfoType ============================== */

/* ----- PrintInfoType ----------------------------------- */
template <class T> PrintInfoType<T>::PrintInfoType(const char* Help, const char* Name, const char* Format, int length, T* toprint) {
  strcpy(fHelpInfo,Help);
  strcpy(fName,Name);
  strcpy(fFormat,Format);
  fInfo = toprint;
  // now compute the length
  fLength = length;
  // truncate the name
  if (fLength<=strlen(Name)) {
    strcpy(fName,Name);
    fName[fLength]='\0';
  } else {
    int nblank = fLength - 1 - strlen(Name);
    int i;
    for (i=0;i<nblank;i++) {
      fName[i]=' ';
    }
    fName[i]='\0';
    sprintf(fName+i,"%s ",Name);
  }
}

/* ----- FillHelp ----------------------------------- */
template <class T> void PrintInfoType<T>::FillHelp(char* HelpString){
  sprintf(HelpString,"%s : %s",fName, fHelpInfo);
}

/* ----- FillLine ----------------------------------- */
template <class T> void PrintInfoType<T>::FillLine(char* LineString){
  char *pos = LineString + strlen(LineString);
  unsigned int i;
  for (i =0 ; i < fLength; i ++)
    *(pos+i) = '-';
  *(pos+i)='\0';
}

/* ----- FillLine ----------------------------------- */
template <class T> void PrintInfoType<T>::FillName(char* NameString){
  
  char tmp[ImageSignatureLineLength];
  strcpy(tmp,NameString);
  sprintf(NameString,"%s%s",tmp,fName);
}

/* ----- FillLine ----------------------------------- */
template <class T> void PrintInfoType<T>::FillValue(char* NameString){
  
  char tmp[ImageSignatureLineLength];
  char format[lg_name+1],formatTmp[lg_name+1] = "%s";
  strcpy(tmp,NameString);
  sprintf(format,"%s%s",formatTmp,fFormat);
  sprintf(NameString,format,tmp,*fInfo);
}



/* ===== ImageSignature ============================== */

ClassImp(ImageSignature)

;

/* ----- constants -----------------------------------------*/
const int ImageSignature::fkNitems=18;

/* ----- constructor -----------------------------------------*/
ImageSignature::ImageSignature(){
  Reset();
  fIsAlive = new int[fkNitems];
  SetPrintItems("100000111111110100");
  fPrintInfo = new (PrintInfo*)[fkNitems];
  int i=0;
  fPrintInfo[i++] = new PrintInfoType<char[80]>(" Name of the exposure","Name","%30.30s",30,&fName);
  fPrintInfo[i++] = new PrintInfoType<int>(" FClass ","FClass","%3d",3,&fFClass);
  fPrintInfo[i++] = new PrintInfoType<int>(" External tag","External","%3d",3,&fExternal);
  fPrintInfo[i++] = new PrintInfoType<int>(" 1=B 2=R 4=P","Channel","%2d",2,&fChannel);
  fPrintInfo[i++] = new PrintInfoType<double>(" Exposure time","ExpTime","%4.1f",4,&fExpTime);
  fPrintInfo[i++] = new PrintInfoType<double>(" Julian date","JD","%14.6f",6,&fJulDate);
  fPrintInfo[i++] = new PrintInfoType<double>(" 5.0 Sigma clipped mean","Mean","%6.1f",6,&fRobustMean);
  fPrintInfo[i++] = new PrintInfoType<double>(" 5.0 Sigma clipped RMS","RMS","%6.1f",6,&fRobustRMS);
  fPrintInfo[i++] = new PrintInfoType<int>(" N pixels under 5.0 Sigma","NUnder","%6d",6,&fNUnder);
  fPrintInfo[i++] = new PrintInfoType<int>(" N pixels over 5.0 Sigma","NOver","%6d",6,&fNOver);
  fPrintInfo[i++] = new PrintInfoType<int>(" Number of saturated pixels","Satu","%6d",6,&fSatu);
  fPrintInfo[i++] = new PrintInfoType<double>(" 0.99 quantile (level for continuum exposure)","Q99","%7.1f",7,&fQuant99);
  fPrintInfo[i++] = new PrintInfoType<double>(" 0.999 quantile (level for arc exposure)","Q999","%7.1f",7,&fQuant999);
  fPrintInfo[i++] = new PrintInfoType<double>(" 0.9999 quantile (level for point source)","Q9999","%7.1f",7,&fQuant9999);
  fPrintInfo[i++] = new PrintInfoType<double>(" magic focus info","Focus","%4.2f",4,&fFocus);
  fPrintInfo[i++] = new PrintInfoType<double>(" RMS of the noise","Noise","%5.2f",5,&fNoise);
  fPrintInfo[i++] = new PrintInfoType<double>(" Ratio of the median (gain stability)","RMed","%5.2f",5,&fGainMed);
  fPrintInfo[i++] = new PrintInfoType<double>(" Ratio of the 0.995 quantiles (gain stability)","RQuant","%4.2f",4,&fGainQuant);
  
}

/* ----- destructor -----------------------------------------*/
ImageSignature::~ImageSignature(){

  for (int i=0;i<fkNitems;i++)
    delete fPrintInfo[i];
  delete[] fPrintInfo;
}


/* ----- Reset -----------------------------------------*/
void ImageSignature::Reset() {
  fExternal = 0;
  fChannel=0;
  fFClass=0;
  fExpTime=0;
  fJulDate=0;
  fSatu =0;
  fQuant99 = 0;
  fQuant999= 0;
  fFocus=0;
  fNoise=0;
  fGainMed=0;
  fGainQuant=0;
  fName[0]='\0';
  
}

/* ----- Reset -----------------------------------------*/
void ImageSignature::SetPrintItems(char* ItemsDescr) {

  for (int i=0;i<fkNitems;i++) {
    if (ItemsDescr[i]=='1')
      fIsAlive[i] = 1;
    else
      fIsAlive[i] = 0;
  }
}

/* ----- FillExternal -----------------------------------------*/
void ImageSignature::FillExternal(int External) {
  fExternal = External;
}

/* ----- FillExternal -----------------------------------------*/
void ImageSignature::Fill(ImageSnifs* Preprocessed) {

  fFClass = Preprocessed->GetFClass();
  fChannel = Preprocessed->GetChannel();

  double expTime;
  Preprocessed->RdDesc("EXPTIME",DOUBLE,1,&expTime);
  fExpTime=expTime;

  // JulianDate
  double jd;
  if (Preprocessed->RdIfDesc("JD",DOUBLE,1,&jd)<=0)
    jd=0;
  fJulDate=jd;

  //
  // Bad (saturating) pixels
  //
  int nAmp;
  if (Preprocessed->RdIfDesc("CCDNAMP",INT,1,&nAmp)<=0)
    nAmp = 2;
  int saturate;
  Preprocessed->RdDesc("SATURATE",INT,1,&saturate);
  // useful definitions for the rest of the method
  Section * dataSec = Preprocessed->DataSec();
  dataSec->SetName("Data");
  Section * subSec[2] ;
  RootAnalyser anaSub[2];
  fSatu =0;
  // analyse only 2 first amps ! 
  for (int iamp =0;iamp<2;iamp++){
    double gain;
    char gainKey[lg_name+1], secName[lg_name+1];
    sprintf(gainKey,"CCD%dGAIN",iamp);
    Preprocessed->RdDesc(gainKey,DOUBLE,1,&gain);

    sprintf(secName,"SubSec%d",iamp);
    int offset = dataSec->XLength()/nAmp;
    subSec[iamp] = new Section(dataSec->X1()+offset*iamp,dataSec->X1()-1+offset*(iamp+1),dataSec->Y1(),dataSec->Y2(),secName);
    anaSub[iamp].SetImage(Preprocessed);
    anaSub[iamp].SetSection(subSec[iamp]);
    // hack : some quantity was already removed from the saturation (overscan)
    fSatu += anaSub[iamp].NPixOver(saturate *0.9 * gain);
  }
  
  // do not know how to do and if it is relevant ...
  //  fOverscanLevel=8;

  RootAnalyser ana(Preprocessed,dataSec);

  TH1F* hData = ana.HistoDataBuild(dataSec->Name(),10000);
  RootUtils::SigmaClip(hData,5.0);
  fRobustMean = hData->GetMean();
  fRobustRMS = hData->GetRMS();
  fNUnder = (int) hData->Integral(0,hData->GetXaxis()->GetFirst()-1);
  fNOver = (int) hData->Integral(hData->GetXaxis()->GetLast()+1,hData->GetNbinsX()+1);
  
  // very slow
  //ana.SigmaClippedInfo(5.0,&fRobustMean, &fRobustRMS, &fNOutliers);
  
  // the quantiles
  double pbSum = 0.99;
  hData->GetQuantiles(1,&fQuant99,&pbSum);
  pbSum=0.999;
  hData->GetQuantiles(1,&fQuant999,&pbSum);
  pbSum=0.9999;
  hData->GetQuantiles(1,&fQuant9999,&pbSum);
  
  fFocus=0;

  // the noise
  double rdnoise[nAmp];
  Preprocessed->RdDesc("RDNOISE",DOUBLE,nAmp,rdnoise);
  fNoise=sqrt(rdnoise[0]*rdnoise[0]/2+rdnoise[1]*rdnoise[1]/2);

  // gain
  // We reduce here the analysis section (quantile is an expensive routine)
  // Well ... not necessary any more
  //
  // if (subSec[0]->XLength()>500) 
  //  subSec[0]->SetX1(subSec[0]->X2() - 500);
  //if (subSec[1]->XLength()>500) 
  //  subSec[1]->SetX2(subSec[1]->X1() + 500);
  //for (int i=0;i<2;i++) {
  //  if (subSec[i]->YLength()>1000){
  //    int yMiddle = subSec[i]->YLength()/2+subSec[i]->Y1();
  //    subSec[i]->SetY1(yMiddle-499);
  //    subSec[i]->SetY2(yMiddle+500);
  //  }
  //  anaSub[i].SetSection(subSec[i]);
  //}
  TH1F* hSub[2];
  double q5[2],q999[2];
  for (int i=0;i<2;i++) {
    hSub[i] = ana.HistoDataBuild(subSec[i]->Name(),10000);
    pbSum = 0.5;
    hSub[i]->GetQuantiles(1,q5+i,&pbSum);
    pbSum = 0.999;
    hSub[i]->GetQuantiles(1,q999+i,&pbSum);
  }
  
  
  fGainMed=q5[0]/q5[1];
  fGainQuant=q999[0]/q999[1];

  strcpy(fName,Preprocessed->Name());
  
  // cleanup
  delete subSec[0];
  delete subSec[1];
  delete dataSec;
  delete hSub[0];
  delete hSub[1];
  delete hData;
}

/* ----- PrintHeader  -----------------------------------------*/
void ImageSignature::PrintBlurb() {

  for (int i=0;i<fkNitems;i++) {
    char blurb[lg_name+1];
    if (fIsAlive[i]) {
      fPrintInfo[i]->FillHelp(blurb);
      print_msg(blurb);
    }
  }
}

/* ----- PrintHeader  -----------------------------------------*/
void ImageSignature::PrintHeader() 
{

  PrintTrailer();
  char Names[ImageSignatureLineLength];
  Names[0]='|';
  Names[1]='\0';
  for (int i=0;i<fkNitems;i++)
    if (fIsAlive[i]) {
      fPrintInfo[i]->FillName(Names);
      size_t len=strlen(Names);
      Names[len]='|';
      Names[len+1]='\0';
    }
  print_msg(Names);
  PrintTrailer();
}


/* ----- FillImageSignature -----------------------------------------*/
void ImageSignature::PrintContent() {

  char Values[ImageSignatureLineLength];
  Values[0]='|';
  Values[1]='\0';
  for (int i=0;i<fkNitems;i++)
    if (fIsAlive[i]) {
      fPrintInfo[i]->FillValue(Values);
      int len=strlen(Values);
      Values[len]='|';
      Values[len+1]='\0';
    }
  print_msg(Values);
}

/* ----- PrintTrailer  -----------------------------------------*/
void ImageSignature::PrintTrailer(){

  char Line[ImageSignatureLineLength];
  Line[0]='+';
  Line[1]='\0';
  for (int i=0;i<fkNitems;i++)
    if (fIsAlive[i]) {
      fPrintInfo[i]->FillLine(Line);
      size_t len = strlen(Line);
      Line[len]='+';
      Line[len+1]='\0'; // the strlen is not defined at this stage
    }
  print_msg(Line);
}

/* ===== AnalImageSignature ============================== */

/* ----- constructor -----------------------------------------*/
AnalImageSignature::AnalImageSignature(TFile * Rfile)
{
  fImage=0;
  fSignature=new ImageSignature();
  fFile = Rfile;
  fTree=0;
  MakeTree();
}

/* ----- destructor -----------------------------------------*/
AnalImageSignature::~AnalImageSignature()
{
  if (fFile) {
    fFile->Write();
    fFile->Purge();
    fFile->Close();
  }
  
    
  delete fSignature;
}

/* ----- SetImage -----------------------------------------*/
void AnalImageSignature::MakeTree() 
{
  if (fFile) {
    fTree = (TTree*) fFile->Get("Signature");
    if (!fTree) {
      fTree = new TTree("Signature","ImageSignatures");
    }
    TBranch * branch;
    char branchName[]={"AllImages"};
    if ( (branch = fTree->GetBranch(branchName)) )
      branch->SetAddress(fSignature);
    else
      branch = fTree->Branch(branchName,"ImageSignature",&fSignature);
  }
}

/* ----- StoreImageSignature -----------------------------------------*/
void AnalImageSignature::StoreImageSignature() 
{
  if (fTree)
    fTree->Fill();
  fSignature->Reset();
}

