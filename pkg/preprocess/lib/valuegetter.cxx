/* === Doxygen Comment ======================================= */
/*! 
 * \file          valuegetter.cxx
 * \copyright     (c) 2003 SNIFS-Supernova Factory Experiment
 * \date          Wed Aug  6 17:44:35 2003
 * \author        Emmanuel Gangler <e.gangler@ipnl.in2p3.fr>
 * \version       0.0
 * \brief         
 *                
 * $Id$
 **/
/* =========================================================== */

//#include "section.hxx"
//#include "image.hxx"
#include "analyser.hxx"
#include "valuegetter.hxx"
#include "darkmodel.hxx"

/* ##### MeanFromAnalyser ################################################# */

/* ===== constructor/destructor ======================================= */

/* ----- constructor  ----------------------------------*/
ValueAnalyserMean::ValueAnalyserMean(Section* Sec){
  fAnalyser = new ImageAnalyser();
  fAnalyser->SetSection(Sec);
}

/* ----- destructor  ----------------------------------*/
ValueAnalyserMean::~ValueAnalyserMean(){
  delete fAnalyser;
}

/* ===== methods ======================================= */

/* ----- GetValue ----------------------------------*/
double ValueAnalyserMean::GetValue(ImageSimple* Image){
  fAnalyser->SetImage(Image);
  return fAnalyser->MeanLevel();
}

/* ##### Value for the dark current ################################################# */

/* ===== constructor/destructor ======================================= */

/* ----- constructor  ----------------------------------*/
ValueDark::ValueDark(DarkModel* darkModel){
  fDarkModel = darkModel;
}

/* ----- destructor  ----------------------------------*/
ValueDark::~ValueDark(){
}

/* ===== methods ======================================= */

/* ----- GetValue ----------------------------------*/
double ValueDark::GetValue(ImageSimple* Image){
  // returns the dark current estimation
  ImageSnifs* image = (ImageSnifs*) Image;

  double timeon,temp,texp;

  image->RdDesc("DETTEMP",DOUBLE,1,&temp); 
  image->RdDesc("DARKTIME",DOUBLE,1,&texp); 

  char timeOnStr[lg_name+1];
  image->RdDesc("TIMEON",CHAR,lg_name+1,timeOnStr); // CAVEAT : TIMEON not defined for all data
  if (strcmp(timeOnStr,"None"))
    image->RdDesc("TIMEON",DOUBLE,1,&timeon); // CAVEAT : TIMEON not defined for all data.
  else
    timeon=-1;
  if (timeon<texp) {
    print_warning("%s has a bad time on",image->Name());
    timeon=-1;
  }

  //Section ** Secs=fDarkModel->GetSections();
  double toremove=0;
  for (int isec=0;isec<fDarkModel->GetNsec();isec++){
    if (timeon>0)
      toremove+=fDarkModel->DarkSub(temp,timeon,texp,isec);
    else
      toremove+=(fDarkModel->GetI0(isec)+fDarkModel->GetI2(isec)*fDarkModel->TempTerm(temp))*texp;
  toremove /= fDarkModel->GetNsec();

  return toremove;
  }
}
