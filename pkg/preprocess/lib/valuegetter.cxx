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

/* ##### ValuesGetterDarkFitter ################################################# */

/* ===== constructor/destructor ======================================= */

/* ----- constructor  ----------------------------------*/
ValuesGetterDarkFitter::ValuesGetterDarkFitter(DarkModel* Model,int * activate){
  if (Model->GetNsec() != 1)
    print_error("ValuesGetterDarkFitter : needs a 1-section model");
  fDarkModel = Model;
  fNParams=0;
  for (i=0;i<3;i++) {
    if (activate && !activate[i]) 
      fActive[i]=0;      
    else {
      fActive[i]=1;
      fNParams+=1;
    }
  }
}

/* ===== methods ======================================= */

/* ----- GetValue ----------------------------------*/
void ValuesGetterDarkFitter::GetValues(ImageSimple* Image, gsl_vector* retValues){
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
    print_error("%s has a bad time on",image->Name());
  }

  int count=0;
  if (fActive[0]) {
    gsl_vector_set(retValues,count,fDarkModel->GetI0(0)*texp);
    count++;
  }
  if (fActive[1]) {
    gsl_vector_set(retValues,count,fDarkModel->GetI1(0)*fDarkModel->DarkTimeTerm(timeon,texp,0));
    count++;
  }
  if (fActive[2]) {
    gsl_vector_set(retValues,count,fDarkModel->GetI0(0)*texp);
    count++;
  }

}
