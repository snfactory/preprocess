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
