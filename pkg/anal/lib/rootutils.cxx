/* === Doxygen Comment ======================================= */
/*! 
 * \file          rootutils.hxx
 * \copyright     (c)  2003 SNIFS-Supernova Factory Experiment
 * \date          Wed Aug  6 18:32:01 2003
 * \author        Emmanuel Gangler <e.gangler@ipnl.in2p3.fr>
 * \version       0.0
 * \brief         
 *                
 * $Id$
 **/
/* =========================================================== */

/* ----- ROOT includes ----- */
#include "TH1.h"

/* ----- local includes ----- */
#include "rootutils.hxx"

/* --------------- HistoSetMinMax ------------------- */
void RootUtils::SigmaClip(TH1 * hist, double sigCut) {

  // reset all range info
  hist->GetXaxis()->SetRange();
  
  double integral, oldIntegral=0;
  while ((integral = hist->Integral())!=oldIntegral) {

    double mean = hist->GetMean();
    double rms = hist->GetRMS();
    hist->GetXaxis()->SetRangeUser(mean-sigCut*rms,mean+sigCut*rms);
    oldIntegral = integral;
  }
  

}

