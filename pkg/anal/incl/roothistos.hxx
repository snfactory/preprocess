/* === Doxygen Comment ======================================= */
/*! 
 * \file          analyser.hxx
 * \copyright     (c)  2003 SNIFS-Supernova Factory Experiment
 * \date          Wed Aug  6 18:32:01 2003
 * \author        Emmanuel Gangler <e.gangler@ipnl.in2p3.fr>
 * \version       0.0
 * \brief         
 *                
 * $Id$
 **/
/* =========================================================== */

#ifndef ROOTHISTOS_H
#define ROOTHISTOS_H

/* ----- ROOT includes ------------------------------ */
class TH1F;
class TH1;

/* ----- image includes ------------------------------ */
class ImageSnifs;
class ImageSimple;
class Section;
#include "analyser.hxx"

/* ===== ROOT HISTOS ============================== */

class RootAnalyser : public ImageAnalyser {
  public :
  
  RootAnalyser() {}
  RootAnalyser(ImageSimple * Image, Section *Sec);
  //  RootAnalyser(ImageSnifs * Image, Section *Sec);
  ~RootAnalyser(){}
 
  /* ----- Tools ------------------------------------------ */
  void HistoSetMinMax(TH1 * histo);

  /* ----- Analysis ---------------------------------------- */

  void FillHistoLine();

  void HorizontalProfile(float sigma=0);

  void VerticalProfile(float sigma=0);

  void OddEvenVerticalProfile(float sigma=0);

  void OverscanError(float sigma=0);

  void HistoData(int nbins=4000);
  TH1F* HistoDataBuild(char * histName ,int nbins=4000);

  void MatrixData();

  void HighFrequency();

  void VertHighFrequencyProf(double nmin=0,double nmax=0, int nbin=1000);

  void Fft();

  void ADCBits(int nBits=16);
  
};

#endif

