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

class ImageSnifs;
class ImageSimple;
class Section;
class TH1;
#include "analyser.hxx"

/* ===== ROOT HISTOS ============================== */

class RootAnalyser : public ImageAnalyser {
  public :
  
  RootAnalyser() {}
  RootAnalyser(ImageSimple * Image, Section *Sec);
  RootAnalyser(ImageSnifs * Image, Section *Sec);
  ~RootAnalyser(){}
 
  /* ----- Tools ------------------------------------------ */
  void HistoSetMinMax(TH1 * histo);

  /* ----- Analysis ---------------------------------------- */

  void FillHistoLine();

  void HorizontalProfile(float sigma=0);

  void VerticalProfile(float sigma=0);

  void OddEvenVerticalProfile(float sigma=0);

  void OverscanError(float sigma=0);

  void HistoData();

  void MatrixData();

  void HighFrequency();

  void Fft();

  void ADCBits(int nBits=16);
  
};

#endif

