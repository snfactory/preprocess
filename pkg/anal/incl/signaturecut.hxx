/* === Doxygen Comment ======================================= */
/*! 
 * \file          signaturecut.hxx
 * \copyright     (c)  2003 SNIFS-Supernova Factory Experiment
 * \date          Wed Aug  6 18:32:01 2003
 * \author        Emmanuel Gangler <e.gangler@ipnl.in2p3.fr>
 * \version       0.0
 * \brief         
 *                
 * $Id$
 **/
/* =========================================================== */

#ifndef SIGNATURECUT_H
#define SIGNATURECUT_H

/* ----- preprocess includes ------------------------------ */
#include "utils.h"

/* ----- local includes ------------------------------ */
#include "analimage.hxx"

/* ===== SignatureCut ============================== */

class SignatureCut {
public:
  SignatureCut(int Fclass, int Channels, PrintInfo* I, int Cut, double Val, double OptVal, char* Message);

  int GetFclass() {return fFclass;}
  int GetChannel() {return fChannel;}
  double GetVariable() {return fInfo->GetValue();}
  int GetCut() {return fCut;}  
  double GetValue() {return fValue;}
  double GetOptValue() {return fOptValue;}
  char* GetMessage() {return fMessage;}
  PrintInfo* Info() {return fInfo;}
  
  

  
protected:
  int fFclass;
  int fChannel;
  PrintInfo* fInfo;
  int fCut;
  double fValue;
  double fOptValue;
  char fMessage[lg_message];
  
};

#endif
