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



/* ----- local includes ------------------------------ */
class PrintInfo;
#include "signaturecut.hxx"


/* ===== SignatureCut ============================== */

/* ----- SignatureCut ----------------------------------- */
SignatureCut::SignatureCut(int Fclass, int Channel,  PrintInfo* I, int Cut, double Val, double OptVal, char* Message )  {

  fFclass = Fclass;
  fChannel = Channel;
  fInfo = I;
  fCut = Cut;
  fValue = Val;
  fOptValue = OptVal;
  strcpy(fMessage, Message);
  
  }
  
