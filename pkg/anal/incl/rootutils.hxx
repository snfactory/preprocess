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

#ifndef ROOTUTILS_H
#define ROOTUTILS_H

/* ----- ROOT includes ------------------------------ */
class TH1;

/* ===== ROOT HISTOS ============================== */

class RootUtils {
  public :

  static void SigmaClip(TH1 * h, double sigCut);
  
  
};

#endif

