/* === Doxygen Comment ======================================= */
/*! 
 * \file          kombinatorfit.cxx
 * \copyright     (c) 2003 IPNL
 * \date          Wed Aug  6 17:44:35 2003
 * \author        Emmanuel Gangler <e.gangler@ipnl.in2p3.fr>
 * \version       0.0
 * \brief         
 *                
 * $Id$
 **/
/* =========================================================== */

#include <math.h>
//#include "IFU_io.h"
#include "utils.h" // ut_big_value
#include "kombinatorfit.hxx"

/* #####  KFGain ################################################# */


/* ===== constructor/destructor ======================================= */

/* ----- KFGain  -------------------------------------------------- */
KFGain::KFGain(double MaxSigmaOutlier) {
  fSigma = MaxSigmaOutlier;
}

/* ===== method ======================================= */

/* ----- Kombine  -------------------------------------------------- */
void KFGain::KombineFit(vector<double> * X,vector<double> * Vals, vector<double> * Vars, double * Val, double * Var) {
  
  fX = X;
  fVals = Vals;
  fVars.assign(fVals->size(),1);
  fVal = Mean();
  // throw-away
  while (RemoveFar()) {
    fVal = Mean();
  }

  *Val = fVal;
}

/* ----- Mean  -------------------------------------------------- */
double KFGain::Mean() {

  double sum=0;
  double sumw=0;
  // First compute the light parameter
  for (unsigned int i=0;i<fVals->size(); i++) {
    sum += (*fVals)[i]/(*fX)[i]/fVars[i] ;
    sumw += 1/fVars[i] ;
  }
  fSlope = sum/sumw;
  // compute now the chi2
  sum=0;
  for (unsigned int i=0;i<fVals->size(); i++) {
    sum += ((*fVals)[i]-fSlope*(*fX)[i])*((*fVals)[i]-fSlope*(*fX)[i])
      / fSlope / (*fX)[i] / fVars[i];
  }
  if (sumw>1)
    fVal = sum/(sumw-1);
  else
    fVal = ut_big_value;
  return fVal;
}

/* ----- RemoveFar  -------------------------------------------------- */
int KFGain::RemoveFar() {

  if (fSigma<=0 || fVals->size()<=2 )
    return 0;

  int iRm=-1;
  double MaxSigma=0;
  double oldVal = fVal;
  double oldSlope = fSlope;

  for ( unsigned int i=0;i<fVals->size();i++) {
    // excluded already
    if (fVars[i]>ut_big_value*0.99)
      continue;
    
    // compute values with current value substracted
    double oldVar = fVars[i];
    fVars[i] = ut_big_value;
    Mean();
    fVars[i] = oldVar;

    double sigma = fabs((*fVals)[i] - fSlope*(*fX)[i])/
                   sqrt( fSlope * (*fX)[i] * fVal);
    if (sigma > MaxSigma && sigma > fSigma) {
      iRm = i;
      MaxSigma = sigma;
    }
  }
  fVal = oldVal;
  fSlope = oldSlope;

  if (iRm>=0) {
    fVars[iRm] = ut_big_value;
    return 1;
  }
  return 0;
}

