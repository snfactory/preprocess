/* === Doxygen Comment ======================================= */
/*! 
 * \file          kombinator.cxx
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
#include "IFU_io.h"
#include "utils.h"
#include "kombinator.hxx"

/* #####  KSpread ################################################# */


/* ===== constructor/destructor ======================================= */

/* ----- KSpread  -------------------------------------------------- */
KSpread::KSpread(double MaxSigmaOutlier,int SpreadMode) {
  fSigma = MaxSigmaOutlier;
  fSpread = SpreadMode;
}

/* ===== method ======================================= */

/* ----- Kombine  -------------------------------------------------- */
void KSpread::Kombine(vector<double> * Vals, vector<double> * Vars, double * Val, double * Var) {
  /* The Kombine assumes the Vals are calibrated in electrons, 
     and the Vars contain gaussian + poisson errors */
  
  fVals = Vals;
  fVars.assign(fVals->size(),1);
  fVal = Mean();
  // throw-away
  while (RemoveFar(fSigma)) {
    fVal = Mean();
  }

  *Val = fVal;
  if (!fSpread) {
    *Var = fVar/fSum;
  } else {
    *Var = fVar;
  }
}

/* ----- Mean  -------------------------------------------------- */
double KSpread::Mean() {
  // fits the fVal given the 
  double sum=0;
  double sumw=0;
  //double weight;
  for (unsigned int i=0;i<fVals->size(); i++) {
    sum += (*fVals)[i]/fVars[i] ;
    sumw += 1/fVars[i] ;
  }
  fVal = sum/sumw;
  fSum = sumw;
  // compute now the estimated variance
  sum=0;
  for (unsigned int i=0;i<fVals->size(); i++) {
    sum += ((*fVals)[i]-fVal)*((*fVals)[i]-fVal)/fVars[i] ;
  }
  if (sumw>1)
    fVar = sum/(sumw-1);
  else
    fVar = ut_big_value;
  return fVal;
}

/* ----- RemoveFar  -------------------------------------------------- */
int KSpread::RemoveFar(double SigmaCut) {

  if (fSigma<=0 || fVals->size()<=2 )
    return 0;

  int iRm=-1;
  double MaxSigma=0;

  for ( unsigned int i=0;i<fVals->size();i++) {
    // excluded already
    if (fVars[i]>ut_big_value*0.99)
      continue;
    
    // compute values with current value substracted
    double newSum = fSum - 1/fVars[i];
    double newVal = (fVal * fSum - (*fVals)[i]/fVars[i]) / newSum;
    double newVar = 
      ( fVar * (fSum-1) 
        + fSum * (fVal-newVal) * (fVal-newVal)
        -  ((*fVals)[i]-newVal)*((*fVals)[i]-newVal)/fVars[i] ) 
      / (newSum -1);
    
    double sigma = fabs((*fVals)[i] - newVal)/sqrt(newVar);
    if (sigma > MaxSigma && sigma > SigmaCut) {
      iRm = i;
      MaxSigma = sigma;
    }
  }
  if (iRm>=0) {
    fVars[iRm] = ut_big_value;
    return 1;
  }
  return 0;
}

/* #####  KGauss ################################################# */


/* ===== constructor/destructor ======================================= */

/* ----- KGaussPoisson  -------------------------------------------------- */
KGauss::KGauss(double MaxSigmaOutlier) {
  fSigma = MaxSigmaOutlier;
  fBrute=0;
}

/* ===== method ======================================= */

/* ----- Kombine  -------------------------------------------------- */
void KGauss::Kombine(vector<double> * Vals, vector<double> * Vars, double * Val, double * Var) {
  /* The Kombine assumes the Vals are calibrated in electrons, 
     and the Vars contain gaussian + poisson errors */
  
  fVals = Vals;
  fVars = Vars;
  if (fVals->size() != fVars->size()) {
      print_error("KGauss::Kombine called with different size vectors\n");
      return;
    }

  fVal = WeightedMean();
  // throw-away
  while (RemoveFar(fSigma)) {
    fVal = WeightedMean();
  }

  *Val = fVal;
  *Var = fVar;
}

/* ----- MakeAGuess  -------------------------------------------------- */
double KGauss::WeightedMean() {
  // fits the fVal given the 
  double sum=0;
  double sumw=0;
  //double weight;
  for (unsigned int i=0;i<fVals->size(); i++) {
    sum += (*fVals)[i] / (*fVars)[i];
    sumw += 1 / (*fVars)[i];
  }
  fVar = 1/sumw;
  return sum/sumw;
}

/* ----- RemoveFar  -------------------------------------------------- */
int KGauss::RemoveFar(double SigmaCut) {

  if (fSigma<=0)
    return 0;

  int iRm=-1;
  double MaxSigma=0;

  for ( unsigned int i=0;i<fVals->size();i++) {
    //we compute it the current point excluded;
    double sigma = fabs((*fVals)[i] - fVal)/sqrt((*fVars)[i]-fVar);
    if (sigma > MaxSigma && sigma > SigmaCut) {
      iRm = i;
      MaxSigma = sigma;
      }      
    if (fBrute && sigma>SigmaCut) {
      (*fVars)[i] = ut_big_value;
    }
  }
  if (iRm>=0) {
    (*fVars)[iRm] = ut_big_value;
    return 1;
  }
  return 0;
}

/* #####  KGaussPoisson ################################################# */


/* ===== constructor/destructor ======================================= */

/* ----- KGaussPoisson  -------------------------------------------------- */
KGaussPoisson::KGaussPoisson(double MaxSigmaOutlier) {
  fSigma = MaxSigmaOutlier;
  fPrecision = 0.005;
  fMaxIter = 10;
}

/* ===== method ======================================= */

/* ----- Kombine  -------------------------------------------------- */
void KGaussPoisson::Kombine(vector<double> * Vals, vector<double> * Vars, double * Val, double * Var) {
  /* The Kombine assumes the Vals are calibrated in electrons, 
     and the Vars contain gaussian + poisson errors */
  
  fVals = Vals;
  fVars = Vars;
  if (fVals->size() != fVars->size()) {
      print_error("KGaussPoisson::Kombine called with different size vectors\n");
      return;
    }
  if (fVals->size() != fWVars.size())
    fWVars.assign(fVals->size(),0);

  fVal = WeightedMean();
  // quick throw-away
  while (RemoveFar(2*fSigma) ) {
    InitWVars();
    fVal = WeightedMean();
  }
  for (;RemoveFar(fSigma);
       InitWVars(),fVal = WeightedMean()) {
    // converge
    int iter=0;
    double oldVal;
    do {
      oldVal = fVal;
      UpdateWVars(fVal);
      fVal = WeightedMean();
    } while(fabs(fVal-oldVal) > sqrt(fVar)*fPrecision && ++iter < fMaxIter);
  }

  *Val = fVal;
  *Var = fVar*(1+fPrecision)*(1+fPrecision);
}

/* ----- WeightedMean  -------------------------------------------------- */
double KGaussPoisson::WeightedMean() {
  // only one iteration
  vector<double>::const_iterator ValCI, VarCI;
  double sum=0;
  double sumw=0;
  for (ValCI=fVals->begin(), VarCI=fWVars.begin();
       ValCI!=fVals->end() && VarCI!=fWVars.end();
       ++ValCI,++VarCI) {
    sum += *ValCI / *VarCI;
    sumw += 1 / *VarCI;
  }
  fVar = 1/sumw;
  return sum/sumw;
}

/* ----- FirstWVars  -------------------------------------------------- */
void KGaussPoisson::InitWVars() {
  // reset of working variance
 
  for ( unsigned int i=0;i<fVals->size();i++)
    fWVars[i] = (*fVars)[i];
}

/* ----- UpdateWVars  -------------------------------------------------- */
void KGaussPoisson::UpdateWVars(double Val) {

  for ( unsigned int i=0;i<fVals->size();i++ ) {
    fWVars[i] = (*fVars)[i] ;
    if ((*fVals)[i]>0)
      fWVars[i] -= (*fVals)[i];
    if (Val > 0)
      fWVars[i] += Val;
  }
}


/* ----- RemoveFar  -------------------------------------------------- */
int KGaussPoisson::RemoveFar(double SigmaCut) {

  if (fSigma<=0)
    return 0;

  int iRm=-1;
  double MaxSigma=0;

  for ( unsigned int i=0;i<fVals->size();i++) {
  double sigma = fabs((*fVals)[i] - fVal)/sqrt(fWVars[i]);
    if (sigma > MaxSigma && sigma > SigmaCut) {
      iRm = i;
      MaxSigma = sigma;
    }
  }
  if (iRm>=0) {
    (*fVars)[iRm] = ut_big_value;
    fWVars[iRm] = ut_big_value;
    return 1;
  }
  return 0;
}


/* #####  KMedian ################################################# */


/* ===== constructor/destructor ======================================= */

/* ----- KMedian  -------------------------------------------------- */
KMedian::KMedian() {
}

/* ===== method ======================================= */

/* ----- Kombine  -------------------------------------------------- */
void KMedian::Kombine(vector<double> * Vals, vector<double> * Vars, double * Val, double * Var) {
  /* The Kombine assumes the Vals are calibrated in electrons, 
     and the Vars contain errors */
  
  double *vals = new double[Vals->size()];
  double sumw=0;
  for (unsigned int i=0;i<Vals->size(); i++) {
    vals[i]=(*Vals)[i];
    sumw += 1/(*Vars)[i];
  }
  
  *Val = ut_median(vals,Vals->size());
  // rough estimate of variance
  *Var = 1./sumw*1.571;
  delete[] vals;
}
  

