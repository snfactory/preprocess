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

#include "kombinatorfit.hxx"

#include <math.h>
//#include "IFU_io.h"
#include "utils.h" // ut_big_value
#define HAVE_INLINE
#include "gsl/gsl_vector.h"
#include "gsl/gsl_matrix.h"


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

/* #####  KFLinear ################################################# */


/* ===== constructor/destructor ======================================= */

/* ----- KFLinear  -------------------------------------------------- */
KFLinear::KFLinear(double MaxSigmaOutlier, int NParams, int AutoRescale) {
  fSigma = MaxSigmaOutlier;
  fNParams = NParams;
  fVal = gsl_vector_calloc(fNParams);
  fCovar = gsl_matrix_calloc(fNParams,fNParams);
  fAutoRescale=AutoRescale;
}

/* ----- KFLinear  -------------------------------------------------- */
KFLinear::~KFLinear() {
  gsl_vector_free(fVal);
  gsl_matrix_free(fCovar);
}


/* ===== method ======================================= */

/* ----- Kombine  -------------------------------------------------- */
void KFLinear::KombineFit(gsl_matrix * X,gsl_vector * Vals, gsl_vector * Vars, gsl_vector * Val, gsl_matrix * Covar) {
  
  fX = X; // X should be a matrix fVals->size x fNParams
  fVals = Vals;
  fWeights = gsl_vector_calloc(fVals->size);
  if (fAutoRescale)
    for (int i=0;i<fNParams; i++)
       gsl_vector_set(fWeights,i,1.);
  else 
    for (int i=0;i<fNParams; i++)
       gsl_vector_set(fWeights,i,1./gsl_vector_get(Vars,i));
	 
  fWspace = gsl_multifit_linear_alloc(fVals->size,fNParams);

  fChi2=Fit();
  // throw-away
  while (RemoveFar()) {
    fChi2=Fit();
  }
  
	   
  gsl_vector_memcpy(Val,fVal);
  gsl_matrix_memcpy(Covar,fCovar);
  gsl_vector_free(fWeights);
  gsl_multifit_linear_free(fWspace);
}

/* ----- Fit  -------------------------------------------------- */
double KFLinear::Fit() {
  // Fitting Y = Sum ( a_i x_i)
  double chi2;
  int ndf=fVals->size - fVal->size;

  gsl_multifit_wlinear(fX,fWeights,fVals,fVal,fCovar,&chi2,fWspace);
  if (fAutoRescale) {
    gsl_matrix_scale(fCovar,chi2/ndf);
  }


  return chi2;
}

/* ----- RemoveFar  -------------------------------------------------- */
int KFLinear::RemoveFar() {

  if (fSigma<=0 || fVals->size <= (unsigned int)fNParams+1 )
    return 0;

  int nRm=-1;
  double MaxSigma=0;

  for ( unsigned int n=0;n<fVals->size;n++) {
    // excluded already
    if (gsl_vector_get(fWeights,n)==0.0)
      continue;
    
    // compute values with current value substracted
    double wsave=gsl_vector_get(fWeights,n);
    gsl_vector_set(fWeights,n,0.0);
    double chi2=Fit();
    double scale=1.;
    if (fAutoRescale)
      scale=chi2/(fVals->size - fVal->size);
    gsl_vector_set(fWeights,n,wsave);

    // Get the model value
    double model,modelerr;
    gsl_vector_view x=gsl_matrix_row(fX,n);
    gsl_multifit_linear_est(&x.vector,fVal,fCovar,&model,&modelerr);

    double sigma = fabs(gsl_vector_get(fVals,n) - model)/sqrt( modelerr*modelerr + scale/gsl_vector_get(fWeights,n));
    if (sigma > MaxSigma && sigma > fSigma) {
      nRm = n;
      MaxSigma = sigma;
    }
  }

  if (nRm>=0) {
    gsl_vector_set(fWeights,nRm,0.0);
    return 1;
  }
  return 0;
}

