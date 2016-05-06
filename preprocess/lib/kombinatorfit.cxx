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
//#define HAVE_INLINE // Seems to be buggy with gsl 1.14
#include "gsl/gsl_sort_vector.h"



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
KFLinear::KFLinear(double MaxSigmaOutlier, int NParams, int AutoRescale, int WeightsOnlyForSelection,int Dark) {
  fSigma = MaxSigmaOutlier;
  fNParams = NParams;
  fVal = gsl_vector_calloc(fNParams);
  fCovar = gsl_matrix_calloc(fNParams,fNParams);
  if (WeightsOnlyForSelection)
    fWeightsOnlyForSelection=fAutoRescale=WeightsOnlyForSelection; // with this option, weights are used only to reject outliers, not to perform the fit
  else
    fAutoRescale=AutoRescale;
  fDark=Dark; // special flag for the dark reduciton problem

  fWeights=0;
  fTmpSort=0;
  fVarX=0;
  fWspace = 0;

}

/* ----- KFLinear  -------------------------------------------------- */
KFLinear::~KFLinear() {
  WorkspaceFree();
  gsl_vector_free(fVal);
  gsl_matrix_free(fCovar);
}

/* ----- AllocWorkspace  -------------------------------------------------- */
void KFLinear::WorkspaceAlloc(int N) {
  if (fWeights != 0 )
    WorkspaceFree();
  fWeights = gsl_vector_calloc(N);
  fTmpSort = gsl_vector_calloc(N);
  fVarX = gsl_vector_alloc(N);
  fWspace = gsl_multifit_linear_alloc(N,fNParams);
}

/* ----- FreeWorkspace  -------------------------------------------------- */
void KFLinear::WorkspaceFree() {
  if (fWeights == 0)
    return;
  gsl_vector_free(fWeights);
  fWeights=0;
  gsl_vector_free(fTmpSort);
  fTmpSort=0;
  gsl_vector_free(fVarX);
  fVarX=0;
  gsl_multifit_linear_free(fWspace);
  fWspace = 0;
}

/* ===== method ======================================= */

/* ----- Kombine  -------------------------------------------------- */
void KFLinear::KombineFit(gsl_matrix * X,gsl_vector * Vals, gsl_vector * Vars, gsl_vector * Val, gsl_matrix * Covar) {
  
  fX = X; // X should be a matrix fVals->size x fNParams
  fVals = Vals;
  fVars = Vars;

  if (fDark) {
    // First compute the median of predicted errors, used to rescale the information

    gsl_vector_memcpy(fTmpSort,fVars);
    gsl_sort_vector(fTmpSort);
    double medianVars;
    if (fVars->size%2)
      medianVars=gsl_vector_get(fTmpSort,fVars->size/2);
    else
      medianVars=(gsl_vector_get(fTmpSort,fVars->size/2-1)+gsl_vector_get(fTmpSort,fVars->size/2)/2);
    // now the median of expected variance

    for (unsigned int i=0;i<fVals->size; i++) {
      double sumx=0;
      for (int j=0;j<fNParams; j++)
	sumx+=gsl_matrix_get(fX,i,j);
      // Dark signal is expected to be sum(Xi)+ron
      gsl_vector_set(fVarX,i,sumx+9.0);
    }
    gsl_vector_memcpy(fTmpSort,fVarX);
    gsl_sort_vector(fTmpSort);
    double medianX;
    if (fVars->size%2)
      medianX=gsl_vector_get(fTmpSort,fVars->size/2);
    else
      medianX=(gsl_vector_get(fTmpSort,fVars->size/2-1)+gsl_vector_get(fTmpSort,fVars->size/2)/2);

    for (unsigned int i=0;i<fVals->size; i++)
      gsl_vector_set(fVars,i,gsl_vector_get(fVarX,i)*medianVars/medianX);

  }

  // need now to store weights since some of them may be killed by early outlier catching... This was experimental, forget about it !
  if (fAutoRescale)
    for (unsigned int i=0;i<fVals->size; i++)
       gsl_vector_set(fWeights,i,1.);
  else 
    for (unsigned int i=0;i<fVals->size; i++)
       gsl_vector_set(fWeights,i,1./gsl_vector_get(Vars,i));
  fNdf = fVals->size - fVal->size;

  //fP = gsl_permutation_calloc(fVals->size);
#ifdef BOFBOF
  if (fDark==2){
    // Experimental : try to check first in order of signal value for outliers
    // Cosmic rate evaluated around 1%(in pixels) for 1H exp.
    for (unsigned int i=0;i<fVals->size; i++)
      gsl_vector_set(fTmpSort,i,fabs(gsl_vector_get(fVals,i)))
    gsl_sort_vector_index(fP,fTmpSort);
  }
#endif


  fChi2=Fit();
  // throw-away
  while (RemoveFar()) {
    fChi2=Fit();
  }
  //gsl_permutation_free(fP);  
	   
  gsl_vector_memcpy(Val,fVal);
  gsl_matrix_memcpy(Covar,fCovar);

}

/* ----- Fit  -------------------------------------------------- */
double KFLinear::Fit() {
  // Fitting Y = Sum ( a_i x_i)
  double chi2;

  gsl_multifit_wlinear(fX,fWeights,fVals,fVal,fCovar,&chi2,fWspace);
  // Bad idea : at the end, all shall correspond to the initial weights
  //if (fAutoRescale) {
  //  gsl_matrix_scale(fCovar,chi2/fNdf);
  //}


  return chi2;
}

/* ----- RemoveFar  -------------------------------------------------- */
int KFLinear::RemoveFar() {

  if (fSigma<=0 || fNdf<=1 )
    return 0;

  int nRm=-1;
  double MaxSigma=0;
  fNdf-=1; // this the hypothesis to be tested on all fits here
  double chi2=-1; // test to know if a fit was performed at all

  for ( unsigned int n=0;n<fVals->size;n++) {
    // excluded already
    if (gsl_vector_get(fWeights,n)==0.0)
      continue;
    if (fDark==2 and fabs(gsl_vector_get(fVals,n)/gsl_vector_get(fVars,n))< fSigma)
      // very unlikely that data compatible with 0 would be in fact a cosmic
      continue;
    
    // compute values with current value substracted
    double wsave=gsl_vector_get(fWeights,n);
    gsl_vector_set(fWeights,n,0.0);
    chi2=Fit();
    double scale=1.;
    if (fAutoRescale)
      scale=chi2/fNdf;
      gsl_matrix_scale(fCovar,scale);
    gsl_vector_set(fWeights,n,wsave);

    // Get the model value
    double model,modelerr;
    gsl_vector_view x=gsl_matrix_row(fX,n);
    gsl_multifit_linear_est(&x.vector,fVal,fCovar,&model,&modelerr);

    double sigma;
    if (fWeightsOnlyForSelection)
      sigma = fabs(gsl_vector_get(fVals,n) - model)/sqrt( modelerr*modelerr + gsl_vector_get(fVars,n));
    else
      sigma = fabs(gsl_vector_get(fVals,n) - model)/sqrt( modelerr*modelerr + scale/gsl_vector_get(fWeights,n));
    if (sigma > MaxSigma && sigma > fSigma) {
      nRm = n;
      MaxSigma = sigma;
    }
  }

  if (nRm>=0) {
    gsl_vector_set(fWeights,nRm,0.0);
    return 1;
  }
  // Restore data
  fNdf += 1;
  if (chi2 != -1)
    Fit(); // not really optimal... better save fVal and fCovar
  return 0;
}

