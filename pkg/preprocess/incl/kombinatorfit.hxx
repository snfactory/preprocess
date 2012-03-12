/* === Doxygen Comment ======================================= */
/*! 
 * \file          kombinatorfit.hxx
 * \copyright     (c) 2003 CRAL-Observatoire de Lyon
 * \date          Wed Aug  6 18:32:01 2003
 * \author        Emmanuel Gangler <e.gangler@ipnl.in2p3.fr>
 * \version       0.0
 * \brief         
 *                
 * $Id$
 **/
/* =========================================================== */

#ifndef KOMBINATORFIT_H
#define KOMBINATORFIT_H
/*   
     Kombinator fit is an abstract class which takes a set of data
     and possibly of estimated variances
     And then makes an estimate of some parameters according to
     an additinnal 'X' coordinate
*/

#include <vector>
using namespace std;
#include "gsl/gsl_vector.h"
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_multifit.h"
#include "gsl/gsl_permutation.h"

/* ===== KOMBINATOR ======================================== */

class KombinatorFit  {
  public :
    virtual ~KombinatorFit(){};
  // the output variables should be allocated by the client code
  // Val is of rank N and CovMat of rank N^2
  virtual void KombineFit(vector<double> *X, vector<double> *Vals, vector<double> *Vars, double * Val, double * CovMat)=0;

    virtual int NParam()=0;  // number of parameters for the fit

    virtual int FillsVarOut() =0;
    virtual int NeedsVarIn() =0;
};

/* ===== GAIN KOMBINATOR ======================================== */

/* This is a gain measurment for the region where poisson noise is dominant.
   The X parameter should be linear with exposure time 
   and the Y parameter is the value of the */

class KFGain : public KombinatorFit {
  public : 
    KFGain(double MaxSigmaOutlier);
  
    void KombineFit(vector<double> * X, vector<double> * Vals, vector<double> * Vars, double * Val, double * Var);
    int NeedsVarIn(){return 0;}
    int FillsVarOut(){return 0;}
    int NParam() {return 1;}

  protected :
    double Mean();
    int RemoveFar();

    double fSigma;
    vector<double> *fX, * fVals,fVars;
    double fVal, fSlope;  
};

/*   
     Kombinator fit ND is an abstract class which takes a set of data
     and possibly of estimated variances
     And then makes an estimate of some parameters according to
     any number of additinnal 'X' coordinates
*/

#include <vector>
using namespace std;

/* ===== KOMBINATOR FIT ND ======================================== */

class KombinatorFitND  {
  public :
    virtual ~KombinatorFitND(){};
  // the output variables should be allocated by the client code
  // Val is of rank N and CovMat of rank N^2
  // X is a vector of rank n vectors of rank m values
  // n being the rank of Vals vector and m the number of input variables
  // in order to use native stuff, we use gsl_ vector format...
    virtual void KombineFit(gsl_matrix *X, gsl_vector *Vals, gsl_vector *Vars, gsl_vector * Val, gsl_matrix * CovMat)=0;
    virtual void WorkspaceAlloc(int N)=0;
    virtual void WorkspaceFree()=0;


    virtual int NParam()=0;  // number of parameters for the fit

    virtual int FillsVarOut() =0;
    virtual int NeedsVarIn() =0;
};

/* ===== DARK KOMBINATOR ======================================== */

/* This is what is needed to fit the residual dark map
   X parameters are : texp, T, Ton
   Dark model residual is D = texp * (d0 + Dark(T) * d3 ) + d1*Dark(ton,texp) +d2*Dark(ton,texp)**2
   */

/* ===== LINEAR KOMBINATOR ======================================== */

/* This is a generic linear combinator Y = sum(a_i x_i)
   */

class KFLinear : public KombinatorFitND {
  public : 
  KFLinear(double MaxSigmaOutlier, int NParams, int AutoRescale=0, int WeightsOnlyForSelection=0, int Dark=0);
    ~KFLinear();
  
  void WorkspaceAlloc(int N);
  void WorkspaceFree();


    void KombineFit(gsl_matrix * X, gsl_vector * Vals, gsl_vector * Vars, gsl_vector * Val, gsl_matrix * CovMat);
    int NeedsVarIn(){return 1-fAutoRescale+fWeightsOnlyForSelection;}
    int FillsVarOut(){return 1;}
    int NParam() {return fNParams;}


  protected :
    double Fit();
    int RemoveFar();

    double fSigma;
    int fNParams;
    gsl_matrix *fX, *fCovar;
    gsl_vector *fVals,*fVars, *fWeights, *fVal;
    double fChi2;  
    int fNdf ;

    gsl_vector *fTmpSort, *fVarX;
    gsl_multifit_linear_workspace * fWspace;
    gsl_permutation * fP;
  
    int fAutoRescale, fWeightsOnlyForSelection, fDark;

};


#endif
