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

#endif
