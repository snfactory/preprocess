/* === Doxygen Comment ======================================= */
/*! 
 * \file          Kombinator.hxx
 * \copyright     (c) 2003 CRAL-Observatoire de Lyon
 * \date          Wed Aug  6 18:32:01 2003
 * \author        Emmanuel Gangler <e.gangler@ipnl.in2p3.fr>
 * \version       0.0
 * \brief         
 *                
 * $Id$
 **/
/* =========================================================== */

#ifndef KOMBINATOR_H
#define KOMBINATOR_H
/*   
     Kombinator is an abstract class which takes a set of data
     and possibly of estimated variances
     And then makes a single number estimate, using a virtual recepee
*/

#include <vector>
using namespace std;

/* ===== KOMBINATOR ======================================== */

class Kombinator  {
  public :
    virtual ~Kombinator(){};
    virtual void Kombine(vector<double> *Vals, vector<double> *Vars, double * Val, double * Var)=0;
    virtual int NeedsVarIn() =0;
  
};

/* ===== SPREAD KOMBINATOR ======================================== */

/* For this class, the noise is computed from the data spread */

class KSpread : public Kombinator {
  public : 
    KSpread(double MaxSigmaOutlier,int SpreadMode);
  
    void Kombine(vector<double> * Vals, vector<double> * Vars, double * Val, double * Var);
    int NeedsVarIn(){return 0;}
  

  protected :
    double Mean();
    int RemoveFar(double SigmaCut);

    double fSigma;
    int fSpread;
    vector<double> * fVals,fVars;
    double fVal, fVar, fSum;
  
};

/* ===== GAUSS KOMBINATOR ======================================== */

/* For this class, the noise is gaussian */

class KGauss : public Kombinator {
  public : 
    KGauss(double MaxSigmaOutlier);
  //    ~KGauss() {}
  
    void Kombine(vector<double> * Vals, vector<double> * Vars, double * Val, double * Var);
    int NeedsVarIn(){return 1;}
  

  protected :
    double WeightedMean();
    int RemoveFar(double SigmaCut);
  

  double fSigma;
  
  vector<double> * fVals, * fVars, fWVars;
  double fVal, fVar;
  
};

/* ===== GAUSS + POISSON KOMBINATOR ======================================== */

/* For this class, the estimated Vars contain estimated gaussian + poisson noise */

class KGaussPoisson : public Kombinator {
  public : 
    KGaussPoisson(double MaxSigmaOutlier);
  //~KGaussPoisson() {}
    void Kombine(vector<double> * Vals, vector<double> * Vars, double * Val, double * Var);
    int NeedsVarIn(){return 1;}
  

  protected :
    double WeightedMean();
    void InitWVars();
    void UpdateWVars(double Val);
    int RemoveFar(double SigmaCut);
  

    double fSigma,fPrecision;
    int fMaxIter;
  
    vector<double> * fVals, * fVars, fWVars;
    double fVal, fVar;
  
};


#endif
