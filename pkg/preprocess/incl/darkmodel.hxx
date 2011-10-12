/* === Doxygen Comment ======================================= */
/*! 
 * \file          darkmodel.hxx
 * \copyright     (c)  2003 SNIFS-Supernova Factory Experiment
 * \date          Wed Aug  6 18:32:01 2003
 * \author        Emmanuel Gangler <e.gangler@ipnl.in2p3.fr>
 * \version       0.0
 * \brief         
 *                
 * $Id$
 **/
/* =========================================================== */

#ifndef DARKMODEL_H
#define DARKMODEL_H
/*
 Dark Model : an utility calss to read and apply dark model parameters
 the DarkFile should read
 int [the nuber of sections]
 section (string)
 4 doubles (I0,I1,Beta,I2)
 and so on
*/

class Section;

/* ===== DARKMODEL ============================== */

class DarkModel {
public :

  DarkModel(char* DarkFile);
  ~DarkModel();
  
  /* ----- Getters ---------------------------------------- */
  double GetI0(int i) {return fI0[i];}
  double GetI1(int i) {return fI1[i];}
  double GetBeta(int i) {return fBeta[i];}
  double GetI2(int i) {return fI2[i];}
  int GetNsec() {return fNsec;}
  Section** GetSections() {return fSections;}

  double DarkSub(double Temp, double Timeon, double Texp,int i);
  double BiasSub(double Temp, double Timeon,int i);

  double DarkTimeTerm(double Timeon,double Texp, int i);
  double BiasTimeTerm(double Timeon, int i);
  double TimeTerm(double Tbeg,double Tend, int i);
  double TempTerm(double Temp);
  
protected :
  Section ** fSections;
  int fNsec;
  double * fI0, *fI1, *fBeta, *fI2;

};

#endif
