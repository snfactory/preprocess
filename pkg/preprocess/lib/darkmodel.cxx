/* === Doxygen Comment ======================================= */
/*! 
 * \file          darkmodel.cxx
 * \copyright     (c) 2003 SNIFS-Supernova Factory Experiment
 * \date          Wed Aug  6 17:44:35 2003
 * \author        Emmanuel Gangler <e.gangler@ipnl.in2p3.fr>
 * \version       0.0
 * \brief         
 *                
 * $Id$
 **/
/* =========================================================== */


/* ----- system includes ----- */
#include <stdio.h>

/* ----- global includes ----- */
#include "IFU_io.h"

/* ----- local includes ----- */
#include "darkmodel.hxx"
#include "section.hxx"

/* ##### DarkModel ################################################# */

/* ===== constructor/destructor ======================================= */

/* ----- DarkModel -------------------------------------------------- */
DarkModel::DarkModel(char* DarkFile) {
  FILE * f=fopen(DarkFile,"r");
  fscanf(f,"%d",&fNsec);
  char secstring[lg_name+1];
  fSections=new Section*[fNsec];
  fI0=new double [fNsec];
  fI1=new double [fNsec];
  fBeta=new double [fNsec];
  fI2=new double [fNsec];

  for (int i=0;i<fNsec;i++){
    fscanf(f,"%s",secstring);
    fSections[i]=new Section(secstring);
    fscanf(f,"%lf,%lf,%lf,%lf",fI0+i,fI1+i,fBeta+i,fI2+i);
  }
}

/* ----- DarkModel -------------------------------------------------- */
DarkModel::~DarkModel() {
  for (int isec=0;isec<fNsec;isec++)
    delete fSections[isec];
  delete[] fSections;
  delete[] fI0;
  delete[] fI1;
  delete[] fBeta;
  delete[] fI2;
}

/* ===== computations ======================================= */

/* ----- DarkSub -------------------------------------------- */
double DarkModel::DarkSub(double Temp, double Timeon, double Texp, int i) {
  return (GetI0(i)+GetI2(i)*TempTerm(Temp))*Texp+GetI1(i)*DarkTimeTerm(Timeon,Texp,i);

}

/* ----- BiasSub -------------------------------------------- */
double DarkModel::BiasSub(double Temp, double Timeon, int i) {
  return (GetI0(i)+GetI1(i)*BiasTimeTerm(Timeon,i)+GetI2(i)*TempTerm(Temp));
}

/* ----- DarkTimeTerm -------------------------------------------- */
double DarkModel::DarkTimeTerm(double Timeon, double Texp, int i) {
  return TimeTerm(Timeon-Texp,Timeon,i);
}

/* ----- BiasTimeTerm -------------------------------------------- */
double DarkModel::BiasTimeTerm(double Timeon, int i) {
  const double kTread=40;
  return TimeTerm(Timeon,Timeon+kTread,i);
}

/* ----- TimeTerm -------------------------------------------- */
double DarkModel::TimeTerm(double Tbeg, double Tend, int i) {
  if ( GetBeta(i) == -1 )
    return  log(Tend/Tbeg);
  else
    return 1.0/(GetBeta(i)+1)*(pow(Tend,(GetBeta(i)+1)) - pow(Tbeg,(GetBeta(i)+1)));
}

/* ----- TempTerm -------------------------------------------- */
double DarkModel::TempTerm(double Temp) {
  const double kBoltz=8.6173e-5; // eV.K-1
  const double kTabs=273.15; // K
  double egap=1.11557-7.021e-4*pow((Temp+kTabs),2)/(1108.+Temp+kTabs);
  return pow((kTabs+Temp),1.5)*exp(-egap/2/kBoltz/(kTabs+Temp));
}

