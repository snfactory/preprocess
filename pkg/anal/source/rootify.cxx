/* === Doxygen Comment ======================================= */
/*! 
 * \file          rootify.cxx
 * \copyright     (c) 2003 SNIFS-Supernova Factory Experiment
 * \date          Wed Aug  6 17:44:35 2003
 * \author        Emmanuel Gangler <e.gangler@ipnl.in2p3.fr>
 * \version       0.0
 * \brief         
 *                
 * $Id$
 **/
/* =========================================================== */

/* ----- ROOT includes ------------------------------ */
#include <TFile.h>

/* ----- image includes ------------------------------ */
#include "image.hxx"
#include "catorfile.hxx"
#include "sectionlist.hxx"

/* ----- local includes ------------------------------ */
#include "roothistos.hxx"


/*------------ main ----------------------*/
int main(int argc, char **argv) {

  char **argval, **arglabel;
  ImageSnifs* in;
  TFile *out;
  char inName[lg_name+1], outName[lg_name+1];
  double sigma;
  RootAnalyser ana;

  set_arglist("-in none -out none -secs Data[1:2048,1:4096] -sigma 0 -range 0,0 -sel null");
  init_session(argv,argc,&arglabel,&argval);

  CatOrFile inCat(argval[0]);
  CatOrFile outCat(argval[1]);
  SectionList secs(argval[2]);
  get_argval(3,"%lf", &sigma);
  double start,end;
  sscanf(argval[4],"%lf,%lf",&start,&end);
  ImageSnifs* selection=0;
  if (is_set(argval[5]) )
    selection=new ImageSnifs(argval[5],"I");
  
  

  /* loop on catalog */
  while(inCat.NextFile(inName) && outCat.NextFile(outName)) {
    printf("Opening now %s\n",inName);
    in = new ImageSnifs(inName);
    ana.SetImage(in);
    out = new TFile(outName,"UPDATE");
    
    /* loop on selected sections */
    for (Section* sec = secs.First(); sec; sec = secs.Next()) {
      ana.SetSection(sec);
      /* Fill the desired histos */ 

      // FillHistoLine(in,sec); // cumbersome... use only with small secs
      
      ana.HorizontalProfile(sigma,selection);
      ana.HorizontalMode(selection);
      ana.VerticalProfile(sigma,selection);
      //ana.OddEvenVerticalProfile(sigma);
      //ana.OverscanError(sigma);
      ana.HistoData(start,end,4000,selection);
      ana.MatrixData ();
      ana.HighFrequency();
      //ana.ADCBits(16);
      ana.VertHighFrequencyProf();
    }
    out->Purge();
    out->Close();
    delete out;
    delete in;
  }
  
  exit_session(0);
  return(0);
}
