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

  set_arglist("-in none -out none -secs Data[4:1027,1:4102]Ovsc[1027:1056,1:4128] -sigma 0");
  init_session(argv,argc,&arglabel,&argval);

  CatOrFile inCat(argval[0]);
  CatOrFile outCat(argval[1]);
  SectionList secs(argval[2]);
  get_argval(3,"%f", &sigma);

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
      ana.HorizontalProfile(sigma);
      ana.VerticalProfile(sigma);
      // OddEvenVerticalProfile(in,sec,sigma);
      ana.OverscanError(sigma);
      ana.HistoData();
      ana.MatrixData ();
      ana.HighFrequency();
    }
    out->Purge();
    out->Close();
    delete out;
    delete in;
  }
  
  exit_session(0);
  return(0);
}



