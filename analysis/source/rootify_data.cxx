/* === Doxygen Comment ======================================= */
/*! 
 * \file          rootify_data.cxx
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
#include <TH1.h>

/* ----- image includes ------------------------------ */
#include "image.hxx"
#include "catorfile.hxx"
#include "sectionlist.hxx"
#include "bichip.hxx"
#include "section.hxx"
#include "preprocessor.hxx"
#include "overscan.hxx"

/* ----- local includes ------------------------------ */
#include "roothistos.hxx"


/*------------ main ----------------------*/
int main(int argc, char **argv) {

  char **argval, **arglabel;
  char inName[lg_name+1],*outName;
  RootAnalyser ana;
  int nlines;

  set_arglist("-in none -out none -nval 5000");
  init_session(argv,argc,&arglabel,&argval);

  CatOrFile inCat(argval[0]);
  outName = argval[1];
  TFile *roout = new TFile(outName,"UPDATE");

    
  get_argval(2,"%d", &nlines);
  Preprocessor P;

  /* loop on catalog */
  while(inCat.NextFile(inName)) {
    print_msg("Opening now %s\n",inName);
    BiChipSnifs *in = new BiChipSnifs(inName);
    BiChipSnifs * out= new BiChipSnifs(*in,"mem://tmp.fits",FLOAT,1);
    delete in;
    out-> HackFitsKeywords();
    P.SetOverscanAuto(out->Chip(0));
    P.Overscan()->Correct(out);
    
    for (int chip=0;chip<2;chip++) {

      char tmp_name[lg_name+1],data[lg_name+1],*begin,*end;
      out->Chip(chip)->RdDesc("DATASEC",CHAR,lg_name+1,data);
      Section* sec=new Section(data,"Data");  
      ana.SetImage(out->Chip(chip));
      ana.SetSection(sec);
      //keeps only core of the name
      if (!(begin = strrchr(inName,'/')))
        begin = inName;
      else begin++;
      if ((end = strrchr(inName,'.')))
        end[0] = '\0';

      sprintf(tmp_name,"%s-c0%d",begin,chip);
      TH1F * hist = ana.HistoDataBuild(tmp_name,nlines);
      hist->Write();
      delete hist;
      
    }
    delete out;
  }
  roout->Purge();
  roout->Close();
  delete roout;
  

  exit_session(0);
  return(0);
}
