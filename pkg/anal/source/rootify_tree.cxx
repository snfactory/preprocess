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

/* ----- image includes ------------------------------ */
#include "imagesnifs.hxx"
#include "bichip.hxx"
#include "catorfile.hxx"

/* ----- local includes ------------------------------ */
#include "analimage.hxx"


/*------------ main ----------------------*/
int main(int argc, char **argv) {

  char **argval, **arglabel;
  char inName[lg_name+1],*outName;
  AnalImageSignature * anal;

  set_arglist("-in none -out none");
  init_session(argv,argc,&arglabel,&argval);

  CatOrFile inCat(argval[0]);
  outName = argval[1];
  TFile *roout = new TFile(outName,"UPDATE");

  anal = new AnalImageSignature (roout);
    
  /* loop on catalog */
  while(inCat.NextFile(inName)) {
    print_msg("Opening now %s\n",inName);
    ImageSnifs *in = new ImageSnifs(inName);
    
    //BiChipSnifs * out= new BiChipSnifs(*in,"mem://tmp.fits",FLOAT,1);

    anal->SetImage(in);
    anal->FillImageSignature();
    anal->StoreImageSignature();
    
    delete in;
  }
  roout->Write();
  roout->Purge();
  roout->Close();
  delete roout;
  delete anal;
  

  exit_session(0);
  return(0);
}
