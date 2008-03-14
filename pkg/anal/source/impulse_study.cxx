/* === Doxygen Comment ======================================= */
/*! 
 * \file          impulse_study.cxx
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
#include <TProfile2D.h>

/* ----- image includes ------------------------------ */
#include "image.hxx"
#include "catorfile.hxx"
#include "sectionlist.hxx"
#include "section.hxx"
#include "utils.h"

/* ----- local includes ------------------------------ */
#include "roothistos.hxx"


/*------------ main ----------------------*/
int main(int argc, char **argv) {

  char **argval, **arglabel;
  ImageSnifs* in;
  TFile *out;
  char inName[lg_name+1], outName[lg_name+1];
  double sigma;

  set_arglist("-in none -out none -secs Data[1:2048,1:4096] -sigma 0 -range 0,0");
  init_session(argv,argc,&arglabel,&argval);

  CatOrFile inCat(argval[0]);
  CatOrFile outCat(argval[1]);
  SectionList secs(argval[2]);
  get_argval(3,"%lf", &sigma);
  double start,end;
  sscanf(argval[4],"%lf,%lf",&start,&end);

  /* loop on selected sections */
  out = new TFile(argval[1],"UPDATE");

  for (Section* sec = secs.First(); sec; sec = secs.Next()) {

    int maxx=20,minx=-20,nbiny=80,miny=1.0,maxy=5.0;
    char histName[lg_name+1];
    int nbinx=maxx-minx+1;
    sprintf(histName,"Impulse%s",sec->Name());
    TProfile2D* prof = new TProfile2D(histName,"Impulse response",nbinx,minx-0.5,maxx+0.5,nbiny,miny,maxy);
    int maxData=20;
    double data[nbinx][nbiny][maxData];
    int n[nbinx][nbiny];
    for(int i=0;i<nbinx;i++)
      for(int j=0;j<nbiny;j++)
        n[i][j]=0;

    /* loop on catalog */
    while(inCat.NextFile(inName)) {
      printf("Opening now %s\n",inName);
      in = new ImageSnifs(inName);
      
      /* Fill the desired histos */
      for (int j=sec->YFirst(); j<sec->YLast();j++){
        for (int i=sec->XFirst()-miny; i<sec->XLast()-maxy-1;i++){
          // selection of events
          if (in->RdFrame(i,j)>in->RdFrame(i+1,j)*2 && in->RdFrame(i,j)>in->RdFrame(i-1,j)*2 && in->RdFrame(i,j)>10 ){
            double y=log10(in->RdFrame(i,j));
            for (int ix=minx;ix<maxx+1;ix++) {
              int binx=prof->GetXaxis()->FindFixBin(ix)-1;
              int biny=prof->GetYaxis()->FindFixBin(y)-1;
              if (binx<0 || biny<0 || binx>=nbinx || biny>=nbiny)
                continue;
              if (n[binx][biny]==maxData) {
                double z=ut_median(data[binx][biny],n[binx][biny]);
                prof->Fill(ix,y,z);
                n[binx][biny]=0;
              }
              data[binx][biny][n[binx][biny]]=in->RdFrame(i+ix,j);
              n[binx][biny]++;
            }
          }
        }
      }

      // FillHistoLine(in,sec); // cumbersome... use only with small secs
      delete in; 
    }
    prof->Write();
    delete prof;
  }
  out->Purge();
  out->Close();
  delete out;
  
  exit_session(0);
  return(0);
}
