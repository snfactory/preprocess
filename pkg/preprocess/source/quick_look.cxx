/* === Doxygen Comment ======================================= */
/*! 
 * \file          quick_look.cxx
 * \copyright     (c) 2004 SNIFS Collaboration
 * \date          Wed Aug  6 17:44:35 2003
 * \author        Emmanuel Gangler <e.gangler@ipnl.in2p3.fr>
 * \version       0.0
 * \brief         
 *                
 * $Id$
 **/
/* =========================================================== */

/* ----- local include ----- */

#include "preprocessor.hxx"
#include "imagesnifs.hxx"
#include "analyser.hxx"
#include "catorfile.hxx"
#include "section.hxx"
#include "overscan.hxx"

/* ----- main ---------------------------------------- */

int main(int argc, char **argv) {

  char **argval, **arglabel;
  char inName[lg_name+1];
  double sigcut, limit;

  set_arglist("-in none -sec [4:1027,1:4096] -sigcut 5.0 -limit 10000");
  init_session(argv,argc,&arglabel,&argval);

  CatOrFile inCat(argval[0]);
  Section *sec = new Section(argval[1]);
  sscanf(argval[2],"%lf",&sigcut);
  sscanf(argval[3],"%lf",&limit);

  printf("Name Mean Sigma PixOut PixOutMean PixOver Max\n");
  Preprocessor P;

  while(inCat.NextFile(inName)) {
    ImageSnifs in(inName);
    char tmp_name[lg_name+1];
    sprintf(tmp_name,"mem://%s",argval[1]);
    ImageSnifs out(inName,tmp_name,FLOAT,1);
    
    int ovsc;
    if (out.RdIfDesc("OVSCDONE",INT,1,&ovsc)<0 || ovsc==0 ) {
      
      out.HackFitsKeywords();
      P.SetOverscanAuto(&out);
      P.Overscan()->Correct(&out);
    }
    
    ImageAnalyser ana(&out,sec);
    double min,max;

    out.Image()->MinMax(sec,&min,&max);
    
    printf("%s %f %f %d %f %d %f \n",in.Name(), 
           ana.MeanLevel(), sqrt(ana.StatsVariance()),ana.NPixOut(sigcut), ana.OutPixMean(sigcut),ana.NPixOver(limit), max);
  }
  delete sec;

  exit_session(0);  
}
