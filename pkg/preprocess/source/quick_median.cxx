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

#include "imagesnifs.hxx"
#include "analyser.hxx"
#include "catorfile.hxx"
#include "sectionlist.hxx"
#include "section.hxx"

/* ----- main ---------------------------------------- */

int main(int argc, char **argv) {

  char **argval, **arglabel;
  char inName[lg_name+1], secname[lg_name+1];

  set_arglist("-in none -secs [1:1024,1:4096] ");
  init_session(argv,argc,&arglabel,&argval);

  CatOrFile inCat(argval[0]);
  SectionList secs(argval[1]);

  while(inCat.NextFile(inName)) {
    ImageSnifs in(inName);
    ImageAnalyser ana(&in,0);
    for (Section* sec = secs.First(); sec; sec = secs.Next()) {

      ana.SetSection(sec);

      sec->GetString(secname);
      printf("%s %s %f\n",in.Name(),secname, ana.Quantile(0.5) );
    }
  }
  

  exit_session(0);  
}
