/* === Doxygen Comment ======================================= */
/*! 
 * \file          section.cxx
 * \copyright     (c) 2003 SNIFS-Supernova Factory Experiment
 * \date          Wed Aug  6 17:44:35 2003
 * \author        Emmanuel Gangler <e.gangler@ipnl.in2p3.fr>
 * \version       0.0
 * \brief         
 *                
 * $Id$
 **/
/* =========================================================== */

#include <stdio.h>
#include "section.hxx"

/* ##### Section ################################################# */

/* ===== constructor/destructor ======================================= */

/* ----- Section -------------------------------------------------- */
Section::Section(char* Descr,char* Name) {
  SetString(Descr);
  if (Name)
    strcpy(fName,Name);
}

/* ----- Section -------------------------------------------------- */
Section::Section() {
  fXFirst=fXLast=fYFirst=fYLast=0;
  fName[0]='\0';
}


/* ----- Section -------------------------------------------------- */
Section::Section(int X1, int X2, int Y1, int Y2) {
  fXFirst=X1-1;
  fXLast=X2;
  fYFirst=Y1-1;
  fYLast=Y2;
  fName[0]='\0';
}

/* ----- SetString -------------------------------------------------- */
void Section::SetString(char * Desc) {
  int x1,x2,y1,y2;
  sscanf(Desc,"[%d:%d,%d:%d]",&x1,&x2,&y1,&y2);
  fXFirst=x1-1;
  fXLast=x2;
  fYFirst=y1-1;
  fYLast=y2;

}
