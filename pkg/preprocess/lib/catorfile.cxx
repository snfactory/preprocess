/* === Doxygen Comment ======================================= */
/*! 
 * \file          image.cxx
 * \copyright     (c) 2003 CRAL-Observatoire de Lyon
 * \date          Wed Aug  6 17:44:35 2003
 * \author        Emmanuel Gangler <e.gangler@ipnl.in2p3.fr>
 * \version       0.0
 * \brief         
 *                
 * $Id$
 **/
/* =========================================================== */

#include "IFU_io.h"

#include "catorfile.hxx"

/* ##### Cat Or File ################################################# */

/* ===== constructor/destructor ============================== */

/* ----- CatOrFile ------------------------------------------- */
CatOrFile::CatOrFile(const char* Name) {
  strcpy(fName,Name);
  fCat = (strstr(fName,".cat") != NULL);
  fFirst=1;
}

/* ===== method ============================== */

/* ----- Next ---------------------------------------- */
int CatOrFile::NextFile(char* FileName) {
  if (fCat) 
    return RD_catalog(fName,FileName);
  else if (fFirst) {
    strcpy(FileName,fName);
    fFirst=0;
    return 1;
  }
  return fFirst;
}

