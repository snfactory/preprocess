/* === Doxygen Comment ======================================= */
/*! 
 * \file          catorfile.hxx
 * \copyright     (c)  2003 SNIFS-Supernova Factory Experiment
 * \date          Wed Aug  6 18:32:01 2003
 * \author        Emmanuel Gangler <e.gangler@ipnl.in2p3.fr>
 * \version       0.0
 * \brief         
 *                
 * $Id$
 **/
/* =========================================================== */

#ifndef CATORFILE_H
#define CATORFILE_H

/*
CatOrFile is an utility to access an input parameter which can be either a 
catalog or a pklain file. The output is a first/next method.
A catalog is a .cat file, a file is the contrary !
*/

/* ===== CATORFILE ============================== */

class CatOrFile {
public :

  CatOrFile(const char* Name);
  ~CatOrFile() {}
  
  int NextFile(char* FileName);

protected :
  char fName[lg_name+1];
  int fCat;
  int fFirst;
};

#endif
