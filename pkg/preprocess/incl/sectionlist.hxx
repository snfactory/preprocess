/* === Doxygen Comment ======================================= */
/*! 
 * \file          sectionlist.hxx
 * \copyright     (c)  2003 SNIFS-Supernova Factory Experiment
 * \date          Wed Aug  6 18:32:01 2003
 * \author        Emmanuel Gangler <e.gangler@ipnl.in2p3.fr>
 * \version       0.0
 * \brief         
 *                
 * $Id$
 **/
/* =========================================================== */

#ifndef SECTIONLIST_H
#define SECTIONLIST_H
/*
 A list of sections. Should ideally be build from a <list>
 but here, only a prototype.
The descriptor in constructor is formatted like  <%s[%d,%d:%d,%d]> n times, 
n being the number of sections.
*/

class Section;

/* ===== SECTION ============================== */

class SectionList {
public :

  SectionList(const char* Descr);
  ~SectionList();
  
  /* ----- Getters ---------------------------------------- */
  Section * First();
  Section * Next();
  
  /* ----- Setters ---------------------------------------- */
  // N

protected :
  Section ** fList;
  int fNSec,fCurSec;

};

#endif
