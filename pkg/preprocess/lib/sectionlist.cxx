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

#include "sectionlist.hxx"
#include "section.hxx"

/* ##### SectionList ################################################# */

/* ===== constructor/destructor ======================================= */

/* ----- SectionList -------------------------------------------------- */
SectionList::SectionList(const char* Descr) {
  const char* psec = Descr;
  fNSec=0;
  // counting the number of tokens
  while (*psec!='\0') {
    fNSec++;
    psec = strchr(psec,']')+1;
  }
  fList = new Section*[fNSec+1];
  psec = Descr;
  for (int isec=0;isec<fNSec;isec++) {
    char bufName[lg_name+1], bufSec[lg_name+1],*pbuf;
    strcpy(bufName,psec);
    *(strchr(bufName,']')+1)='\0'; // cut the word at the end
    strcpy(bufSec,bufName);
    *(strchr(bufName,'['))='\0';   // suppresses the '['
    pbuf= strchr(bufSec,'[');   // beginning of the [...]
    fList[isec] = new Section(pbuf,bufName);
    psec = strchr(psec,']')+1;
  }
  fList[fNSec]=0;
}

/* ----- ~SectionList -------------------------------------------------- */
SectionList::~SectionList() {
  for (int i=0;i<fNSec;i++)
    delete fList[i];
  delete[] fList;
}

/* ===== getters ======================================= */

/* ----- First -------------------------------------------------- */

Section* SectionList::First(){
  fCurSec=0;
  return Next();
}

/* ----- Next -------------------------------------------------- */

Section* SectionList::Next(){
  if (fCurSec<fNSec)
    return fList[fCurSec++];
  else
    return 0;
}

