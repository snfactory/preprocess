/* === Doxygen Comment ======================================= */
/*! 
 * \file          section.hxx
 * \copyright     (c)  2003 SNIFS-Supernova Factory Experiment
 * \date          Wed Aug  6 18:32:01 2003
 * \author        Emmanuel Gangler <e.gangler@ipnl.in2p3.fr>
 * \version       0.0
 * \brief         
 *                
 * $Id$
 **/
/* =========================================================== */

#ifndef SECTION_H
#define SECTION_H
/*
 Section :
 a small utility to go from [x1:x2,y1:y2] format, 
 with x1/y1 starting at 1 and x2/y2 is the last valid item
 to C-style
 XFirst, XLast,YFirst,YLast, starting at 0

 By convention, X1 starts at 1, XFirst at 0

*/

/* ===== SECTION ============================== */

class Section {
public :

  Section();
  Section(char* Descr);
  Section(int X1, int X2, int Y1, int Y2);
  ~Section(){}
  
  /* ----- Getters ---------------------------------------- */
  int XFirst() {return fXFirst;}
  int XLast() {return fXLast;}
  int X1() {return fXFirst+1;}
  int X2() {return fXLast;}
  int XLength() {return fXLast-fXFirst;}
  
  int YFirst() {return fYFirst;}
  int YLast() {return fYLast;}
  int Y1() {return fYFirst+1;}
  int Y2() {return fYLast;}
  int YLength() {return fYLast-fYFirst;}
  
  void GetString(char* ToFill)
    {sprintf(ToFill,"[%d:%d,%d:%d]",fXFirst+1,fXLast,fYFirst+1,fYLast);}

  /* ----- Setters ---------------------------------------- */
  void SetXFirst(int X) {fXFirst = X;}
  void SetXLast(int X) {fXLast = X;}
  void SetYFirst(int Y) {fYFirst = Y;}
  void SetYLast(int Y) {fYLast = Y;}

  void SetX1(int X) {fXFirst = X-1;}
  void SetX2(int X) {fXLast = X;}
  void SetY1(int Y) {fYFirst = Y-1;}
  void SetY2(int Y) {fYLast = Y;}
  void SetString(char* Desc);

protected :
  int fXFirst,fXLast,fYFirst,fYLast;
  

};

#endif
