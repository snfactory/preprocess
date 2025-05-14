# /* === Doxygen Comment ======================================= */
# /*! 
#  * \file          section.cxx
#  * \copyright     (c) 2003 SNIFS-Supernova Factory Experiment
#  * \date          Wed Aug  6 17:44:35 2003
#  * \author        Emmanuel Gangler <e.gangler@ipnl.in2p3.fr>
#  * \version       0.0
#  * \brief         
#  *                
#  * $Id$
#  **/
# /* =========================================================== */

#include <stdio.h>
#include "section.hxx"

# /* ##### Section ################################################# */

# /* ===== constructor/destructor ======================================= */
class Section:
# /* ----- Section -------------------------------------------------- */
    def __init__(self, *args):
        self.fName = ""
        if len(args) == 2:
            self.SetString(args[0])
            self.SetName(args[1])
        elif len(args) == 0:
            self.fXFirst = 0
            self.fXLast= 0
            self.fYFirst=0
            self.fYLast=0
            self.fName = ""
        elif len(args) == 5:
            self.fXFirst=args[0]-1
            self.fXLast=args[1]
            self.fYFirst=args[2]-1
            self.fYLast=args[3]
            self.SetName(args[4])

# /* ----- SetString -------------------------------------------------- */
    def SetString(self, Desc):
        words = Desc.split()
        if (len(words)!=4):
            #print_error what is this function TODO  should we error out?
            print("Section::SetString : %s is ill-formatted. Expecting [%%d:%%d,%%d:%%d]" % Desc)
    
        x1 = float(Desc[0])
        x2 = float(Desc[1])
        y1 = float(Desc[2])
        y2 = float(Desc[3])
        
        self.fXFirst=x1-1
        self.fXLast=x2
        self.fYFirst=y1-1
        self.fYLast=y2



# /* ----- SetString -------------------------------------------------- */
    def SetName(self, Name):
        self.fName = Name
    def XFirst(self): return self.fXFirst
    def XLast(self): return self.fXLast
    def X1(self): return self.fXFirst+1
    def X2(self): return self.fXLast
    def XLength(self): return self.fXLast-self.fXFirst

    def YFirst(self): return self.fYFirst
    def YLast(self): return self.fYLast
    def Y1(self): return self.fYFirst+1
    def Y2(self): return self.fYLast
    def YLength(self): return self.fYLast-self.fYFirst

    def GetString(self):
        return f"[{self.fXFirst+1}:{self.fXLast},{self.fYFirst+1}:{self.fYLast}]"
    def Name(self): 
        return self.fName



    #   /* ----- Setters ---------------------------------------- */
    def SetXFirst(self, X):
        self.fXFirst = X
    def SetXLast(self, X) :
        self.fXLast = X
    def SetYFirst(self,Y):
        self.fYFirst = Y
    def SetYLast(self, Y):
        self.fYLast = Y

    def SetX1(self,X):
        self.fXFirst = X-1
    def SetX2(self,X):
        self.fXLast = X
    def SetY1(self, Y):
        self.fYFirst = Y-1
    def SetY2(self,Y):
        self.fYLast = Y

