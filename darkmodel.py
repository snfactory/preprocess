#/* === Doxygen Comment ======================================= */
#/*! 
# * \file          darkmodel.cxx
# * \copyright     (c) 2003 SNIFS-Supernova Factory Experiment
# * \date          Wed Aug  6 17:44:35 2003
# * \author        Emmanuel Gangler <e.gangler@ipnl.in2p3.fr>
# * \version       0.0
#  * \brief         
#  *                
#  * $Id$
#  **/
# /* =========================================================== */


# /* ----- system includes ----- */
# include <stdio.h>

# /* ----- global includes ----- */
# include "IFU_io.h"
import IFU_io
# /* ----- local includes ----- */
from section import Section
import math

# /* ##### DarkModel ################################################# */

# /* ===== constructor/destructor ======================================= */

#/* ----- DarkModel -------------------------------------------------- */
class DarkModel:
  def __init__(self, DarkFile):
      with open(DarkFile,"r") as f:
        # fscanf(f,"%d",&fNsec);
        line = f.readLine().strip().split()
        word = line[0]
        nextWord = 1
        self.fNsec = float(word)
        # char secstring[lg_name+1];
        #all arrays are fNsec long
        self.fSections=[None]*self.fNsec #an array of section 
        self.fI0=[0]*self.fNsec # an array of doubles
        self.fI1= [0]*self.fNsec  #array of doubles
        self.fBeta= [0]*self.fNsec # array of doubles
        self.fI2= [0]*self.fNsec # array of doubles

        for i in range(self.fNsec):
          if nextWord == len(line):
            line = f.readLine().strip().split()
            word = line[0]
            nextWord=1
          else:
            word = line[nextWord]
            nextWord+=1
          secstring = str(word)
          # fscanf(f,"%s",secstring)
          self.fSections[i] = Section(secstring)

          if nextWord == len(line):
            line = f.readLine().strip().split()
            word = line[0]
            nextWord=1
          else:
            word = line[nextWord]
            nextWord+=1
          self.fI0[i] = float(word)
          if nextWord == len(line):
            line = f.readLine().strip().split()
            word = line[0]
            nextWord=1
          else:
            word = line[nextWord]
            nextWord+=1
          self.fI1[i] = float(word)
          if nextWord == len(line):
            line = f.readLine().strip().split()
            word = line[0]
            nextWord=1
          else:
            word = line[nextWord]
            nextWord+=1
          self.fBeta[i] = float(word)
          if nextWord == len(line):
            line = f.readLine().strip().split()
            word = line[0]
            nextWord=1
          else:
            word = line[nextWord]
            nextWord+=1
          self.fI2[i] = float(word)

        
  def GetI0(self, i):return self.fI0[i]
  def GetI1(self, i): return self.fI1[i]
  def GetBeta(self, i): return self.fBeta[i]
  def GetI2(self, i): return self.fI2[i]
  def GetNsec(self): return self.fNsec
  def GetSections(self): return self.fSections
        

# /* ----- DarkModel -------------------------------------------------- */
# DarkModel::~DarkModel() {
#   for (int isec=0;isec<fNsec;isec++)
#     delete fSections[isec];
#   delete[] fSections;
#   delete[] fI0;
#   delete[] fI1;
#   delete[] fBeta;
#   delete[] fI2;
# }

  # /* ===== computations ======================================= */

  # /* ----- DarkSub -------------------------------------------- */
  def DarkSub(self, Temp, Timeon, Texp, i):
    return (self.GetI0(i)+self.GetI2(i)*self.TempTerm(Temp))*Texp+self.GetI1(i)*self.DarkTimeTerm(Timeon,Texp,i);



  # /* ----- BiasSub -------------------------------------------- */
  def BiasSub(self, Temp,Timeon,i):
    return (self.GetI0(i)+self.GetI1(i)*self.BiasTimeTerm(Timeon,i)+self.GetI2(i)*self.TempTerm(Temp));


  # /* ----- DarkTimeTerm -------------------------------------------- */
  def DarkTimeTerm(self, Timeon,  Texp,  i):
    return self.TimeTerm(Timeon-Texp,Timeon,i)


  # /* ----- BiasTimeTerm -------------------------------------------- */
  def BiasTimeTerm(self, Timeon,  i):
    kTread=40
    return self.TimeTerm(Timeon,Timeon+kTread,i)


  # /* ----- TimeTerm -------------------------------------------- */
  def TimeTerm(self, Tbeg,  Tend,  i):
    if ( self.GetBeta(i) == -1 ):
      return  math.log(Tend/Tbeg)
    else:
      return 1.0/(self.GetBeta(i)+1)*(pow(Tend,(self.GetBeta(i)+1)) - pow(Tbeg,(self.GetBeta(i)+1)))


  # /* ----- TempTerm -------------------------------------------- */
  def TempTerm(self, Temp):
    kBoltz=8.6173e-5# // eV.K-1
    kTabs=273.15 # // K
    egap=1.11557-7.021e-4*pow((Temp+kTabs),2)/(1108.+Temp+kTabs)
    return pow((kTabs+Temp),1.5)*pow(math.e, -egap/2/kBoltz/(kTabs+Temp))


