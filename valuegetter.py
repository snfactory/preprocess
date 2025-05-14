# /* === Doxygen Comment ======================================= */
# /*! 
#  * \file          valuegetter.cxx
#  * \copyright     (c) 2003 SNIFS-Supernova Factory Experiment
#  * \date          Wed Aug  6 17:44:35 2003
#  * \author        Emmanuel Gangler <e.gangler@ipnl.in2p3.fr>
#  * \version       0.0
#  * \brief         
#  *                
#  * $Id$
#  **/
# /* =========================================================== */

from analyser import ImageAnalyser
import darkmodel

#numpy import to replace GSL
import numpy as np

lg_name = 80
#include "analyser.hxx"
#include "valuegetter.hxx"
#include "darkmodel.hxx"

# /* ##### MeanFromAnalyser ################################################# */

# /* ===== constructor/destructor ======================================= */

# /* ----- constructor  ----------------------------------*/
class ValueAnalyserMean:
    def __init__(self, Sec):
        self.fAnalyser = ImageAnalyser()
        self.fAnalyser.SetSection(Sec)

# /* ----- destructor  ----------------------------------*/
# ValueAnalyserMean::~ValueAnalyserMean(){
#   delete fAnalyser;
# }

# /* ===== methods ======================================= */

# /* ----- GetValue ----------------------------------*/
    def GetValue(self,Image):
        self.fAnalyser.SetImage(Image)
        return self.fAnalyser.MeanLevel()


# /* ##### Value for the dark current ################################################# */

# /* ===== constructor/destructor ======================================= */

# /* ----- constructor  ----------------------------------*/
class ValueDark:
    def __init__(self,darkModel):
        self.fDarkModel = darkModel

# /* ----- destructor  ----------------------------------*/
# ValueDark::~ValueDark(){
# }

# /* ===== methods ======================================= */

# /* ----- GetValue ----------------------------------*/
    def GetValue(self, Image):
    # // returns the dark current estimation
    # ImageSnifs* image = (ImageSnifs*) Image;

    # double timeon,temp,texp;
        image = Image
        temp = image.RdDesc("DETTEMP",'DOUBLE',1); 
        texp = image.RdDesc("DARKTIME",'DOUBLE',1); 

        # char timeOnStr[lg_name+1];
        timeOnStr = image.RdDesc("TIMEON",'CHAR',lg_name+1, timeOnStr)  # // CAVEAT : TIMEON not defined for all data
        if (timeOnStr=="None"):
            timeon = image.RdDesc("TIMEON",'DOUBLE',1,timeon)  # // CAVEAT : TIMEON not defined for all data.
        else:
            timeon=-1
        if (timeon<texp):
            print("WARNING: %s has a bad time on" % image.Name())
            timeon=-1
        

        # //Section ** Secs=fDarkModel->GetSections();
        toremove=0
        for isec  in range(self.fDarkModel.GetNsec()):
            if (timeon>0):
                toremove+=self.fDarkModel.DarkSub(temp,timeon,texp,isec)
            else:
                toremove+=(self.fDarkModel.GetI0(isec)+self.fDarkModel.GetI2(isec)*self.fDarkModel.TempTerm(temp))*texp
        toremove /= self.fDarkModel.GetNsec()

        return toremove
      

# /* ##### ValuesGetterDarkFitter ################################################# */

# /* ===== constructor/destructor ======================================= */

# /* ----- constructor  ----------------------------------*/

class ValuesGetterDarkFitter:
    def __init__(self, Model,activate,  offseton, offsetT):#darkmodel, a int array, double, double
        if (Model.GetNsec() != 1):
            print("ERROR: ValuesGetterDarkFitter : needs a 1-section model")
        self.fDarkModel = Model
        self.fNParams=0
        for i in range(3):
            if (activate and not activate[i]):
                self.fActive[i]=0;      
            else:
                self.fActive[i]=1
                self.fNParams+=1



# /* ===== methods ======================================= */

# /* ----- GetValue ----------------------------------*/
    def GetValues(self, Image, retValues): #ret values is a GSL vector, perhaps use np array?
    #   // returns the dark current estimation

    #cast from simple image to image snifs, I may rewrite the image code to handle this in python
    #   ImageSnifs* image = (ImageSnifs*) Image 
        image = Image

        timeon=0
        temp=0
        texp=0

        temp = image.RdDesc("DETTEMP",'DOUBLE',1,temp); 
        texp = image.RdDesc("DARKTIME",'DOUBLE',1,texp); 

        timeOnStr = image.RdDesc("TIMEON",'CHAR',lg_name+1,timeOnStr)  # // CAVEAT : TIMEON not defined for all data
        if (timeOnStr == "None"): #maybe change this later
            timeon = image.RdDesc("TIMEON",'DOUBLE',1,timeon)  #// CAVEAT : TIMEON not defined for all data.
        else:
            timeon=-1
        if (timeon<texp):
            print("Error, %s has a bad time on" % image.Name())
        

        count=0
        if (self.fActive[0]):
            # gsl_vector_set(retValues,count,fDarkModel->GetI0(0)*texp);
            retValues[count] = self.fDarkModel.getI0(0)*texp
            count+=1
        if (self.fActive[1]):
            # gsl_vector_set(retValues,count,fDarkModel->GetI1(0)*fDarkModel->DarkTimeTerm(timeon,texp,0));
            retValues[count] = self.fDarkModel.GetI1(0)*self.fDarkModel.DarkTimeTerm(timeon,texp,0)
            count+=1
        
        if (self.fActive[2]):
            # gsl_vector_set(retValues,count,fDarkModel->GetI2(0)*fDarkModel->TempTerm(temp)*texp);
            retValues[count] = self.fDarkModel.GetI2(0)*self.fDarkModel.TempTerm(temp)*texp
            count+=1

        #im assuming the C++ is assigning directly to the vector so to update in py we need to return
        return retValues
        

