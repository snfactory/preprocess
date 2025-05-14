# /* === Doxygen Comment ======================================= */
# /*! 
#  * \file          image.cxx
#  * \copyright     (c) 2003 CRAL-Observatoire de Lyon
#  * \date          Wed Aug  6 17:44:35 2003
#  * \author        Emmanuel Gangler <e.gangler@ipnl.in2p3.fr>
#  * \version       0.0
#  * \brief         
#  *                
#  * $Id$
#  **/
# /* =========================================================== */

from IFU_io import RD_catalog  # TODO where do I find IFU_io?

#/* ##### Cat Or File ################################################# */

#/* ===== constructor/destructor ============================== */

#/* ----- CatOrFile ------------------------------------------- */
class CatOrFile:
    def __init__(self, name):
        if (name==None):
            self.fCat = 0
            self.fFirst = 0
        else:
            self.fName = name
            self.fCat = self.fName.find(".cat") != -1
            self.fFirst=1



#/* ===== method ============================== */

#/* ----- Next ---------------------------------------- */
    def NextFile(self, FileName):
        if (self.fCat):
            return RD_catalog(self.fName,FileName)
        elif (self.fFirst):
            FileName = self.fName #TODO change input string
            self.fFirst=0
            return self.fName #used to be return 1 for a truthy value, intead return string to update input string
    
        return self.fFirst
        

