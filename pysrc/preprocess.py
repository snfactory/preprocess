#/* === Doxygen Comment ======================================= */
#/*! 
# * \file          preprocess.cxx
# * \copyright     (c) 2004 SNIFS Collaboration
# * \date          Wed Aug  6 17:44:35 2003
# * \author        Emmanuel Gangler <e.gangler@ipnl.in2p3.fr>
# * \version       0.0
# * \brief         
# *                
# * $Id: preprocess.cxx,v 1.15 2018/07/10 00:25:16 kponder Exp $
# **/
#/* =========================================================== */

#/* ----- local imports ----- */
import bichip
from imagesnifs import ImageSnifs
from catorfile import CatOrFile
from preprocessor import Preprocessor
from darkmodel import DarkModel
from binaryoffsetmodel import BinaryOffsetModel
from imagestacksnifs import ImageStackSnifs

import gc

#/* ----- main ---------------------------------------- */

if __name__ == "__main__": #(int argc, char **argv) 
  import sys
  args = sys.argv

#  char **argval, **arglabel;
  
 # set_arglist("-in none -out none -bias null -dark null -flat null -fast -all -bm null -dm null -bom null");
  
  #init_session(argv,argc,&arglabel,&argval);
  
  arglabels = ["-in", "-out", "-bias", "-dark", "-flat", "-fast", "-all", "-bm", "-dm", "-bom"]
  argval =    [None,    None,    None,    None,    None,    None,   None,  None,  None,  None]
  for i in range(args): 
    j = arglabels.index(args[i])
    if (j!=-1):
      argval[j] = args[i+1]


  inName = None
  outName = None


  inCat = CatOrFile(argval[0])
  outCat = CatOrFile(argval[1])
  


# Load once auxilliary files
  if ((argval[2])):
    bias = ImageSnifs(argval[2]);
 # // See later
 # //if (is_set(argval[3]))
 # //  dark = new ImageSnifs(argval[3]);
  if ((argval[4])):
    flat = ImageSnifs(argval[4]);

  P = Preprocessor()

  if ((argval[5])):
    P.SetFastMode(1)
  if ((argval[6])):
    P.SetAllImage(1)


  if ((argval[7])):
    biasModel =  DarkModel(argval[7]);
  if ((argval[8])):
    darkModel =  DarkModel(argval[8]);
  if ((argval[9])):
    binaryOffsetModel = BinaryOffsetModel(argval[9]);

#  // Switch between dark subtraciton versions
  dark=None
  darkStack=None
  hduName = None
  if ((argval[3])):
    hduName = f"{argval[3]}[image{0}]"
    if (hduName is not None and darkModel):
      darkStack = ImageStackSnifs(argval[3])
    else:
      dark=  ImageSnifs(argval[3])
    
  inName = inCat.NextFile(inName)
  outName = outCat.NextFile(outName)
  while (inName and outName):
    print(f'Processing {inName}')
    out = P.Preprocess(inName,outName,bias,dark,flat,biasModel,darkModel,darkStack,binaryOffsetModel)
    del out
    inName = inCat.NextFile(inName)
    outName = outCat.NextFile(outName)
  

  if (flat): del flat
  if (dark): del dark
  if (bias): del bias

  gc.collect()

  sys.exit(0)
  
