#/* === Doxygen Comment ======================================= */
#/*! 
# * \file          preprocessor.cxx
#  * \copyright     (c) 2003 SNIFS-Supernova Factory Experiment
#  * \date          Wed Aug  6 17:44:35 2003
#  * \author        Emmanuel Gangler <e.gangler@ipnl.in2p3.fr>
#  * \version       0.0
#  * \brief         
#  *                
#  * $Id$
#  **/
# /* =========================================================== */

#/* ----- local includes ----- */
from bichip import BiChipSnifs
from imagesnifs import ImageSnifs
from section import Section
import utils
import algocams
from overscan import OverscanSnifs, OverscanFromData
import IFU_io
import fits_primary_hd
from items import Anyfile

#/* Note about the fast mode : it suppresses and Variance computation */

#/* ##### Preprocessor ################################################# */

#/* ===== constructor/destructor ======================================= */

#/* ----- constructor ---------------------------------------- */
class Preprocessor:
  def __init__(self):
      self.fMode = 'kIoPlain' #im just using strings instead of enums
      self.fOverscanSnifs = OverscanSnifs()
      self.fOverscanRescue = OverscanFromData()
      self.fOverscan = self.fOverscanSnifs # default current overscan
      self.fFast=0
      self.fAllImage=0

  def SetIoMethod(self,mode) :self.fMode = mode
  def Overscan(self): return self.fOverscan
  def SetFastMode(self, IsFast):self.fFast = IsFast
  def FastMode(self): return  self.fFast
  def SetAllImage(self, IsAll): self.fAllImage = IsAll
  def GetAllImage(self): return self.fAllImage



#/* ----- destructor ---------------------------------------- */
#NOT needed because py has a garbage collector
#Preprocessor::~Preprocessor () {
#  if (fOverscanSnifs)
#    delete fOverscanSnifs;
#  if (fOverscanRescue)
#    delete fOverscanRescue;
#}

#/* ===== Setters ======================================= */

#/* ----- SetOverscanAuto -------------------------------------------------- */
  def SetOverscanAuto(self, Image):
#/* auto - toggling of the overscan between regular / rescue */
      self.fOverscan = self.fOverscanSnifs

      if (Image.HasOverscan()):
          self.fOverscan = self.fOverscanRescue
  
  



#/* ===== Methods ======================================= */

#/* ----- BuildRawBiChip -------------------------------------------------- */
  def BuildRawBiChip(self, name,  outName):
  #returns the original bichip for detcom image
  #reconstructs a bichip from raw data for otcom image.

  #detcom bichips are with extensions [chip00] or must contain a %d in name
  #so we check the %
  #by default, the bichip is loaded read only in totality
  #build the names
        # imname[lg_name+1] 
    imname = ""

    if ( not outName): 
        imname = name
        imname = utils.ut_build_tmp_name(imname,"bichip")
    else:
        imname = outName

    bichip = None

      # // First case : DETCOM

      # // We could have made a copy of the bichip at this stage,
      # // but we want to allow for an utilitary routine converting
      # // otcom.detcom without time overhead in case of an original
      # // detcom image
    if (utils.ut_is_bichip_detcom(name)):
        bichip = BiChipSnifs(name,"I",self.fMode,'kIoAll');
        # // not allowed because I mode
        # //    bichip.SetAlgo("DETCOM");
        # // no new image asked for
        if (not outName or name != outName):
          return bichip
        else:
          out = BiChipSnifs(bichip,outName,'FLOAT',1,self.fMode,'kIoAll')
          del bichip

          # // now try to know if we have CFHT of SNFactory detcom
          primary_name = "" 
          acqswv = ""
          primary_name = utils.ut_primary_header_name(name,primary_name)
          primary_header = Anyfile()
          primary_header = fits_primary_hd.open_primary_hd(primary_header,primary_name,"I")
          IFU_io.disable_user_warnings()
          retCode = IFU_io.RD_desc(primary_header,"ACQSWV",'CHAR',80+1,acqswv)
          IFU_io.restore_user_warnings()
          fits_primary_hd.close_primary_hd(primary_header)
          if (retCode > 0):
              out.SetAlgo("SNFDETCOM")
          else:
              out.SetAlgo("DETCOM")

          return out
            
      
        
        # // Default : OTCOM
        # //
        # // otcom image is already assembled, with bias put somewhere at the end
        # // so I assume :
        # // chip0 1.1024 chip1 1024.1 bias0 1025.1057 bias 1 1057.1025

    image = ImageSnifs(name)

    # // read necessary values for the image
    lg_name = 80
    dataSecString = ''
    biasSecString = ''
    image.RdDesc("DATASEC",'CHAR',lg_name+1,dataSecString)
    image.RdDesc("BIASSEC",'CHAR',lg_name+1,biasSecString)
    dataSec = Section(dataSecString)
    biasSec = Section(biasSecString)

    # // detect if it is a raster : 
    # // in case of a raster, substract 1 column of bad data. 
    # // (last one 'before'(=in 1 channel readout sequence) the overscan )
    isRaster=0
    nAmp = 0
    ccdSecString = ''
    image.RdDesc("CCDSEC",'CHAR',lg_name+1,ccdSecString)
    image.RdDesc("CCDNAMP",'INT',1,nAmp)
    ccdSec = Section(ccdSecString)
    if (ccdSec.XLength()!=1024*nAmp):
        print("Preprocessor::BuildRawBiChip detected a Raster image");
        isRaster=0#;// =0 in case of a bias. Or we just keep all
        # // the reason is : we always keep the last column, even if this is garbadge
        # // it will be garbadge each time the raster doesn't include the last column
        # // before the overscan
        # // BUT : the garbadge there is the only evidence after every computation 
        # // that the image is a non-continuous image, so we keep it
    
    # // build the bichip
    newNamp=nAmp
    if (nAmp != 2 and self.GetAllImage()):
        bichip =  BiChipSnifs(nAmp)
    else:
        bichip = BiChipSnifs(2)
        newNamp=2
    
    im = [None] * bichip.NChips() #imagesnifs

    # // build now the images and the headers
    for chip in range(bichip.NChips()):
        extName = ''
        newDataXLength=dataSec.XLength()/nAmp - isRaster
        newBiasXLength=biasSec.XLength()/nAmp

        im[chip] = ImageSnifs(self.fMode,'kIoAll')
        print(extName,"%s[chip0%d]" % (imname,chip))
        im[chip].CreateFrame(extName,newDataXLength + newBiasXLength, dataSec.YLength())
        im[chip].ImportHeader(image)
        gain = 0.0
        gainKey = ""
        print(gainKey,"CCD%dGAIN" % chip)
        gain = image.RdDesc(gainKey,'DOUBLE',1,gain) #TODO make sure we waant to assign to gain
        im[chip].WrDesc("GAIN",'DOUBLE',1,gain)

        if (newNamp != nAmp):
          im[chip].WrDesc("CCDNAMP",'INT',1,newNamp)
        
        print(dataSecString,"[%d:%d,%d:%d]"%(1,newDataXLength,1,dataSec.YLength()))
        im[chip].WrDesc("DATASEC",'CHAR',lg_name+1,dataSecString);
        print(biasSecString,"[%d:%d,%d:%d]"%(newDataXLength+1,newDataXLength+newBiasXLength,1,dataSec.YLength()))
        
        im[chip].WrDesc("BIASSEC",'CHAR',lg_name+1,biasSecString)
        bichip.SetChip(chip,im[chip])

        # // various other keywords to uniformize
        aKey = ''
        aString = ''
        print(aKey,"CCDSEC%d" %chip)
        if (image.RdIfDesc(aKey,'CHAR',lg_name+1,aString)>0):
          im[chip].WrDesc("CCDSEC",'CHAR',lg_name+1,aString);
        else:# // old data, we have to get sync
          image.RdDesc("CCDSEC",'CHAR',lg_name+1,aString)
          codSec = Section(aString);
          print(aString,"[%d:%d,%d:%d]"%(ccdSec.X1()+newDataXLength*chip,ccdSec.X1()+newDataXLength*(chip+1)-1,ccdSec.Y1(),ccdSec.Y2()))
          im[chip].WrDesc("CCDSEC",'CHAR',lg_name+1,aString);      
        

        print(aKey,"AMPSEC%d" % chip)
        if (image.RdIfDesc(aKey,'CHAR',lg_name+1,aString)>0):
          im[chip].WrDesc("AMPSEC",'CHAR',lg_name+1,aString)
          print(aKey,"DETSEC%d" % chip)
          image.RdDesc(aKey,'CHAR',lg_name+1,aString)
          im[chip].WrDesc("DETSEC",'CHAR',lg_name+1,aString)
        
        

        # // fill the image
        sec = Section() #section
        if (chip%2==0):
          sec.SetX1(dataSec.X1()+(chip*dataSec.XLength())/nAmp)
          sec.SetX2(sec.X1() + newDataXLength - 1 )
          sec.SetY1(dataSec.Y1())
          sec.SetY2(dataSec.Y2())
          im[chip].ImportSection(image,&sec,1,1,1,1)
          sec.SetX1(biasSec.X1()+(chip*biasSec.XLength())/nAmp)
          sec.SetX2(sec.X1() + newBiasXLength - 1)
          sec.SetY1(biasSec.Y1())
          sec.SetY2(biasSec.Y2())
          im[chip].ImportSection(image,&sec,newDataXLength+1,1,1,1)
        else:
          sec.SetX1(dataSec.X1()+(chip*dataSec.XLength())/nAmp + isRaster)
          sec.SetX2(sec.X1() + newDataXLength - 1)
          sec.SetY1(dataSec.Y1())
          sec.SetY2(dataSec.Y2())
          im[chip].ImportSection(image,&sec,newDataXLength,1,-1,1)
          sec.SetX1(biasSec.X1()+(chip*biasSec.XLength())/nAmp)
          sec.SetX2(sec.X1() + newBiasXLength - 1)
          sec.SetY1(biasSec.Y1())
          sec.SetY2(biasSec.Y2())
          im[chip].ImportSection(image,&sec,newDataXLength+newBiasXLength,1,-1,1)
        

    // needs the images 
    bichip.SetAlgo("OTCOM");
    delete image;
    
    return bichip;


# /* ----- PreprocessOverscan------------------------------------------------- */
BiChipSnifs * Preprocessor::PreprocessOverscan(char* name, char* outName){
  /* handles all the information relevant to the overscan */

  // Preliminary : getting a bichip in IO mode

  // build a temporary name if needed
  char imName[lg_name+1];
  if (!outName[0]) {
    strcpy(imName,name);
    ut_build_tmp_name(imName,"bias");
  } else
    strcpy(imName,outName);

  BiChipSnifs * out = BuildRawBiChip(name,imName);

  # // keyword hacking

  # // Detcom image . make a special header hack
  if (out.Chip(0).Algo().GetId() != kOtcom ) {
    char primary_name[lg_name+1];
    # // for some reason, it is not possible to open the temporary image ...
    # // I guess it is because it was not written properly yet...
    # // also, I pretty bet the primary header was not copied anyway ...
    # //    ut_primary_header_name(out.Chip(0).Name(),primary_name);
    ut_primary_header_name(name,primary_name);
    out.HackFitsKeywords(primary_name);
  } else { 
    # // keywords hacking
    out.HackFitsKeywords();
  }
    
  # //
  # // OK, we have now a working copy of the image !
  # //


  # // variance creation
  if (!FastMode()) {
    out.CreateVarianceFrame();
    out.HandleSaturation();
  }
  
  # // overscan substraction
  if (out.Chip(0).HasOverscan()) {
    # // normal exposure with an overscan
    # // we have to substract odd-even, as it has deep impacts on the fit_background, which is very sensitive to local minimas
    fOverscanSnifs.SetOddEven(1);
    fOverscanSnifs.Correct(out);
  }
  # // rescue procedure
  else {
    fOverscanRescue.Correct(out);
  }

  out.UpdateFClass();
  return out;
}

# /* ----- PreprocessAssemble ------------------------------------------------ */
ImageSnifs * Preprocessor::PreprocessAssemble(char* name, char* outName){

  // simply returns the debiased bichip
  // first builds teh debiassed bichip
  // The opening mode shall be the standard one (temporary creation)
  IoMethod_t mode = fMode;
  SetIoMethod(kIoPlain);
  
  BiChipSnifs * bichip = PreprocessOverscan(name);
  // no hack for the variating gain - no good hack found anyway!
  //bichip.HackGainRatio();

  SetIoMethod(mode);
  ImageSnifs *out = bichip.Assemble(outName,fMode,kIoAll);
  delete bichip;

  return out;  
}

# /* ----- PreprocessFlat ------------------------------------------------ */
ImageSnifs * Preprocessor::PreprocessDark(char* name, char* outName,ImageSnifs* bias, DarkModel *biasModel, DarkModel *darkModel){

  ImageSnifs* out = Preprocess(name,outName,bias,0,0,biasModel,darkModel);
  out.BuildDark();
  return out;
}

# /* ----- PreprocessFlat ------------------------------------------------ */
ImageSnifs * Preprocessor::PreprocessFlat(char* name, char* outName,ImageSnifs* bias, ImageSnifs* dark){

  ImageSnifs* out = Preprocess(name,outName,bias,dark);
  out.BuildFlat();
  return out;
}

# /* ----- Preprocess ------------------------------------------------ */
ImageSnifs* Preprocessor::Preprocess(char* name, char* outName,ImageSnifs *bias,ImageSnifs *dark,ImageSnifs* flat,DarkModel * biasModel, DarkModel* darkModel, ImageStackSnifs* darkStack) {

  ImageSnifs* out = PreprocessAssemble(name, outName);
  if (bias) 
    out.SubstractBias(bias);

  out.AddPoissonNoise();

  if (biasModel)
    out.SubstractBiasModel(biasModel);
  if (dark) 
    out.SubstractDark( dark);
  if (darkStack)
    out.SubstractDarkMap(darkStack,darkModel);
  if (darkModel)
    out.SubstractDarkModel(darkModel);

  // do it before multiplicative issues.
  out.HandleCosmetics();

  if (flat) 
    out.ApplyFlat(flat);
  else if (!FastMode())
    out.CustomFlat();
  return out;
}

