/* === Doxygen Comment ======================================= */
/*! 
 * \file          imagestacksnifs.cxx
 * \copyright     (c) 2003 IPNL
 * \date          Wed Aug  6 17:44:35 2003
 * \author        Emmanuel Gangler <e.gangler@ipnl.in2p3.fr>
 * \version       0.0
 * \brief         
 *                
 * $Id$
 **/
/* =========================================================== */

/*
This class is a handler for processing a stack of images.
It may be seen as a steering for the whole processing.
But for the moment, only a few functionnalities are handled.
The input catalog is supposed to be of homogeneous
*/

/* #####  ImageStackSnifs ################################################# */

#include "utils.h"
#include "catorfile.hxx"
#include "imagesnifs.hxx"
#include "bichip.hxx"
#include "imagestack.hxx"
#include "kombinator.hxx"
#include "kombinatorfit.hxx"
#include "imagestacksnifs.hxx"


/* ===== constructor/destructor ======================================= */

/* ----- ImageStackSnifs  -------------------------------------------------- */
ImageStackSnifs::ImageStackSnifs(CatOrFile* Cat, char * Mode){
  /* first loads the headers - to check the homogeneity of the catalog */
  /* the ioslice is used, this means the numer of loaded lines has to be 
   reworked for image-based processing */
  
  char fileName[lg_name+1];
  if (!Cat->NextFile(fileName)) {
    print_error("ImageStackSnifs::ImageStackSnifs : unable to read %s",Cat->Name());
    return;
  }

  int nDone=0;
  ImageSnifs * imRef;
  do {
    if (ut_is_bichip(fileName)) {
      print_error("ImageStackSnifs::ImageStackSnifs can not read bichips %s",fileName);
        return; 
      }

    // loads the image
    ImageSnifs * image = new ImageSnifs(fileName,Mode,kIoSlice,1);
    fImages.push_back(image);
    if (!nDone) {
      imRef = image;
    } else 
      if (!imRef->CanBeStackedWith(image)) {
        print_error("ImageStackSnifs::ImageStackSnifs %s is incompatible with $s",image->Name(),imRef->Name());
        return;
      }
    nDone++;  
  } while(Cat->NextFile(fileName));
  
  
}

/* ----- ~ImageStackSnifs  -------------------------------------------------- */

ImageStackSnifs::~ImageStackSnifs(){
  vector<ImageSnifs*>::iterator iter;
  for (iter = fImages.begin();iter != fImages.end();++iter){
    delete *iter;
  }
}


/* ===== method ======================================= */

/* ----- Kombine ---------------------------------------- */
ImageSnifs* ImageStackSnifs::Kombine(char* OutName,  Kombinator * K) {

  ImageSnifs * out = new ImageSnifs(*fImages[0],OutName,FLOAT,0,kIoSlice,fNLinesMem);
  vector<ImageSimple * > imList;
  for (unsigned int i=0;i<fImages.size();i++) {
    imList.push_back(fImages[i]);
  }
  if (!out->Variance()) {
      out->CreateVarianceFrame();
    }
    
  ImageStack images(imList);
  images.SetKombinator(K);
  images.Kombine(out);
  imList.clear();

  return out;
}



/* #####  BiChipStackSnifs ################################################# */


/* ===== constructor/destructor ======================================= */

/* ----- BiChipStackSnifs  ------------------------------------------------ */
BiChipStackSnifs::BiChipStackSnifs(CatOrFile* Cat, char * Mode,int NLines){
  /* first loads the headers - to check the homogeneity of the catalog */
  /* the ioslice is used, this means the numer of loaded lines has to be 
   reworked for image-based processing */
  
  fNLinesMem=NLines;

  char fileName[lg_name+1];
  if (!Cat->NextFile(fileName)) {
    print_error("ImageStackSnifs::ImageStackSnifs : unable to read %s",Cat->Name());
    return;
  }

  int nDone=0;
  ImageSnifs * imRef[2];
  do {
    if (!ut_is_bichip(fileName)) {
      print_error("BiChipStackSnifs::BiChipStackSnifs can not read plain image %s",fileName);
        return; 
      }

    // loads the image
    BiChipSnifs * bichip = new BiChipSnifs(fileName,Mode,kIoSlice,fNLinesMem);
    fBiChips.push_back(bichip);
    if (!nDone) {
      imRef[0] = bichip->Chip(0);
      imRef[1] = bichip->Chip(1);
    } else 
      if (!imRef[0]->CanBeStackedWith(bichip->Chip(0)) ||
          !imRef[1]->CanBeStackedWith(bichip->Chip(1))) {
        
        print_error("ImageStackSnifs::ImageStackSnifs %s is incompatible with teh first one\n",fileName);
        return;
      }
    nDone++;
  } while(Cat->NextFile(fileName));
}

/* ----- BiChipStackSnifs  -------------------------------------------------- */
BiChipStackSnifs::~BiChipStackSnifs(){
  vector<BiChipSnifs*>::iterator iter;
  for (iter = fBiChips.begin();iter != fBiChips.end();++iter){
    delete *iter;
  }
}


/* ===== method ======================================= */

/* ----- PreprocessBias ---------------------------------------- */
BiChipStackSnifs* BiChipStackSnifs::PreprocessBias(CatOrFile * Out, int NLines) {

  char fileName[lg_name+1];
  unsigned int iFile=0,nlines;
  BiChipStackSnifs* outStack = new BiChipStackSnifs(NLines);
  while(Out->NextFile(fileName) && iFile < fBiChips.size()) {
    // load all the file for fast processing
    if (fBiChips[iFile]->Chip(0)->Ny() > fBiChips[iFile]->Chip(1)->Ny())
      nlines = fBiChips[iFile]->Chip(0)->Ny();
    else
      nlines = fBiChips[iFile]->Chip(1)->Ny();
    BiChipSnifs * out = new BiChipSnifs(*fBiChips[iFile],fileName,FLOAT,1,kIoSlice,nlines);
    // free some memory
    fBiChips[iFile]->SetNLines(1);
    out->PreprocessBias();
    outStack->AddBiChip(out);
    iFile++;
  }
  return outStack;
  
}

/* ----- PreprocessDark ---------------------------------------- */
ImageStackSnifs* BiChipStackSnifs::PreprocessDark(CatOrFile * Out, BiChipSnifs* bias, int NLines) {

  char fileName[lg_name+1];
  unsigned int iFile=0,nlines;
  ImageStackSnifs* outStack = new ImageStackSnifs(NLines);
  while(Out->NextFile(fileName) && iFile < fBiChips.size()) {
    // load all the file for fast processing
    if (fBiChips[iFile]->Chip(0)->Ny() > fBiChips[iFile]->Chip(1)->Ny())
      nlines = fBiChips[iFile]->Chip(0)->Ny();
    else
      nlines = fBiChips[iFile]->Chip(1)->Ny();
    BiChipSnifs * tmpBi = new BiChipSnifs(*fBiChips[iFile],"mem://tmp.fits",FLOAT,1,kIoPlain);
    // free some memory
    fBiChips[iFile]->SetNLines(1);
        ImageSnifs * out = tmpBi->PreprocessAssemble(fileName,bias);
    delete tmpBi;
    out->Assembled2Dark();
    outStack->AddImage(out);
    iFile++;
  }
  return outStack;
  
}



/* ----- Kombine ---------------------------------------- */
BiChipSnifs* BiChipStackSnifs::Kombine(char* OutName,  Kombinator * K) {

  BiChipSnifs * out = new BiChipSnifs(*fBiChips[0],OutName,FLOAT,0,kIoSlice,fNLinesMem);
  vector<ImageSimple * > imList;
  for (int iChip=0;iChip<2;iChip++) {
    for (unsigned int i=0;i<fBiChips.size();i++) {
      imList.push_back(fBiChips[i]->Chip(iChip));
      
    }
    if (!out->Chip(iChip)->Variance()) {
      out->Chip(iChip)->CreateVarianceFrame();
    }
    
    ImageStack images(imList);
    images.SetKombinator(K);
    images.Kombine(out->Chip(iChip));
    imList.clear();
  }
  return out;
}


/* ----- KombineFit ---------------------------------------- */
void BiChipStackSnifs::KombineFit(BiChipSnifs ** Out, char** OutName, KombinatorFit * K, ValueGetter * V) {

  ImageSimple ** outSimple;
  outSimple = new ImageSimple*[K->NParam()];
  for (int iOut=0;iOut<K->NParam();iOut++)
    Out[iOut] = new BiChipSnifs(*fBiChips[0],OutName[iOut],FLOAT,0,kIoSlice,fNLinesMem);
    
  vector<ImageSimple * > imList;
  for (int iChip=0;iChip<2;iChip++) {
    for (unsigned int i=0;i<fBiChips.size();i++) {
      imList.push_back(fBiChips[i]->Chip(iChip));
    }
    for (int iOut=0;iOut<K->NParam();iOut++)
      outSimple[iOut] = Out[iOut]->Chip(iChip);
    //    if (!Out[i]->Chip(iChip)->Variance()) {
    //  Out->Chip(iChip)->CreateVarianceFrame();
    //}
    
    ImageStack images(imList);
    images.SetKombinatorFit(K);
    images.SetValueGetter(V);
    images.KombineFit(outSimple,0);
    imList.clear();
  }
  return;
}


/* ----- MakeBiasFrame ---------------------------------------- */
BiChipSnifs* BiChipStackSnifs::MakeBiasFrame(CatOrFile * tmpOut, char* BiasName, double SigCut) {
  
  BiChipStackSnifs* tmp = PreprocessBias(tmpOut,fNLinesMem);
  KGauss gauss(SigCut);
  BiChipSnifs* out = tmp->Kombine(BiasName,&gauss);
  delete tmp;
  return out;
  
}
