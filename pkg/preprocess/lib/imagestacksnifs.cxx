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

#include "gsl/gsl_linalg.h"
#include "gsl/gsl_blas.h"
#include "gsl/gsl_eigen.h"

#include "utils.h"
#include "catorfile.hxx"
#include "imagesnifs.hxx"
#include "bichip.hxx"
#include "imagestack.hxx"
#include "kombinator.hxx"
#include "kombinatorfit.hxx"
#include "imagestacksnifs.hxx"


/* ===== constructor/destructor ======================================= */

/* ----- imagestacksnifs  -------------------------------------------------- */
ImageStackSnifs::ImageStackSnifs(CatOrFile* Cat, char * Mode,  int NLines){
  /* first loads the headers - to check the homogeneity of the catalog */
  /* the ioslice is used, this means the numer of loaded lines has to be 
   reworked for image-based processing */
  fNLinesMem = NLines;

  char fileName[lg_name+1];
  if (!Cat->NextFile(fileName)) {
    print_error("ImageStackSnifs::ImageStackSnifs : unable to read %s",Cat->Name());
    return;
  } 

  int nDone=0;
  ImageSnifs * imRef=0;
  do {
    if (ut_is_bichip_detcom(fileName)) {
      print_error("ImageStackSnifs::ImageStackSnifs can not read bichips %s",fileName);
        return; 
      }

    // loads the image
    ImageSnifs * image = new ImageSnifs(fileName,Mode,kIoSlice,fNLinesMem);
    fImages.push_back(image);
    if (!nDone) {
      imRef = image;
    } else 
      if (!imRef->CanBeStackedWith(image)) {
        print_error("ImageStackSnifs::ImageStackSnifs %s is incompatible with %s",image->Name(),imRef->Name());
        return;
      }
    nDone++;  
  } while(Cat->NextFile(fileName));
  
  
}


/* ----- ImageStackSnifs  -------------------------------------------------- */
ImageStackSnifs::ImageStackSnifs(char* Name, char * Mode, int NLines){
  /* open as a stack a frame made of imageNNN varianceNNN */

  fNLinesMem = NLines; 

  char hduName[lg_name+1];
  ImageSnifs * out;
  int i;
  for (i=0 ; ; i++) {
    sprintf(hduName,"%s[image%03d]",Name,i);
    if (exist(hduName)) {
      out = new ImageSnifs(hduName,Mode,kIoSlice,fNLinesMem);
      fImages.push_back(out);
    }
    else
      break;
  }
  if (i == 0 ) {
    print_error("ImageStackSnifs::ImageStackSnifs : expect at least an image000 extension in %s", Name);
    return;
  }
  int nImages=i; // the really relevant number
  

  for (i=0 ; ; i++) {
    sprintf(hduName,"%s[variance%03d]",Name,i);
    if (exist(hduName)) {
      out = new ImageSnifs(hduName,Mode,kIoSlice,fNLinesMem);
      fImages.push_back(out);
    }
    else
      break;
  }
  if (i != 0 && i != (nImages * (nImages +1) / 2)) {
    print_error("ImageStackSnifs::ImageStackSnifs : expect 0 or %d variance terms", (nImages * (nImages +1) / 2));
    return;
  }
}

/* ----- ImageStackSnifs  -------------------------------------------------- */
ImageStackSnifs::ImageStackSnifs(ImageSnifs *Image, char* Name, int NImages, int NLines){
  /* creates a stack a frame made of imageNNN varianceNNN 
   image is used as a template when crating new files*/

  fNLinesMem = NLines; 

  char hduName[lg_name+1];
  ImageSnifs * out;
  int i;
  for (i=0 ; i<NImages; i++) {
    sprintf(hduName,"%s[image%03d]",Name,i);
    out = new ImageSnifs(*Image,hduName,FLOAT,0,kIoSlice,fNLinesMem);
    fImages.push_back(out);
  }

  for (i=0 ; i< NImages*(NImages+1)/2; i++) {
    sprintf(hduName,"%s[variance%03d]",Name,i);
    out = new ImageSnifs(*Image,hduName,FLOAT,0,kIoSlice,fNLinesMem);
    fImages.push_back(out);
  }
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

/* ----- KombineFitND ---------------------------------------- */
ImageStackSnifs* ImageStackSnifs::KombineFitND(char* OutName,  KombinatorFitND * K, ValuesGetter * VG) {

#ifdef OLD
  ImageStackSnifs * outStack = new ImageStackSnifs(fNLinesMem);
  char hduOutName[lg_name+1];
  for (int i=0;i<K->NParam();i++) {
    sprintf(hduOutName,"%s[image%03d]",OutName,i);
    ImageSnifs * out = new ImageSnifs(*fImages[0],hduOutName,FLOAT,0,kIoSlice,fNLinesMem);
    outStack->AddImage(out);
  }
  int count=0;
  for (int i=0;i<K->NParam();i++) {
    for (int j=i;j<K->NParam();j++) {
      sprintf(hduOutName,"%s[variance%03d]",OutName,count);
      count++;
      ImageSnifs * out = new ImageSnifs(*fImages[0],hduOutName,FLOAT,0,kIoSlice,fNLinesMem);
      outStack->AddImage(out);
    }
  }
#endif
  ImageStackSnifs * outStack = new ImageStackSnifs(fImages[0],OutName,K->NParam(),fNLinesMem);

  vector<ImageSimple * > imList, outList, outvarList;
  for (unsigned int i=0;i<fImages.size();i++) {
    imList.push_back(fImages[i]);
  }
  for (int i=0;i<K->NParam();i++) {
    outList.push_back(outStack->GetImages()[i]);
  }
  for (int i=0;i<K->NParam()*(K->NParam()+1)/2;i++) {
    outvarList.push_back(outStack->GetImages()[i+K->NParam()]);
  }
    
  ImageStack images(imList);
  ImageStack outImages(outList);
  ImageStack outvarImages(outvarList);

  images.SetKombinatorFitND(K);
  images.SetValuesGetter(VG);
  images.KombineFitND(&outImages,&outvarImages);

  imList.clear();
  outList.clear();
  outvarList.clear();

  return outStack;
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
    if (!ut_is_bichip_detcom(fileName)) {
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

#ifdef OLD
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
#endif



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

/* #####  ImageStackStack ################################################# */

/* ===== constructor/destructor ======================================= */

/* ----- imagestackstack  -------------------------------------------------- */
ImageStackStack::ImageStackStack(CatOrFile* Cat, char * Mode,  int NLines){
  /* first loads the headers - to check the homogeneity of the catalog */
  /* the ioslice is used, this means the numer of loaded lines has to be 
   reworked for image-based processing */
  fNLinesMem = NLines;

  char fileName[lg_name+1];
  if (!Cat->NextFile(fileName)) {
    print_error("imagestacksnifs::imagestacksnifs : unable to read %s",Cat->Name());
    return;
  } 

  int nDone=0;
  //ImageSnifs * stackRef=0;
  do {
    //if (ut_is_bichip_detcom(fileName)) {
    //print_error("ImageStackSnifs::ImageStackSnifs can not read bichips %s",fileName);
    //  return; 
    //}

    // loads the stack/image
    ImageStackSnifs * stack = new ImageStackSnifs(fileName,Mode,fNLinesMem);
    fStacks.push_back(stack);
    //if (!nDone) {
    //imRef = image;
    //} else 
    //if (!imRef->CanBeStackedWith(image)) {
    //  print_error("ImageStackSnifs::ImageStackSnifs %s is incompatible with %s",image->Name(),imRef->Name());
    //  return;
    //}
    nDone++;  
  } while(Cat->NextFile(fileName));
  
  
}

/* ----- ~ImageStackStack  -------------------------------------------------- */

ImageStackStack::~ImageStackStack(){
  vector<ImageStackSnifs*>::iterator iter;
  for (iter = fStacks.begin();iter != fStacks.end();++iter){
    delete *iter;
  }
}


/* ===== method ======================================= */

/* ----- Kombine ---------------------------------------- */
ImageStackSnifs * ImageStackStack::Kombine(char* OutName) {

  unsigned int nimages= int (rint( (sqrt(9+8*fStacks[0]->GetImages().size()) - 3) /2 ));
  ImageStackSnifs* outStack = new ImageStackSnifs(fStacks[0]->GetImages()[0],OutName,nimages,fNLinesMem);

  gsl_vector * vals=gsl_vector_calloc(nimages);
  gsl_matrix * vars=gsl_matrix_calloc(nimages,nimages);
  gsl_vector * sumwx=gsl_vector_calloc(nimages);
  gsl_matrix * sumw=gsl_matrix_calloc(nimages,nimages);

  gsl_eigen_symm_workspace * workspace= gsl_eigen_symm_alloc(nimages);
  gsl_vector * eigvals=gsl_vector_calloc(nimages);
  gsl_matrix * copy=gsl_matrix_alloc(nimages,nimages);
  int *badimages = new int[fStacks.size()];

  // loop on pixels row-wise
  for (int j=0;j<fStacks[0]->GetImages()[0]->Ny();j++) {
    if (VERBOSE)
      print_progress("Kombining",(j+1.0)*100.0/fStacks[0]->GetImages()[0]->Ny(),1.0);  
    for (int i=0;i<fStacks[0]->GetImages()[0]->Nx();i++) {

      // these vectors have to be set to 0...
      gsl_vector_set_zero(sumwx);
      gsl_matrix_set_zero(sumw);

      for (unsigned int n=0;n<fStacks.size();n++) {
	int count=0;
	for (unsigned int k=0;k<nimages;k++) {
	  gsl_vector_set(vals,k,fStacks[n]->GetImages()[k]->RdFrame(i,j));
	  count++;
	}
	for (unsigned int ki=0;ki<nimages;ki++) {
	  for (unsigned int kj=ki;kj<nimages;kj++) {
	    gsl_matrix_set(vars,ki,kj,fStacks[n]->GetImages()[count]->RdFrame(i,j));
	    gsl_matrix_set(vars,kj,ki,fStacks[n]->GetImages()[count]->RdFrame(i,j));
	    count++;
	  }
	}

	// but we do want sum W and sum WX
	// Unfortunately, due to some numerical traps, vars may be non pos def !
	gsl_matrix_memcpy(copy,vars);
	gsl_eigen_symm(copy,eigvals,workspace);
	unsigned int k;
	for (k=0;k<nimages;k++) {
	  if (gsl_vector_get(eigvals,k)<= 1e-12) { // rounding error
	    if (VERBOSE and badimages[n]==0) {
	      print_warning("dropping image %d, pixel %d,%d",n,i,j);
	      badimages[n]=1;
	    }
	    break;
	  }
	}
	// ok, we can invert the matrix
	if (k==nimages) {
	  gsl_linalg_cholesky_decomp(vars);
	  gsl_linalg_cholesky_invert(vars);
	  gsl_blas_dsymv(CblasUpper, 1.0, vars, vals, 1.0, sumwx);
	  gsl_matrix_add(sumw,vars);
	}
      }	// loop on stacks
       
      // here : invert the solution and store it.
      gsl_linalg_cholesky_decomp(sumw);
      gsl_linalg_cholesky_invert(sumw);
      gsl_blas_dsymv(CblasUpper, 1.0, sumw, sumwx, 0.0, vals);
	
      // store result
      int count=0;
      for (unsigned int k=0;k<nimages;k++) {
	outStack->GetImages()[count]->WrFrame(i,j,gsl_vector_get(vals,k));
	count++;
      }
      for (unsigned int ki=0;ki<nimages;ki++) {
	for (unsigned int kj=ki;kj<nimages;kj++) {
	  outStack->GetImages()[count]->WrFrame(i,j,gsl_matrix_get(sumw,ki,kj));
	  count++;
	}
      }
    }
  } // loop on pixels

  gsl_vector_free(vals);
  gsl_vector_free(sumwx);
  gsl_matrix_free(vars);
  gsl_matrix_free(sumw);

  gsl_eigen_symm_free(workspace);
  gsl_vector_free(eigvals);
  gsl_matrix_free(copy);
  delete[] badimages;

  return outStack;
}
