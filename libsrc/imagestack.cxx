/* === Doxygen Comment ======================================= */
/*! 
 * \file          imagestack.cxx
 * \copyright     (c) 2003 IPNL
 * \date          Wed Aug  6 17:44:35 2003
 * \author        Emmanuel Gangler <e.gangler@ipnl.in2p3.fr>
 * \version       0.0
 * \brief         
 *                
 * $Id$
 **/
/* =========================================================== */

#include "catorfile.hxx"
#include "image.hxx"
#include "kombinator.hxx"
#include "kombinatorfit.hxx"
#include "imagestack.hxx"


/* #####  ImageStack ################################################# */


/* ===== constructor/destructor ======================================= */

/* ----- ImageStack -------------------------------------------------- */
ImageStack::ImageStack() {
  fKombinatorFit=0;
  fKombinator=0;
  fIsOwner=0;
}

/* ----- ImageStack -------------------------------------------------- */
ImageStack::ImageStack(vector<ImageSimple*> Images) {
  /* check the list */
  fKombinator=0;
  fIsOwner=0;
  fImageList = Images;
  vector<ImageSimple*>::const_iterator iter;
  iter = fImageList.begin();
  fNx = (*iter)->Nx();
  fNy = (*iter)->Ny();
  
  for (;iter !=fImageList.end();++iter) {
    if ((*iter)->Nx()!=fNx || (*iter)->Ny()!=fNy) {
      print_error("ImageStack::ImageStack : images are not identical\n");
      fImageList.clear();
    }
  }
}

/* ----- ImageStack -------------------------------------------------- */
ImageStack::ImageStack(CatOrFile *Cat, char * Mode){
  // The files are open in slice mode with 1 line buffer : 
  // the buffering is negligible with respect to the ios
  char fileName[lg_name+1];
  while(Cat->NextFile(fileName)) {
    ImageSimple * image = new ImageSimple(fileName,Mode,kIoSlice,1);
    fImageList.push_back(image);
  }
  fIsOwner=1;
  vector<ImageSimple*>::const_iterator iter = fImageList.begin();
  fNx = (*iter)->Nx();
  fNy = (*iter)->Ny();
  
  for (;iter !=fImageList.end();++iter) {
    if ((*iter)->Nx()!=fNx || (*iter)->Ny()!=fNy) {
      print_error("ImageStack::ImageStack : images are not identical\n");
      fImageList.clear();
    }
  }
}


/* ----- ImageStack -------------------------------------------------- */
ImageStack::~ImageStack() {
  vector<ImageSimple*>::iterator iter;
  if (fIsOwner) {
    for (iter = fImageList.begin();iter !=fImageList.end();iter++) {
      delete *iter;
      *iter =0;
    }
  }
}


/* ===== method ======================================= */

/* ----- Combine  -------------------------------------------------- */
void ImageStack::Kombine(ImageSimple *ToFill, int FillsVarOut, int UpdateInitialVar) {
  vector<ImageSimple*>::const_iterator iter;
  vector<double> vals;
  vector<double> vars;
  double retVar, retVal;
  vals.assign(fImageList.size(),0);
  vars.assign(fImageList.size(),0);

  // check it is OK
  if (GetKombinator()->NeedsVarIn()) {
    for (iter = fImageList.begin();iter !=fImageList.end();iter++) {
      if (!(*iter)->Variance()) {
        print_error(" ImageStack::Kombine %s needs a variance",(*iter)->Name() );
      }
    }
  } else {
    UpdateInitialVar=0;
  }
  // check images are all of the same size
  for (unsigned int n=1;n<fImageList.size();n++) {
    if (fImageList[n]->Nx() != Nx() ||  fImageList[n]->Ny() != Ny())
      print_error("ImageStack::Kombine %s and %s not of the same size",fImageList[0]->Name(),fImageList[n]->Name());
  }

  if (FillsVarOut && ! ToFill->Variance())
    print_error(" ImageStack::Kombine %s needs a variance frame",ToFill->Name() );

  reset_print_progress();
  for (int j=0;j<Ny();j++) {
    if (VERBOSE)
      print_progress("Kombining",(j+1.0)*100.0/Ny(),1.0);  
    for (int i=0;i<Nx();i++) {
      for (unsigned int n=0;n<fImageList.size();n++) {
        vals[n] = fImageList[n]->RdFrame(i,j);
        if (GetKombinator()->NeedsVarIn())
          vars[n] = fImageList[n]->Variance()->RdFrame(i,j);
          }
      GetKombinator()->Kombine(&vals,&vars,&retVal,&retVar);
      ToFill->WrFrame(i,j,retVal);
      if (FillsVarOut) {
        ToFill->Variance()->WrFrame(i,j,retVar);
      }
      if (UpdateInitialVar) {
        for (unsigned int n=0;n<fImageList.size();n++) {
          fImageList[n]->Variance()->WrFrame(i,j,vars[n]);
        }
      }
    }
  }
}

/* ----- Combine  -------------------------------------------------- */
void ImageStack::KombineFit(ImageSimple **ToFill, int FillsVarOut, int UpdateInitialVar) {
  vector<ImageSimple*>::const_iterator iter;
  vector<double> x;
  vector<double> vals;
  vector<double> vars;
  double retVar, *retVal;
  int nParam = GetKombinatorFit()->NParam();
  retVal = new double[nParam];
  vals.assign(fImageList.size(),0);
  vars.assign(fImageList.size(),0);
  x.assign(fImageList.size(),0);

  // check it is OK
  if (GetKombinatorFit()->NeedsVarIn()) {
    for (iter = fImageList.begin();iter !=fImageList.end();iter++) {
      if (!(*iter)->Variance()) {
        print_error(" ImageStack::KombineFit %s needs a variance",(*iter)->Name() );
      }
    }
  } else {
    UpdateInitialVar=0;
  }
  
  // Covariance matrix is not yet implemented
  // because I do not know how to do !
  //  if (FillsVarOut && GetKombinatorFit()->FillsVarOut() && ! ToFill->Variance())
  //  print_error(" ImageStack::KombineFit %s needs a variance frame",ToFill->Name() );
  
  for (unsigned int n=0;n<fImageList.size();n++) {
    x[n]= GetValue(fImageList[n]);
  }

  for (int j=0;j<Ny();j++) {
    for (int i=0;i<Nx();i++) {
      for (unsigned int n=0;n<fImageList.size();n++) {
        vals[n] = fImageList[n]->RdFrame(i,j);
        if (GetKombinatorFit()->NeedsVarIn())
          vars[n] = fImageList[n]->Variance()->RdFrame(i,j);
      }
      GetKombinatorFit()->KombineFit(&x,&vals,&vars,retVal,&retVar);
      for (int nout=0;nout<nParam;nout++)
        ToFill[nout]->WrFrame(i,j,retVal[nout]);
      //      if (FillsVarOut) {
      //ToFill->Variance()->WrFrame(i,j,retVar);
      
      if (UpdateInitialVar) {
        for (unsigned int n=0;n<fImageList.size();n++) {
          fImageList[n]->Variance()->WrFrame(i,j,vars[n]);
        }
      }
    }
  }
}

/* ----- KombineFitND  -------------------------------------------------- */
void ImageStack::KombineFitND(ImageStack *outImages, ImageStack *outvarImages) {
  vector<ImageSimple*>::const_iterator iter;
  gsl_vector* vals = gsl_vector_alloc(fImageList.size());
  gsl_vector* vars = gsl_vector_alloc(fImageList.size());
  int nparams=GetKombinatorFitND()->NParam();
  gsl_matrix* X = gsl_matrix_alloc(fImageList.size(),nparams);
  gsl_vector* Xvals = gsl_vector_alloc(nparams);
  gsl_vector* retVal = gsl_vector_alloc(nparams);
  gsl_matrix* retCovar = gsl_matrix_alloc(nparams,nparams);

  // check it is OK
  if (GetKombinatorFitND()->NeedsVarIn()) {
    for (iter = fImageList.begin();iter !=fImageList.end();iter++) {
      if (!(*iter)->Variance()) {
        print_error(" ImageStack::Kombine %s needs a variance",(*iter)->Name() );
      }
    }
  } 

  // check images are all of the same size
  for (unsigned int n=1;n<fImageList.size();n++) {
    if (fImageList[n]->Nx() != Nx() ||  fImageList[n]->Ny() != Ny())
      print_error("ImageStack::Kombine %s and %s not of the same size",fImageList[0]->Name(),fImageList[n]->Name());
  }

  // assumes the value doesn't depend on i,j ...
  for (unsigned int n=0;n<fImageList.size();n++) {
    GetValues(fImageList[n],Xvals);
    gsl_matrix_set_row(X,n,Xvals);
  }

  reset_print_progress();
  GetKombinatorFitND()->WorkspaceAlloc(fImageList.size());
  for (int j=0;j<Ny();j++) {
    if (VERBOSE)
      print_progress("Kombining",(j+1.0)*100.0/Ny(),1.0);  
    for (int i=0;i<Nx();i++) {
      for (unsigned int n=0;n<fImageList.size();n++) {
	gsl_vector_set(vals,n,fImageList[n]->RdFrame(i,j));
        if (GetKombinatorFitND()->NeedsVarIn()) {
          gsl_vector_set(vars,n,fImageList[n]->Variance()->RdFrame(i,j));
	}
      }

      GetKombinatorFitND()->KombineFit( X , vals,vars, retVal, retCovar);

      for (int n=0;n<nparams;n++)
	(*outImages->List())[n]->WrFrame(i,j,gsl_vector_get(retVal,n));
      int count=0;
      for (int ni=0;ni<nparams;ni++) {
	for (int nj=ni;nj<nparams;nj++) {
	  (*outvarImages->List())[count]->WrFrame(i,j,gsl_matrix_get(retCovar,ni,nj));
	  count++;
	}
      }
    }
  }
  gsl_vector_free(vals);
  gsl_vector_free(vars);
  gsl_matrix_free(X);
  gsl_vector_free(Xvals);
  gsl_vector_free(retVal);
  gsl_matrix_free(retCovar);
}
