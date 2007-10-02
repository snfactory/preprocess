/* === Doxygen Comment ======================================= */
/*! 
 * \file          image.cxx
 * \copyright     (c) 2004 SNIFS Collaboration
 * \date          Wed Aug  6 17:44:35 2003
 * \author        Emmanuel Gangler <e.gangler@ipnl.in2p3.fr>
 * \version       0.0
 * \brief         
 *                
 * $Id$
 **/
/* =========================================================== */

/* ------ Library includes ----- */
#include "gsl/gsl_statistics.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>

/* ----- Local includes ----- */
#include "utils.h"
#include "image.hxx"
#include "section.hxx"
#include "ioplain.hxx"
#include "ioslice.hxx"

/* ##### IMAGE SIMPLE ################################################# */

/* ===== constructor/destructor ======================================= */

/* ----- ImageSimple -------------------------------------------------- */
ImageSimple::ImageSimple(IoMethod_t Method, int MParam) {

  SetIoMethod(Method,MParam);
  fVariance=0;
}

/* ----- ImageSimple -------------------------------------------------- */
ImageSimple::ImageSimple(char* name, char* mode,IoMethod_t Method, int MParam) {
  SetIoMethod(Method,MParam);
  OpenFrame(name,mode);
  fVariance=0;
}


/* ----- ImageSimple -------------------------------------------------- */
ImageSimple::ImageSimple(const ImageSimple &image,char* newname,short newtype,int copydata, IoMethod_t Method, int MParam) {
  // Copy constructor
  SetIoMethod(Method,MParam);
  fVariance =0;
  int npix[2];
  double start[2];
  double step[2];

  npix[0] = image.Nx();
  npix[1] = image.Ny();
  start[0] = image.Frame()->startx;
  start[1] = image.Frame()->starty;
  step[0] = image.Frame()->stepx;
  step[1] = image.Frame()->stepy;
  newtype = (newtype==0) ? image.Frame()->data_type : newtype;

  CreateFrameFull(newname,npix,start,step,newtype,image.Frame()->ident,image.Frame()->cunit);

  CP_non_std_desc(image.Frame(),Frame());

  if (copydata)
    for (int line=0 ; line<image.Ny(); line ++)
      for (int col=0 ; col<image.Nx(); col ++)
        WrFrame(col,line,image.RdFrame(col,line));
}  


/* ----- ~ImageSimple ------------------------------------------------- */
ImageSimple::~ImageSimple(){
  if (fVariance)
    delete fVariance;
  if (strstr(Name(),"mem://"))
    DeleteFrame();
  else
    CloseFrame();
  delete fIo;
}

/* ===== Wrappers to IMAGE2D methods ================================== */

/* ----- RdDesc ------------------------------------------------------- */
int ImageSimple::RdDesc(char* Descr, short Type, int NbElements, void* Values) const {
  // Wrapper to RdDesc
  return RD_desc(Frame(), Descr,Type,NbElements,  Values); 
}

/* ----- WrDesc ------------------------------------------------------- */
int ImageSimple::WrDesc(char* Descr, short Type, int NbElements, const void* Values) {
  // Wrapper to WrDesc
  // necessary because C library doesn't know about constness
  void * values = const_cast<void*>(Values);
  return WR_desc(Frame(), Descr,Type,NbElements,  values);
}

/* ----- RdIfDesc ----------------------------------------------------- */
int ImageSimple::RdIfDesc(char* Descr, short Type, int NbElements, void* Values) const {
  // Variant of RdDesc which does nothing if the descriptor does not exists
  // the content of Values in that case undefined
  disable_user_warnings();
  int retCode = RD_desc(Frame(), Descr,Type,NbElements,  Values);
  restore_user_warnings();
  return retCode;
}

/* ----- DeleteDesc ----------------------------------------------------- */
int ImageSimple::DeleteDesc(char* Descr) {
  // Deletes teh descriptor
  return delete_desc(Frame(),Descr);
}

/* ===== Setter ================================================== */
void ImageSimple::SetIoMethod(IoMethod_t Method, int MParam) { 
  if (Method == kIoSlice)
    fIo = new IoSlice(MParam);
  else
    fIo = new IoPlain();
}

void ImageSimple::SetNLines(int NLines) {
  Io()->SetMaxLines(NLines);
  if (Variance()) {
    Variance()->Io()->SetMaxLines(NLines);
  }
}

/* ===== Utility Methods ================================================== */



/* ----- Inside  -------------------------------------------------- */
int ImageSimple::Inside(int Xc,int Yc) const {
  // 1 if xc,yc is inside the frame,
  // 0 else
  // BEWARE : x and y are C-style
  if (Xc>=0 && Xc<Nx() && Yc>=0 && Yc<Ny())
    return 1;
  else
    return 0;
}

/* ----- MinMax  -------------------------------------------------- */
void ImageSimple::MinMax(Section* Sec, double * min, double * max) {
  // sets the minimum and maximum value of the frame

  *min=ut_big_value;
  *max=-ut_big_value;
  for (int iy = Sec->YFirst(); iy<Sec->YLast();iy++)
    for (int ix = Sec->XFirst(); ix<Sec->XLast();ix++) {
      if (RdFrame(ix,iy)<*min)
        *min = RdFrame(ix,iy);
      if (RdFrame(ix,iy)>*max)
        *max = RdFrame(ix,iy);
    }
}



/* ----- ImportSection  -------------------------------------------------- */
void ImageSimple::ImportSectionFrame(ImageSimple * From, Section* Sec, int X1Start, int Y1Start,int XDir,int YDir,double ZScale) {
  // puts the section Sec at starting pos X1Start (1=first),
  // and clips the image without warnings if it goes outside bounds
  // X/Y Dir are set to -1 to flip the image
  // ZScale will scale the image imported

  
  if (this==From) {
    print_error("ImageSimple::ImportSection, images shall be different");
    return;
  }

  double val;
  /* OLD 
  for (int iy = Sec->YFirst(); iy<Sec->YLast();iy++)
    for (int ix = Sec->XFirst(); ix<Sec->XLast();ix++)
      if ( From->Inside(ix,iy) 
           && Inside ( X1Start-1+(ix-Sec->XFirst())*XDir,
                       Y1Start-1+(iy-Sec->YFirst())*YDir)) {
 
        val = From->RdFrame(ix,iy)*ZScale;
        WrFrame( X1Start-1+(ix-Sec->XFirst())*XDir,
                 Y1Start-1+(iy-Sec->YFirst())*YDir,val);
      }
  */
  int iylow = Sec->YFirst();
  if (iylow <0) 
    iylow=0;
  int iyup = Sec->YLast();
  if (iyup > From->Ny())
    iyup = From->Ny();
  if (YDir==1) {
    if ((Y1Start-1+(iylow-Sec->YFirst())) <0)
      iylow = -Y1Start+1+Sec->YFirst();
    if ((Y1Start-1+(iyup-Sec->YFirst())) >Ny() )
      iyup = Ny()-Y1Start+1+Sec->YFirst();
  } else {
    if ((Y1Start-1-(iylow-Sec->YFirst())) >=Ny())
      iylow = Y1Start + Sec->YFirst() - Ny();
    if ((Y1Start-1-(iyup-Sec->YFirst())) <-1 )
      iyup = Y1Start + Sec->YFirst();
  }
  int ixlow = Sec->XFirst();
  if (ixlow<0)
    ixlow = 0;
  int ixup = Sec->XLast();
  if (ixup>From->Nx())
    ixup = From->Nx();
  if (XDir==1) {
    if ((X1Start-1+(ixlow-Sec->XFirst())) <0)
      ixlow = -X1Start+1+Sec->XFirst();
    if ((X1Start-1+(ixup-Sec->XFirst())) >Nx() )
      ixup = Nx()-X1Start+1+Sec->XFirst();
  } else {
    if ((X1Start-1-(ixlow-Sec->XFirst())) >=Nx())
      ixlow = X1Start + Sec->XFirst() - Nx();
    if ((X1Start-1-(ixup-Sec->XFirst())) <-1 )
      ixup = X1Start + Sec->XFirst();
  }

  // inlining
  int SecYFirst = Sec->YFirst();
  int SecXFirst = Sec->XFirst();

  for (int iy = iylow; iy<iyup;iy++)
    for (int ix = ixlow; ix<ixup;ix++) {
      //      if ( From->Inside(ix,iy) 
      //     && Inside ( X1Start-1+(ix-Sec->XFirst())*XDir,
      //                 Y1Start-1+(iy-SecYFirst)*YDir)) {
 
        val = From->RdFrame(ix,iy)*ZScale;
        WrFrame( X1Start-1+(ix-SecXFirst)*XDir,
                 Y1Start-1+(iy-SecYFirst)*YDir,val);
    }

}

/* ----- ImportSection  -------------------------------------------------- */
void ImageSimple::ImportSection(ImageSimple * From, Section* Sec, int X1Start, int Y1Start,int XDir,int YDir,double ZScale) {
  // pilots the ImportSectionFrame for both the image and the variance 
  // if it exists
  ImportSectionFrame(From,Sec,X1Start,Y1Start,XDir,YDir,ZScale);
  if (Variance()&&From->Variance()) {
    Variance()->ImportSection(From->Variance(),Sec,X1Start,Y1Start,XDir,YDir,ZScale*ZScale);
    }
}

/* ----- ImportFlip  -------------------------------------------------- */
void ImageSimple::ImportFlip(ImageSimple * From,int XDir,int YDir,double ZScale) {
  // pilots the ImportSection for a simple image flip
  Section sec(1,From->Nx(),1,From->Ny());
  int xstart=1,ystart=1;
  if (XDir == -1)
    xstart = From->Nx();
  if (YDir == -1)
    ystart = From->Ny();
  ImportSection(From,&sec,xstart,ystart,XDir,YDir,ZScale);
}


/* ----- Add ----------------------------------------------------- */
void ImageSimple::Add(ImageSimple* ToAdd, double Scale){
  //Adds ToAdd, scaled by Scale 

  for (int j=0;j<Ny();j++)
    for (int i=0;i<Nx();i++) {
      double val = RdFrame(i,j)+ToAdd->RdFrame(i,j)*Scale;
      WrFrame(i,j,val);
    }
  if (Variance()&&ToAdd->Variance()) {
    for (int j=0;j<Ny();j++)
      for (int i=0;i<Nx();i++) {
        double val = Variance()->RdFrame(i,j)+ToAdd->Variance()->RdFrame(i,j)*Scale*Scale;
        Variance()->WrFrame(i,j,val);
      }
  }  
}

/* ----- Add ----------------------------------------------------- */
void ImageSimple::Add(double Constant){
  //Adds Constant to all image

  for (int j=0;j<Ny();j++)
    for (int i=0;i<Nx();i++) {
      double val = RdFrame(i,j)+Constant;
      WrFrame(i,j,val);
    }
}

/* ----- Scale ----------------------------------------------------- */
void ImageSimple::Scale(double Scale){
  //Scales the image by Scale

  for (int j=0;j<Ny();j++)
    for (int i=0;i<Nx();i++) {
      double val = RdFrame(i,j)*Scale;
      WrFrame(i,j,val);
    }
  if (Variance()) {
    for (int j=0;j<Ny();j++)
      for (int i=0;i<Nx();i++) {
        double val = Variance()->RdFrame(i,j)*Scale*Scale;
        Variance()->WrFrame(i,j,val);
      }
  }

}

/* ----- Set to ----------------------------------------------------- */
void ImageSimple::SetTo(double Value){
  //Sets the content to

  for (int j=0;j<Ny();j++)
    for (int i=0;i<Nx();i++) {
      WrFrame(i,j,Value);
    }
  if (Variance()) {
    for (int j=0;j<Ny();j++)
      for (int i=0;i<Nx();i++) {
        Variance()->WrFrame(i,j,0.);
      }
  }
}

/* ----- CutLow ----------------------------------------------------- */
void ImageSimple::CutLow(double LowBound){
  // Cuts the image by a lower bound
  // The variance after the cut is questionnable !

  for (int j=0;j<Ny();j++)
    for (int i=0;i<Nx();i++) {
      if (RdFrame(i,j)<LowBound) {
          if (Variance()) {
          // We have a systematic problem !
          // Choice is to do add something systematic :
          // Due to the bounding, statistics are anyway screwed-up

          double outside = (RdFrame(i,j)-LowBound);
          double newvar = Variance()->RdFrame(i,j) + outside*outside;
          Variance()->WrFrame(i,j,newvar );
          }
        WrFrame(i,j,LowBound);
      } 
    }
}

/* ----- CutLow ----------------------------------------------------- */
int ImageSimple::AbsThreshold(double Threshold){
  // Cuts the lowest significant image values
  // The variance is not touched after the cut
  int nout=0;
  for (int j=0;j<Ny();j++)
    for (int i=0;i<Nx();i++) {
      if (fabs(RdFrame(i,j))<Threshold) {
        WrFrame(i,j,0);
        nout++;
      } 
    }
  return nout;
}

/* ----- Divide ----------------------------------------------------- */
void ImageSimple::Divide(ImageSimple* Denom){
  //Divides by Denom 

  if (Variance()&&Denom->Variance()) {
    for (int j=0;j<Ny();j++)
      for (int i=0;i<Nx();i++) {
        double num=RdFrame(i,j);
        double den=Denom->RdFrame(i,j);
        double Vnum=Variance()->RdFrame(i,j);
        double Vden=Denom->Variance()->RdFrame(i,j);
        double val = Vnum/den/den + num*num/den/den/den/den * Vden;
        Variance()->WrFrame(i,j,val);
      }
  }  
  for (int j=0;j<Ny();j++)
    for (int i=0;i<Nx();i++) {
      double val = RdFrame(i,j)/Denom->RdFrame(i,j);
      WrFrame(i,j,val);
    }

}

/* ----- MeanValue ----------------------------------------------------- */
double ImageSimple::MeanValue(Section* Sec,int step) {
  //Computes the mean value of the section data
  // the step factor is here to speed-up the computation
  
  double sum=0;
  int i,j,n=0;

  i = Sec->XFirst();
  for ( j=Sec->YFirst() ; j<Sec->YLast() ; ) {
    for ( ; i<Sec->XLast() ; i+=step ) {
      sum += RdFrame(i,j);
      n++;
    }
    while(i>=Sec->XLast()) {
      i-=Sec->XLength();
      j++;
    }
  }
  return sum/n;
}
/* ----- GetSignificance ----------------------------------------------------- */
double ImageSimple::GetSignificance(int step) {
  //Computes the mean value of the section data
  // the step factor is here to speed-up the computation
  // it assumes a gaussian distribution an returns the significance of the mean.
  
  if (!Variance()) {
    print_error("ImageSimple::GetSignificance no associated variance");
    return 0;
  }
  

  double sumX=0;
  double sumW=0;
  int i,j,n=0;

  i = 0;
  for ( j=0 ; j<Ny() ; ) {
    for ( ; i<Nx() ; i+=step ) {
      sumX += fabs(RdFrame(i,j))/sqrt(Variance()->RdFrame(i,j));
      sumW += 1;
      n++;
    }
    while(i>=Nx()) {
      i-=Nx();
      j++;
    }
  }
  return sumX/sumW;
}


/* ===== Methods ====================================================== */


/* ----- SectionRms --------------------------------------------------- */
double ImageSimple::SectionRms(Section* S,double sigcut) {
  // Computes the RMS of readout (i.e. readout noise if everything
  // goes fine) from the overscan strip.

  double lineVal[S->XLength()],*line;
  int remain;
  double varinc=0;
  double weight=0;

  for (int j=S->YFirst();j<S->YLast();j++) {

    for (int i=S->XFirst();i<S->XLast();i++)
      lineVal[i-S->XFirst()] = RdFrame(i,j);
    remain = S->XLength();
    line = lineVal;
    if (sigcut>0) {
      ut_trunc_sigma_unknown(&line,&remain,sigcut);
    }
    if (remain>1) {
      // carefull statistical analysis gives N as the weight
      // of estimated variance from sample
      // even if in the estimation, there is a 1/(N-1) factor
      varinc += remain * gsl_stats_variance(line,1,remain);
      weight += remain;
    }  
  }
  return sqrt(varinc/weight);
}

/* ----- FitBy -------------------------------------------------- */
gsl_vector* ImageSimple::FitBy(vector<ImageSimple *> Refs, double cut) {

  if (!Variance())
    print_error("ImageSimple::FitBy no variance in frame");

  // allocate some space
  gsl_matrix * XX = gsl_matrix_alloc(Refs.size(),Refs.size());
  gsl_vector * XY = gsl_vector_alloc(Refs.size());
  gsl_vector * A = gsl_vector_alloc(Refs.size());
  // Y=XA
  vector<ImageSimple*>::iterator iter, iter2;

  int nbad;
  int niter,niter2;
  do {
    gsl_matrix_set_zero(XX);
    gsl_vector_set_zero(XY);
    
    for (iter = Refs.begin(),niter=0;iter != Refs.end();++iter,niter++){
      for (int j=0;j<Ny();j++) {
        for (int i=0;i<Nx();i++) {
          double var = Variance()->RdFrame(i,j);
          if (var>ut_big_value/10)
            var=ut_big_value;
          else
            var=1;
          (*gsl_vector_ptr(XY,niter))+=(*iter)->RdFrame(i,j)*RdFrame(i,j)/var;
        }
      }
      for (iter2 = iter,niter2=niter;iter2 != Refs.end();++iter2,niter2++){
        for (int j=0;j<Ny();j++) {
          for (int i=0;i<Nx();i++) {
            double var = Variance()->RdFrame(i,j);
            if (var>ut_big_value/10)
              var=ut_big_value;
            else
              var=1;
            
            (*gsl_matrix_ptr(XX,niter,niter2))+=(*iter)->RdFrame(i,j)*(*iter2)->RdFrame(i,j)/var;      
          }
          
        }
      }
    }
    for (unsigned int n1=0;n1<Refs.size();n1++) {
      for (unsigned int n2=n1+1;n2<Refs.size();n2++) {
        gsl_matrix_set(XX,n2,n1,gsl_matrix_get(XX,n1,n2));
      } 
    }
    
    gsl_linalg_cholesky_decomp (XX);
    gsl_linalg_cholesky_solve (XX,XY,A);
    
    for (unsigned int n=0;n<Refs.size();n++)
      print_msg("Coefficient %f",gsl_vector_get(A,n));
    // removing far
    nbad=0;
    for (int j=0;j<Ny();j++) {
      for (int i=0;i<Nx();i++) {
        double variance = Variance()->RdFrame(i,j);
        double delta=RdFrame(i,j);
        for (iter = Refs.begin(),niter=0;iter != Refs.end();++iter,niter++)
          delta -= gsl_vector_get(A,niter)*(*iter)->RdFrame(i,j);
        if (delta > cut*sqrt(variance)) {
          Variance()->WrFrame(i,j,ut_big_value);
          nbad++;
        }
      }
    }
    print_msg("removed %d pixels",nbad);
  }
  while(nbad>0);
  gsl_matrix_free(XX);
  gsl_vector_free(XY);
  
  return A;

}



/* ----- AddPoissonNoise-------------------------------------------------- */
void ImageSimple::AddPoissonNoise() {
  
  double poisNois;
  for (int j=0;j<Ny();j++)
    for (int i=0;i<Nx();i++) {
      if ((poisNois = RdFrame(i,j))>0) { 
        double val = Variance()->RdFrame(i,j)+poisNois;
        Variance()->WrFrame(i,j,val);
      }
    }
}

/* ----- HandleSaturation -------------------------------------------------- */
int ImageSimple::HandleSaturation(double Level) {
  // there is an additional trick, in order to remove ptential bleeding 
  int nsat=0;
  for (int j=0;j<Ny();j++)
    for (int i=0;i<Nx();i++) {
      if (RdFrame(i,j)>=Level) {
        nsat++;
        Variance()->WrFrame(i,j,ut_big_value);
        if (j>0) Variance()->WrFrame(i,j-1,ut_big_value);
        if (j+1<Ny()) Variance()->WrFrame(i,j+1,ut_big_value);
      }
    }
  return nsat;
}


/* ----- CleanWith -------------------------------------------------- */
int ImageSimple::CleanWith(ImageSimple * Ref, double SigCut) {
  // this method replaces pixels incompatible with the reference image
  // with the reference image - variance set to infinity.
  // 
  // The main purpose is to remove cosmics
  if (!Variance() || !Ref->Variance())
    print_error("ImageSimple::CleanWith method needs variances");
  int nbad=0;
  for (int j=0;j<Ny();j++)
    for (int i=0;i<Nx();i++) {
      if (fabs(RdFrame(i,j)-Ref->RdFrame(i,j)) / sqrt(Variance()->RdFrame(i,j)+Ref->Variance()->RdFrame(i,j)) > SigCut) {
        nbad++;
        WrFrame(i,j,Ref->RdFrame(i,j));
        Variance()->WrFrame(i,j,ut_big_value);
      }
    }
  return nbad;
}


