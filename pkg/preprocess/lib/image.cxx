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
#include "gsl/gsl_statistics.h"
//#include "gsl/gsl_fit.h"

//#include "IFU_io.h"
//#include "IFU_math.h"

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
int ImageSimple::WrDesc(char* Descr, short Type, int NbElements, void* Values) {
  // Wrapper to RdDesc
  return WR_desc(Frame(), Descr,Type,NbElements,  Values);
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

/* ----- Inside  -------------------------------------------------- */
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
  // puts the section Sec at starting pos X1Start,
  // and clips the image without warnings if it goes outside bounds
  // X/Y Dir are set to -1 to flip the image
  // ZScale will scale the image imported

  
  if (this==From) {
    print_error("ImageSimple::ImportSection, images shall be different");
    return;
  }

  double val;

  for (int iy = Sec->YFirst(); iy<Sec->YLast();iy++)
    for (int ix = Sec->XFirst(); ix<Sec->XLast();ix++)
      if ( From->Inside(ix,iy) 
           && Inside ( X1Start-1+(ix-Sec->XFirst())*XDir,
                       Y1Start-1+(iy-Sec->YFirst())*YDir)) {
 
        val = From->RdFrame(ix,iy)*ZScale;
        WrFrame( X1Start-1+(ix-Sec->XFirst())*XDir,
                 Y1Start-1+(iy-Sec->YFirst())*YDir,val);
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
  
  double sum;
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


/* ===== Methods ====================================================== */

/* ----- Overscan ----------------------------------------------------- */

#ifdef OLD
void ImageSimple::SubstractOverscan(Section* S) {
  // Substracts the overscan from whole image, taking the datas from 
  // the supplied window

  // 2 loops for in-place substrction

  // Build the list of overscan values
  double ov[S->YLength()];
  ComputeLinesOverscan(S,ov);
  
  
  // Substract the overscan - makes a ramp between left and right
  // overscans. Left overscan = right line before. The fixed point
  // for the ramp are the middle ov overscan strip
  // An small extrapolation is made on one line after the fixed point
  // As it amounts to about 14 points, it is considered as a negligible
  // approximation 
  for (int col=0;col<Ny();col++) {
    for (int line=0;line<Nx();line++) {
      double ovGuess;
      // no data before - be careful
      if (col<S->YFirst())
        ovGuess=ov[0];
      // no data after - no info
      else if (col>=S->YLast())
        ovGuess = ov[S->YLength()-1];
      // first line : extrapolate the slope (best guess)
      else if (col==S->YFirst()) {
        ovGuess = (ov[S->YFirst()+1]-ov[S->YFirst()])
          * (line- (S->X1()+S->X2())/2.0+1+Nx() ) / (Nx()*1.0) 
          + 2*ov[S->YFirst()] - ov[S->YFirst()+1];
      }
      // Ok, do it.
      else {
        ovGuess = (ov[col-S->YFirst()]-ov[col-S->YFirst()-1])
          * (line- (S->X1()+S->X2())/2.0+1+Nx() ) / (Nx()*1.0) 
          + ov[col-S->YFirst()-1];

      }
      double val = RdFrame(line,col) -ovGuess;
      WrFrame(line,col,val);   
    } // for line
  }

  // And now the effect on the variance
  if (Variance()) {
    double rms,addedvar;
    RdDesc("RDNOISE",DOUBLE,1,&rms);
    double linevar = rms*rms / S->XLength();

    for (int col=0;col<Ny();col++) {
      for (int line=0;line<Nx();line++) {
        // no data before or after : we don't know anything
        // choice is to get the first parameter
        // we could have put the variance to infinity
        if (col<=S->YFirst() || col>=S->YLast() )
          addedvar = linevar;
        // Ok, do it.
        else {
          // nlines of the 2 windows
          int nLines1 = nLinesDefault;
          if (col-S->YFirst() < nLinesDefault) 
            nLines1 = col-S->YFirst();
          if (S->YLast()-1-col < nLines1)
            nLines1 = S->YLast()-1-col;
          int nLines2 = nLinesDefault;
          if (col-1-S->YFirst() < nLinesDefault) 
            nLines2 = col-1-S->YFirst();
          if (S->YLast()-col < nLines2)
            nLines2 = S->YLast()-col;
          // the math are computed as if it were the median of
          // a square distribution (which is false, and only a good
          // approximation with big N)
          // As the lines are fully correlated, we shall 
          // add also a small something due to the interpolation
          // but the error we do here is nevertheless negligible
          // with respect to readout error.
          double addedvar1 = linevar / (2*nLines1+3) * 3 ;
          double addedvar2 = linevar / (2*nLines2+3) * 3 ;
          addedvar = (addedvar1-addedvar2)
            * (line- (S->X1()+S->X2())/2.0+1+Nx() ) / (Nx()*1.0)
            +addedvar2;  
        }
        double val = Variance()->RdFrame(line,col) + addedvar;
        Variance()->WrFrame(line,col,val);   
      } // for line
    }
  } //if variance
}

void ImageSimple::ComputeLinesOverscan(Section* S,double *ov) {
  // Computes the value of the overscan for 1 line

  // Various methods exists. 
  // The overscan frame can not be treated as a bias frame.
  // Polynomial fit was proven not to be satisfactory.
  // The sauron method was : median of the mean of the lines
  // Mean of the line is not fully adequate as we saw that the first 
  // lines of the 
  // overscan are dangerous to use for the blueCCD, As they amount for 1% of
  // last data column. But median of the line needs a removed odd-even effect
  // if present.
  // Median of n lines is perhaps adequate in case of sudden bumps. Variants 
  // of the Lee filter may perhaps also be applied. 

  // for the moment : median window of 5 lines, mean of the line.
  // the window is restricted to symetric window near the edges

  //
  // First the value line per line
  //

  double ovColRaw[S->YLength()];
  double ovLine[S->XLength()];
  int col;

  for (col=S->YFirst();col<S->YLast();col++) {

    for (int line=S->XFirst();line<S->XLast();line++)
      ovLine[line-S->XFirst()] = RdFrame(line,col);

    // implies the overscan is pure from spurious...
    // but we do not get enough data in order to know
    // what happens in the beginning of the overscan strip.
    // so the best would be the truncated mean with known variance method
    // -> which cuts ?
    ovColRaw[col-S->YFirst()] = gsl_stats_mean(ovLine,1,S->XLength());
  }

  //
  // now compute something smarter -> window
  //


  // whith this choice, this makes a theoretical error mean,median of 
  // less than
  // 1/10 of the readout error for our CCDs, which is affordable !
  // except of course at the beginning, but there the systematics
  // effects are unavoidable due to the strip structure (essentially at 
  // the beginning) 
  // const int nLinesDefault = 5;
  double buffer[nLinesDefault*2+1];
  int index[nLinesDefault*2+1];

  for (col=S->YFirst() ; col<S->YLast() ; col++) {
    int nLines = nLinesDefault;
    // we restrict the window near the edges in order to get a symetric
    // window around the value. Bad overscan guess else occurs.
    if (col-S->YFirst() < nLinesDefault) 
      nLines = col-S->YFirst();
    if (S->YLast()-1-col < nLines)
      nLines = S->YLast()-1-col;
    

    // now the loop for the median
    for (int iCol=col-nLines; iCol<=col+nLines;iCol++) {
      buffer[iCol-col+nLines] = ovColRaw[iCol];
    }
    ov[col]=median(buffer,nLines*2+1,index);
  } // for col
}

#endif

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

#ifdef OLD
/* ----- odd-Even ----------------------------------------------------- */
void ImageSimple::OddEvenCorrect(Section* S,double* param, double sigcut ) {
  //   Odd-even substraction works as follows
  //   1) compute odd-even for each line in the overscan strip
  // 2) interpolate a linear fit
  // 3) substract it
  // NOTE : if there is no odd-even effect, the total error added shall be 
  // negligible
  //
  // X and Y define the window, sigcut the level of oulier removal (in sigmas)
  

  // fill the array
  double oddEven[S->XLength()/2];
  double array[S->YLength()];
  double weight[S->YLength()];
  int i,j,remain;
  
  double sigma=0;
  if (sigcut)
    sigma = OverscanRms(S,ut_fraction_sigcut(1.0/(S->XLength()-0.5)));

  for (j= S->YFirst() ;j< S->YLast();j++) {
    /* fill the intermediate odd-even array to get rid of cosmics */
    for (i= S->XFirst() ; i<S->XLast() ; i+=2) /* 'odd' */
      oddEven[(i-S->XFirst())/2] = RdFrame(i,j)-RdFrame(i+1,j);
    
    remain = S->XLength()/2;

    double * oddEvenRet = oddEven;
    if (sigcut>0)
      ut_trunc_sigma_known(&oddEvenRet,&remain,sigma,sigcut);
    
    /* now compute the mean */
    array[j- S->YFirst() ]=0;
    for(i=0;i<remain;i++)
      array[j- S->YFirst() ] += oddEvenRet[i];
    array[j- S->YFirst() ] /= remain;  
    weight[j- S->YFirst() ] = remain;
  }

  // linear fit to the array (with wrapper to gsl_fit_wlinear)
  double X[ S->YLength() ];
  double cov00,cov01,cov11,chi2;
  for (i=0;i< S->YLength() ;i++)
    X[i]=i;
  gsl_fit_wlinear (X, 1, weight, 1, array, 1, S->YLength(), 
                   param, param+1, &cov00, &cov01, &cov11, 
                   &chi2);
  
  // perform the substraction (on all image)
  double sub,val;
  for (j=0;j<Ny();j++) {
    sub = (param[0] + (j- S->YFirst() )*param[1])/2;
    for (i=( S->XFirst() )%2;i<Nx();i+=2) {
      val = RdFrame(i,j) -sub;
      WrFrame(i,j,val);
    }
    for (i=( S->XFirst() +1)%2;i<Nx();i+=2) {
      val = RdFrame(i,j) +sub;
      WrFrame(i,j,val);
    }
  }
}

#endif

/* ----- AddPoissonNoise-------------------------------------------------- */
void ImageSimple::AddPoissonNoise() {
  
  for (int j=0;j<Ny();j++)
    for (int i=0;i<Nx();i++) {
      if (RdFrame(i,j)>0) { 
        double val = Variance()->RdFrame(i,j)+RdFrame(i,j);
        Variance()->WrFrame(i,j,val);
      }
    }
}

/* ----- AddPoissonNoise-------------------------------------------------- */
void ImageSimple::HandleSaturation(double Level) {
  for (int j=0;j<Ny();j++)
    for (int i=0;i<Nx();i++) {
      if (RdFrame(i,j)>=Level)
        Variance()->WrFrame(i,j,1.e10);
    }
}


