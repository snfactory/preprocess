/* === Doxygen Comment ======================================= */
/*! 
 * \file          image.cxx
 * \copyright     (c) 2003 CRAL-Observatoire de Lyon
 * \date          Wed Aug  6 17:44:35 2003
 * \author        Emmanuel Gangler <e.gangler@ipnl.in2p3.fr>
 * \version       0.0
 * \brief         
 *                
 * $Id$
 **/
/* =========================================================== */
#include "gsl/gsl_statistics.h"
#include "gsl/gsl_fit.h"

#include "utils.h"

#include "image.hxx"
// for median
#include "IFU_math.h"
#include "section.hxx"

/* ##### IMAGE SIMPLE ################################################# */

const int ImageSimple::nLinesDefault = 5;


/* ===== constructor/destructor ======================================= */

/* ----- ImageSimple -------------------------------------------------- */
ImageSimple::ImageSimple() {
  fFrame = new IMAGE2D;  
  fLoaded = false;
  fVariance=0;
}

/* ----- ImageSimple -------------------------------------------------- */
ImageSimple::ImageSimple(char* name, char* mode) {
  fFrame = new IMAGE2D;
  OpenFrame(name,mode);
  fVariance=0;
}


/* ----- ImageSimple -------------------------------------------------- */
ImageSimple::ImageSimple(const ImageSimple &image,char* newname,short newtype,int copydata) {
  // Copy constructor
  fFrame = new IMAGE2D;
  fVariance =0;
  int npix[2];
  double start[2];
  double step[2];

  npix[0] = image.Nx();
  npix[1] = image.Ny();
  start[0] = image.fFrame->startx;
  start[1] = image.fFrame->starty;
  step[0] = image.fFrame->stepx;
  step[1] = image.fFrame->stepy;
  newtype = (newtype==0) ? image.fFrame->data_type : newtype;

  create_frame(fFrame,newname,npix,start,step,newtype,image.fFrame->ident,image.fFrame->cunit);
  fLoaded=true;

  CP_non_std_desc(image.fFrame,fFrame);

  if (copydata)
    for (int col=0 ; col<image.Nx(); col ++)
      for (int line=0 ; line<image.Ny(); line ++)
        WrFrame(col,line,image.RdFrame(col,line));
}  


/* ----- ~ImageSimple ------------------------------------------------- */
ImageSimple::~ImageSimple(){
  if (fVariance)
    delete fVariance;
  if (fLoaded)
    CloseFrame();
  delete fFrame;
}

/* ===== Wrappers to IMAGE2D methods ================================== */

/* ----- WrFrame ------------------------------------------------------ */
void ImageSimple::WrFrame(int line, int col, double value){
  //Wrapper to WR_frame
  WR_frame(fFrame,line,col,value);
}

/* ----- WrFrame ------------------------------------------------------ */
void ImageSimple::WrFrame(int line, int col, long value){
  //Wrapper to WR_frame
  WR_frame(fFrame,line,col,value);
}

/* ----- RdFrame ------------------------------------------------------ */
double ImageSimple::RdFrame(int line, int col) const {
  //Wrapper to RD_frame
  return (double) RD_frame(fFrame,line,col);
}

/* ----- RdDesc ------------------------------------------------------- */
int ImageSimple::RdDesc(char* Descr, short Type, int NbElements, void* Values) const {
  // Wrapper to RdDesc
  return RD_desc(fFrame, Descr,Type,NbElements,  Values); 
}

/* ----- WrDesc ------------------------------------------------------- */
int ImageSimple::WrDesc(char* Descr, short Type, int NbElements, void* Values) {
  // Wrapper to RdDesc
  return WR_desc(fFrame, Descr,Type,NbElements,  Values);
}

/* ----- RdIfDesc ----------------------------------------------------- */
int ImageSimple::RdIfDesc(char* Descr, short Type, int NbElements, void* Values) const {
  // Variant of RdDesc which does nothing if the descriptor does not exists
  // the content of Values in that case undefined
  disable_user_warnings();
  int retCode = RD_desc(fFrame, Descr,Type,NbElements,  Values);
  restore_user_warnings();
  return retCode;
}

/* ----- DeleteDesc ----------------------------------------------------- */
int ImageSimple::DeleteDesc(char* Descr) {
  // Deletes teh descriptor
  return delete_desc(fFrame,Descr);
}

/* ----- OpenFrame ---------------------------------------------------- */
int ImageSimple::OpenFrame(char *name, char *mode){
  // Wrapper to open_frame
  fLoaded=true;
  return open_frame(fFrame,name,mode);
}

/* ----- CloseFrame --------------------------------------------------- */
int ImageSimple::CloseFrame(){
  // Wrapper to close_frame
  fLoaded=false;
  return close_frame(fFrame);
}


/* ----- CreateFrame -------------------------------------------------- */
int ImageSimple::CreateFrame(char *Name,int Nx, int Ny, short Type ){
  // this CreateFrame is less powerfull than teh original, but
  // this is because I don't understand all parameters !
  int npix[2];
  double start[2];
  double step[2];

  start[0]=npix[0]=Nx;
  start[1]=npix[1]=Ny;
  step[0]=step[1]=1;
  fLoaded=true;

  return create_frame(fFrame,Name,npix,start,step,Type,"","");
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
void ImageSimple::MinMax(Section* Sec, double * min, double * max) const {
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

void ImageSimple::ImportSection(ImageSimple * From, Section* Sec, int X1Start, int Y1Start,int XDir,int YDir,double ZScale) {
  // pilots the ImportSectionFrame for both the image and the variance 
  // if it exists
  ImportSectionFrame(From,Sec,X1Start,Y1Start,XDir,YDir,ZScale);
  if (Variance()&&From->Variance()) {
    Variance()->ImportSection(From->Variance(),Sec,X1Start,Y1Start,XDir,YDir,ZScale*ZScale);
    }
  
    

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
double ImageSimple::MeanValue(Section* Sec,int step) const{
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

/* ----- OverscanRms --------------------------------------------------- */
double ImageSimple::OverscanRms(Section* S,double sigcut) {
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
    sigma = OverscanRms(S,ut_fraction_sigcut(1.5/(S->XLength())));

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


/* ##### IMAGE SNIFS ################################################## */

/* ===== constructor/Destructor ======================================= */
ImageSnifs::ImageSnifs () {
  SetParanoMode(true);
}

/* ----- ImageSnifs copy ------------------------------ */
ImageSnifs::ImageSnifs(const ImageSnifs &image,char* newname,short newtype,int copydata) 
  : ImageSimple(image, ut_create_check_name(newname), newtype, copydata) {
  SetParanoMode(true);

  // if variance, try to copy if a meaningful naem can be done
  if (image.Variance()) {
    char varname[lg_name+1];
    ut_varname_from_imname(newname,varname);
    if (varname[0]) {
      ImageSimple* var = new ImageSimple(*image.Variance(),varname,newtype,copydata);
      SetVarianceFrame(var);
    }
  }
}

/* ----- ImageSnifs open ------------------------------ */
ImageSnifs::ImageSnifs(char* name, char* mode)
  : ImageSimple( ut_open_check_name(name), mode) {

  SetParanoMode(true);

  // load varianece if available
  char varname[lg_name+1] ;
  ut_varname_from_imname(Name(),varname);
  if (varname[0] && exist(varname)) {
    ImageSimple* var = new ImageSimple(varname,mode);
    SetVarianceFrame(var);
  }
}


/* ----- ~ImageSnifs ------------------------------ */
ImageSnifs::~ImageSnifs(){  
}

/* ===== Methods ======================================= */

/* ----- BuildSubImage  ------------------------------ */
ImageSnifs* ImageSnifs::BuildSubImage(Section* Sec,char *Name){
  // Build a raster image

  char varname[lg_name+1], secString[lg_name+1] ;
  ut_varname_from_imname(Name,varname);
  
  // raw contructor (create_frame has to be thought about)
  ImageSnifs* sub = new ImageSnifs();
  if (Variance() && varname[0]!='\0') {
    ImageSnifs* variance = new ImageSnifs();
    sub->SetVarianceFrame(variance);
  }

  // create the frames
  sub->CreateFrame(Name, Sec->XLength() ,Sec->YLength());
  sub->ImportHeader(this);
  if (sub->Variance()) {
    sub->Variance()->CreateFrame(varname, Sec->XLength() ,Sec->YLength());
    sub->Variance()->ImportHeader(Variance());
  }
  sub->ImportSection(this,Sec,1,1);

  // update the keys : RASTER
  Sec->GetString(secString);
  sub->WrDesc("RASTER",CHAR,lg_name+1,secString);
  if (sub->Variance())
    sub->Variance()->WrDesc("RASTER",CHAR,lg_name+1,secString);
  // DATASEC
  sprintf(secString,"[%d:%d,%d:%d]",1,Sec->XLength() ,1,Sec->YLength());
  sub->WrDesc("DATASEC",CHAR,lg_name+1,secString);
  if (sub->Variance())
    sub->Variance()->WrDesc("DATASEC",CHAR,lg_name+1,secString);

  return sub;
}



/* ----- SubstractOverscan ------------------------------ */
void ImageSnifs::SubstractOverscan() {

  int overscanDone;

  if (ParanoMode()) {
    // check fcalsses !

    // check overscan was not already substracted
    if ( RdIfDesc("OVSCDONE",INT, 1, &overscanDone) >= 0 
         && overscanDone ){
      print_error("Overscan already substracted for %s\n",Name());
      print_error("Nothing Done%s\n",Name());
      return;
    }
  }
  
  // Get the overscan bounds from file
  char bias[lg_name+1];
  RdDesc("BIASSEC",CHAR,lg_name+1,bias);
  Section * Sec = new Section(bias);
  ImageSimple::SubstractOverscan(Sec);  
  delete Sec;

  overscanDone=1;
  WrDesc("OVSCDONE",INT,1,&overscanDone);

  // set the biasframe if allowed
  char obstype[lg_name+1];
  RdDesc("OBSTYPE",CHAR,lg_name+1,obstype);
  if (! strcmp(obstype,"BIAS")) {
    int biasFrame=1;
    WrDesc("BIASFRAM",INT,1,&biasFrame);
  }

}

/* ----- odd-Even ----------------------------------------------------- */
void ImageSnifs::OddEvenCorrect() {
  // odd-even substraction
  //

  double param[2];
  
  if (ParanoMode()) {
    // check fcalsses !
    
    // check type of CCD (is operation allowed ?)

    // check overscan was not already substracted
    // returns 0 if not done !
    if ( RdIfDesc("OEPARAM",DOUBLE, 2, param) > 0 ){
      print_error("Odd-Even already substracted for %s\n",Name());
      print_error("Nothing Done%s\n",Name());
      return;
    }
  }

  // Get the overscan bounds from file
  char bias[lg_name+1];
  RdDesc("BIASSEC",CHAR,lg_name+1,bias);
  Section* Sec=new Section(bias);

  // the magic 2.0 parameter is because of CTE on CCD EEV#1 (blue)
  ImageSimple::OddEvenCorrect(Sec,param,2.0);
  delete Sec;

  WrDesc("OEPARAM",DOUBLE, 2, param);
}

/* -----  SubstractBias ------------------------------------------------- */
void ImageSnifs::SubstractBias(ImageSnifs* Bias) {
  // substracts a bias image
  int biasDone,overscanDone,biasFrame;

  if (ParanoMode()) {
    // check fclasses
    Bias->RdDesc("BIASFRAM",INT,1,&biasFrame);
    if (!biasFrame) {
      print_error("ImageSnifs:SubstractBias %s is not a bias frame",Bias->Name());
      return;
    }

    // check bias was not already substracted
    if ( RdIfDesc("BIASDONE",INT, 1, &biasDone) >= 0 
         && biasDone ){
      print_error("ImageSnifs::SubstractBias already substracted for %s",Name());
      return;
    }

    // check overscan was already substracted
    RdDesc("OVSCDONE",INT, 1, &overscanDone);
    if (!overscanDone ){
      print_error("Overscan was not substracted yet for %s\n",Name());
      return;
    }

    // check images are of the same size
    if (Nx()!=Bias->Nx() || Ny()!=Bias->Ny()) {
      print_error("ImageSnifs::SubstractBias : images are not of same size");
      return;
    }    
  }
  
  // do the real job
  Add(Bias,-1);

  biasDone=1;
  WrDesc("BIASDONE",INT,1,&biasDone);

  // set the darkframe if allowed
  char obstype[lg_name+1];
  RdDesc("OBSTYPE",CHAR,lg_name+1,obstype);
  if (! strcmp(obstype,"DARK")) {
    int darkFrame=1;
    WrDesc("DARKFRAM",INT,1,&darkFrame);
  }

}

/* -----  SubstractDark ------------------------------------------------- */
void ImageSnifs::SubstractDark(ImageSnifs* Dark) {
  // substracts a dark image
  int darkDone,overscanDone,darkFrame;
  double exptime,darktime;

  if (ParanoMode()) {
    // check fclasses
    Dark->RdDesc("DARKFRAM",INT,1,&darkFrame);
    if (!darkFrame) {
      print_error("ImageSnifs::SubstractDark %s is not a dark frame",Dark->Name());
      return;
    }

    // check dark was not already substracted
    if ( RdIfDesc("DARKDONE",INT, 1, &darkDone) >= 0 
         && darkDone ){
      print_error("ImageSnifs::SubstractDark already substracted for %s",Name());
      return;
    }

    // check overscan was already substracted
    RdDesc("OVSCDONE",INT, 1, &overscanDone);
    if (!overscanDone ){
      print_error("ImageSnifs::SubstractDark Overscan was not substracted yet for %s\n",Name());
      return;
    }

    // check images are of the same size
    if (Nx()!=Dark->Nx() || Ny()!=Dark->Ny()) {
      print_error("ImageSnifs::SubstractDark : images are not of same size");
      return;
    }
  }

  // do the real job
  RdDesc("EXPTIME",DOUBLE,1,&exptime);
  Dark->RdDesc("EXPTIME",DOUBLE,1,&darktime);
  Add(Dark,-exptime/darktime);

  darkDone=1;
  WrDesc("DARKDONE",INT,1,&darkDone);


}


/* -----  ApplyFlat     ------------------------------------------------- */
void ImageSnifs::ApplyFlat(ImageSnifs* Flat) {
  // apply a flat frame
  int flatDone,biasDone,flatFrame;

  if (ParanoMode()) {
    // check fclasses
    RdDesc("BIASDONE",INT, 1, &biasDone);
    if ( !biasDone ){
      print_error("ImageSnifs::ApplyFlat : Bias was not substracted yet for %s\n",Name());
      return;
    }

    Flat->RdDesc("FLATFRAM",INT, 1, &flatFrame);
    if ( !flatFrame ){
      print_error("ImageSnifs::ApplyFlat : Bias was not substracted yet for %s\n",Name());
      return;
    }

    // check bias (so also overscan) was already substracted
    RdDesc("BIASDONE",INT, 1, &biasDone);
    if ( !biasDone ){
      print_error("ImageSnifs::ApplyFlat : Bias was not substracted yet for %s\n",Name());
      return;
    }

    // check dark was not already substracted
    if ( RdIfDesc("FLATDONE",INT, 1, &flatDone) >= 0 
         && flatDone ){
      print_warning("ImageSnifs::ApplyFlat already done for %s",Name());
    }

    // check images are of the same size
    if (Nx()!=Flat->Nx() || Ny()!=Flat->Ny()) {
      print_error("ImageSnifs::ApplyFlat : images are not of same size");
      return;
    }
  }

  // do the real job
  Divide(Flat);

  flatDone=1;
  WrDesc("FLATDONE",INT,1,&flatDone);
  int eleCalib=0;
  WrDesc("ELECALIB",INT,1,&eleCalib);
}

/* -----  BuildFlat     ------------------------------------------------- */

void ImageSnifs::BuildFlat() {
  // transforms the image into a flat frame
  int biasDone,flatFrame;

  if (ParanoMode()) {
    // check fclasses
    char obstype[lg_name+1];
    RdDesc("OBSTYPE",CHAR,lg_name+1,obstype);
    if ( strcmp(obstype,"OBJECT")) {
       print_error("ImageSnifs::BuildFlat : %s is not a flat frame",Name());
      return;
    }
    // check it is not already a flat frame
    if ( RdIfDesc("FLATFRAM",INT, 1, &flatFrame) >= 0 
         && flatFrame ){
      print_warning("ImageSnifs::buildFlat %s is already a flat frame %s",Name());
    }

    // check bias was substracted
    RdDesc("BIASDONE",INT, 1, &biasDone);
    if ( !biasDone ){
      print_error("ImageSnifs::BuildFlat : Bias was not substracted yet for %s",Name());
      return;
    }
  }

  // Get the data bounds from file
  char data[lg_name+1];
  RdDesc("DATASEC",CHAR,lg_name+1,data);
  Section* Sec=new Section(data);
  // do the real job
  double norm = MeanValue(Sec,100);
  Scale(1/norm);
  delete Sec;
  
  flatFrame=1;
  WrDesc("FLATFRAM",INT,1,&flatFrame);
  int eleCalib=0;
  WrDesc("ELECALIB",INT,1,&eleCalib);
}

/* -----  AddPoissonNoise ----------------------------------------------- */

void ImageSnifs::AddPoissonNoise() {
  // Adds the poisson noise to the image.
  //
  // Side effect : image is bounded to reject negative values.

  int poisNoise=0;

  // ParanoChecks
  if (ParanoMode()) {
    // check if in electrons
    // and if noise was not already added
    int eleCalib=0;
    if ( RdIfDesc("ELECALIB",INT, 1, &eleCalib) <= 0 
         || !eleCalib) {
      print_error("ImageSnifs::AddPoissonNoise %s is not in electrons",Name());
      return;
    }
    if ( Variance()->RdIfDesc("POISNOIS",INT, 1, &poisNoise) > 0 
         && poisNoise) {
      print_error("ImageSnifs::AddPoissonNoise already performed for %s",Name());
      return;
    }
  }
  if (!Variance()) {
    print_error("ImageSnifs::AddPoissonNoise but %s has no variance",Name());
    return;  
  }
  // End of checks

  // Go
  ImageSimple::AddPoissonNoise();

  // Set variables
  poisNoise=1;
  Variance()->WrDesc("POISNOIS",INT, 1, &poisNoise);
}

/* -----  AddPoissonNoise ----------------------------------------------- */
void ImageSnifs::HandleSaturation() {
  // puts variance to infinity for saturated pixels
  int saturate;
  RdDesc("SATURATE",INT,1,&saturate);
  ImageSimple::HandleSaturation(saturate-0.5);
}



/* ===== Hacks ======================================= */

/* ----- HackFitsKeywords ---------------------------------------------- */
void ImageSnifs::HackFitsKeywords() {

  // The following is not valid for any raster yet !
  char key[lg_name+1];
  RdDesc("RASTER",CHAR,lg_name+1,key);
  if (strcmp(key,"FULL")){
    print_error("ImageSnifs::HackFitsKeywords : only works for FULL raster");
    return;
  }

  // The BiasSec is not correct
  int x1,x2,y1,y2;
  RdDesc("BIASSEC",CHAR,lg_name+1,key);
  sscanf(key,"[%d:%d,%d:%d]",&x1,&x2,&y1,&y2);
  sprintf(key,"[%d:%d,%d:%d]",x1+2,x2,1,Ny());
  WrDesc("BIASSEC",CHAR,lg_name+1,key);
  
  // the DATASEC needs help too
  RdDesc("DATASEC",CHAR,lg_name+1,key);
  sscanf(key,"[%d:%d,%d:%d]",&x1,&x2,&y1,&y2);
  sprintf(key,"[%d:%d,%d:%d]",4,x2+2,1,y2+5);
  WrDesc("DATASEC",CHAR,lg_name+1,key);

}


/* ===== Processing Managment ======================================= */

/* ----- CreateVarianceFrame ------------------------------ */
void ImageSnifs::CreateVarianceFrame(char* name) {

  if (ParanoMode()) {
    if (Variance()) {
      print_error(" ImageSnifs::CreateVarianceFrame : already a variance in %s",Name());
      return;
    }
  
    if (!name[0]) {
      ut_varname_from_imname(Name(),name);
      if (!name[0]) {
         print_error(" ImageSnifs::CreateVarianceFrame : no name provided",Name());
      return;
      }
    }

  } // ParanoMode

  // Builds the variance frame
  ImageSimple * var = new ImageSimple(*this , name,FLOAT,0);
  SetVarianceFrame(var);
  HandleSaturation();
  
}


/* ----- AddOverscanVariance ------------------------------ */
void ImageSnifs::AddOverscanVariance() {

  int ovscNoise;
  if (ParanoMode()) {
    if (Variance()->RdIfDesc("OVSCNOIS",INT, 1, &ovscNoise) > 0 
     && ovscNoise ) {
       print_error(" ImageSnifs::AddOverscanVariance already done in %s",Name());
      return;
      }
  }

  // Fills with analysis of RMS strip
  char bias[lg_name+1];
  RdDesc("BIASSEC",CHAR,lg_name+1,bias);
  Section * Sec = new Section(bias);
  // we allow 1 outlier...
  double rms = OverscanRms(Sec,ut_fraction_sigcut(1.0/(Sec->XLength()-0.5)));
  WrDesc("RDNOISE",DOUBLE,1,&rms);
  Variance()->Add(rms*rms);
  delete Sec;

}

/*====== Specific photometry tools ================== */

/* ----- BuildFilterSection ------------------------------ */
void ImageSnifs::FilterSectionFill(Section * Sec,int Filter) {

  // Sec will be overwritten in the process

  switch (Filter) {
  case 0 : 
    Sec->SetString("[1:2048,1:686]");
    break;
  case 1 :
    Sec->SetString("[1:772,1192:2393]");
    break;
  case 2:
    Sec->SetString("[1277:2048,1192:2393]");
    break;
  case 3:
    Sec->SetString("[1:772,2899:4096]");
    break;
  case 4:
    Sec->SetString("[1:772,2899:4096]");
    break;
  default:
    print_error("ImageSnifs::BuildFilterSection : paramter out of 0-4 range");
  }
}

/* ----- SplitMultifilter ------------------------------ */
void ImageSnifs::SplitMultiFilter(char* NameRecipee,ImageSnifs** Subs){
  // the NameRecipee shall contain one %s which will
  // be transformated with 'a' ... 'e'
  // The Subs should contain 5 allocated ImageSnifs*

  // as usual, paranochecks
  if (ParanoMode()) {
  }
  // and checks
  if (!strstr(NameRecipee,"%c")){
    print_error("ImageSnifs::SplitMultifilter Name shall contain %c");
    return;
  }

  // Do it :
  char imageName[lg_name+1];
  Section sec;

  for (int i=0;i<5;i++) {
    // names for the image
    sprintf(imageName,NameRecipee,'a'+i);
    FilterSectionFill(&sec,i);
    Subs[i] = BuildSubImage(&sec,imageName);
  }
}


