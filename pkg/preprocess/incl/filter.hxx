/* === Doxygen Comment ======================================= */
/*! 
 * \file          filter.hxx
 * \copyright     (c) 2004 SNIFS Collaboration
 * \date          Wed Aug  6 18:32:01 2003
 * \author        Emmanuel Gangler <e.gangler@ipnl.in2p3.fr>
 * \version       0.0
 * \brief         
 *                
 * $Id$
 **/
/* =========================================================== */

#ifndef FILTER_H
#define FILTER_H
/*
 Filter: a simple filtering interface
*/   

/* ----- local includes ----- */
class Section;
class ImageSimple;
class ImageAnalyser;

/* ===== Filter ============================== */

class ImageFilter {
  public :
  enum Bound_t { kNoData, kShrinks, kSymetric};
  
  virtual ~ImageFilter() {};

  void SetWindowSize(int X,int Y) {fXsize=X;fYsize=Y;} // the HALF
  // size, i.e. 3x3 = (1,1)

  void SetNoDataAnswer(double Val, double Var) {fNoData=Val;fNoDataVar=Var;}
  //  void  GetWindowBounds(int Xcenter, int Ycenter, int* Xlow, int* Xup, int* Ylow, int* Yup);
  void SetBoundaryWay(Bound_t B) {fBound = B;}
  virtual void SetInputImage(ImageSimple* I) {fInput = I;}
  void SetOutputImage(ImageSimple* I) {fOutput = I;}
  virtual void Filter(int X, int Y, Section* S)=0;
  void Filter();

protected:
  void WriteNoData(int I, int J);
  
  int fXsize, fYsize;
  Bound_t fBound;
  double fNoData;
  double fNoDataVar;
  // running paramters
  ImageSimple* fInput;
  ImageSimple* fOutput;
};

/* ===== FilterHF ==================== */

// principle : computes the Pixel divided by the median of the window
//    only for pixel with value > threshold
//    only if result has a signicicance (in sigma) above threshold
//    the median variance is set to mean variance * sqrt(3/N)

class ImageFilterHF : public ImageFilter {
  public :

  ImageFilterHF(int Xsize, int Ysize, Bound_t B=kNoData);
  ~ImageFilterHF();

  void SetThreshold(double Value) {fThreshold = Value;}
  void SetSignificance(double Sigma) {fSignificance=Sigma;}
  virtual void SetInputImage(ImageSimple* I);
  virtual void Filter(int X, int Y, Section* S);
  

protected:
  double fThreshold;
  double fSignificance;
  ImageAnalyser* fAnal;

};

/* ===== FilterMax ==================== */

// returns the maximum in the window

class ImageFilterMax : public ImageFilter {
  public :

  ImageFilterMax(int Xsize, int Ysize, Bound_t B=kNoData);
  ~ImageFilterMax();

  //  virtual void SetInputImage(ImageSimple* I);
  virtual void Filter(int X, int Y, Section* S);
  

protected:
  ImageAnalyser* fAnal;

};

/* ===== FilterMedian ==================== */

// returns the maximum in the window

class ImageFilterMedian : public ImageFilter {
  public :

  ImageFilterMedian(int Xsize, int Ysize, Bound_t B=kNoData);
  ~ImageFilterMedian();

  virtual void SetInputImage(ImageSimple* I);
  virtual void Filter(int X, int Y, Section* S);
  

protected:
  ImageAnalyser* fAnal;

};

#endif
