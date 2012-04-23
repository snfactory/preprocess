/* === Doxygen Comment ======================================= */
/*! 
 * \file          valuegetter.hxx
 * \copyright     (c) 2003 CRAL-Observatoire de Lyon
 * \date          Wed Aug  6 18:32:01 2003
 * \author        Emmanuel Gangler <e.gangler@ipnl.in2p3.fr>
 * \version       0.0
 * \brief         
 *                
 * $Id$
 **/
/* =========================================================== */

#ifndef VALUEGETTER_H
#define VALUEGETTER_H
/*   
     ValueGetter is a simple service class whose purpose is provide a value
     for a given image
*/

class ImageSimple;
class ImageAnalyser;
class Section;
class DarkModel;

/* ===== ValueGetter ======================================== */

class ValueGetter  {
  public :
  virtual ~ValueGetter() {};
  virtual double GetValue(ImageSimple* image)=0;
};

/* ===== ValueAnalyserMean ======================================== */

class ValueAnalyserMean : public ValueGetter {
  public :

    ValueAnalyserMean(Section* Sec);
    virtual ~ValueAnalyserMean();
    virtual double GetValue(ImageSimple* image);

  protected :
    ImageAnalyser* fAnalyser;
  
};

/* ===== ValueDark ======================================== */

class ValueDark : public ValueGetter {
  public :

    ValueDark(DarkModel* darkModel);
    virtual ~ValueDark();
    virtual double GetValue(ImageSimple* image);

  protected :
    ImageAnalyser* fAnalyser;
    DarkModel* fDarkModel;

};


/* ##### Get multiple values ################################################# */

/* ===== ValuesGetter ======================================== */

class ValuesGetter  {
  public :
  virtual ~ValuesGetter() {};
  virtual void GetValues(ImageSimple* image, gsl_vector* retValues)=0;
  virtual int NParams()=0;
};

/* ===== ValuesGetterDarkFitter ======================================== */

class ValuesGetterDarkFitter : public ValuesGetter {
  // Dark fitter needs 0=1 1=timeonterm 2 =timeonterm^2 3=tempterm
  public :
  ValuesGetterDarkFitter(DarkModel* Model, int* activate=0);
  virtual ~ValuesGetterDarkFitter() {};
  virtual void GetValues(ImageSimple* image, gsl_vector* retValues);
  virtual int NParams() {return fNParams;};

  protected:
  // for the moment, only the 1st section of the model is considered
  DarkModel* fDarkModel;
  int fNParams;
  int fActive[3];
};


#endif
