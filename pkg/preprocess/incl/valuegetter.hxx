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


#endif
