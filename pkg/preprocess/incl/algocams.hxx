/* === Doxygen Comment ======================================= */
/*! 
 * \file          algocams.hxx
 * \copyright     (c) 2003 CRAL-Observatoire de Lyon
 * \date          Wed Aug  6 18:32:01 2003
 * \author        Emmanuel Gangler <e.gangler@ipnl.in2p3.fr>
 * \version       0.0
 * \brief         
 *                
 * $Id$
 **/
/* =========================================================== */

#ifndef ALGOCAMS_H
#define ALGOCAMS_H

/* ----- Local includes and definitions ----- */
class ImageSnifs;
class Section;


/*   
 Algocams contains the algorithms specific for cameras, which have to be 
 highly configurable, without disturbing the flow of the preprocessing.
*/

/* ===== ALGOCAMS ======================================== */

enum AlgoMethod_t {kNone,kDetcom,kOtcom};
/* The virtual class*/

class AlgoCams {
  public :
  
  
   
    // Constructors/Destructors
    virtual ~AlgoCams() {}

    // Id
    virtual AlgoMethod_t GetId()   {return kNone;}

    // method

    virtual void HackFitsKeywords(ImageSnifs * Image)=0;
  //  virtual void HackGainRatio(ImageSnifs * Image)=0;
  //    virtual void OddEvenCorrect(ImageSnifs* Image)=0;
  //  virtual void AddOverscanVariance(ImageSnifs* image)=0;
  //  virtual void SubstractOverscan(ImageSnifs* image)=0;
  virtual Section * SafeOverscanStrip(ImageSnifs* image)=0;
  
};

class AlgoDetcom : public AlgoCams{
  public :
   
    // Constructors/Destructors
    virtual ~AlgoDetcom() {}
  
    virtual AlgoMethod_t GetId()   {return kDetcom;}

    // method

    virtual void HackFitsKeywords(ImageSnifs * Image);
  // virtual void OddEvenCorrect(ImageSnifs* Image);
  //  virtual void AddOverscanVariance(ImageSnifs* Image);
  //  virtual void SubstractOverscan(ImageSnifs* Image);
    virtual Section* SafeOverscanStrip(ImageSnifs* Image);
  
};

class AlgoOtcom :public AlgoCams {
  public :
   
    // Constructors/Destructors
    virtual ~AlgoOtcom() {}
  
    virtual AlgoMethod_t GetId()   {return kOtcom;}
    // method

    virtual void HackFitsKeywords(ImageSnifs * Image);
  //virtual void OddEvenCorrect(ImageSnifs* Image);
  //  virtual void AddOverscanVariance(ImageSnifs* Image);
  //  virtual void SubstractOverscan(ImageSnifs* Image);
    virtual Section* SafeOverscanStrip(ImageSnifs* Image);
  
};

#endif
