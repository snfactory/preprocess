/* === Doxygen Comment ======================================= */
/*! 
 * \file          imagestack.hxx
 * \copyright     (c) 2003 CRAL-Observatoire de Lyon
 * \date          Wed Aug  6 18:32:01 2003
 * \author        Emmanuel Gangler <e.gangler@ipnl.in2p3.fr>
 * \version       0.0
 * \brief         
 *                
 * $Id$
 **/
/* =========================================================== */

#ifndef IMAGESTACKSNIFS_H
#define IMAGESTACKSNIFS_H
/*   
     ImageStackSnifs is the SNIFS processor for ImageStack.
*/


/* ===== ImageStackSnifs ======================================== */

class ImageStackSnifs  {
  public :

    ImageStackSnifs(CatOrFile * Cat, char * Mode="I");
  
  protected :
    ImageStack fImages;

};

#endif
