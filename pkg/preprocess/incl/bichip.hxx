/* === Doxygen Comment ======================================= */
/*! 
 * \file          bichip.hxx
 * \copyright     (c) 2003 SNIFS-Supernova Factory Experiment
 * \date          Wed Aug  8 18:32:01 2003
 * \author        Emmanuel Gangler <e.gangler@ipnl.in2p3.fr>
 * \version       0.0
 * \brief         
 *                
 * $Id$
 **/
/* =========================================================== */

#ifndef BICHIP_H
#define BICHIP_H
/*
  BiChip : methods that are relevant to a CCd image split in 2 readout preamps.
  We suppose here the SNIFS environment
*/

class ImageSnifs;
class Section;

/* ===== Bi Chip Snifs ============================== */

class BiChipSnifs {
  public :

    BiChipSnifs();
    BiChipSnifs(char* ChipNameRecipee,char* mode="Input");
    BiChipSnifs(const BiChipSnifs &image,char* NewNameRecipee,short newtype = 0,int copydata=0);
    ~BiChipSnifs();

    // Utilities
    void CheckNameRecipee(char* ChipNameRecipee);  
  ImageSnifs* Chip(int chip)  {return fChip[chip];}
  

    // Methods

    void SubstractOverscan();
    void OddEvenCorrect();
    void SubstractBias(BiChipSnifs* Bias);

    ImageSnifs* Assemble(char* ImageName);
    double GuessGainRatio(Section* S);

    void CreateVarianceFrame(char* VarianceNameRecipee="");
    void AddOverscanVariance();
    void PreprocessBias();
    ImageSnifs* PreprocessDark(char* OutName,BiChipSnifs *bias);
    ImageSnifs* PreprocessFlat(char* OutName,BiChipSnifs *bias,ImageSnifs *dark);
    ImageSnifs* Preprocess(char* OutName,BiChipSnifs *bias,ImageSnifs *dark,ImageSnifs *flat);

    // Hacks

    void HackFitsKeywords();
    void HackGainRatio();

  protected :
    ImageSnifs* fChip[2];
};


#endif
