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

/* ----- local includes -----*/
class ImageSnifs;
class Section;
class AlgoCams;
#include "image.hxx"

/* ===== Bi Chip Snifs ============================== */

class BiChipSnifs {
  public :

    BiChipSnifs();
    BiChipSnifs(char* ChipNameRecipee,char* mode="Input", IoMethod_t Method=kIoPlain,int MParam=0);
    BiChipSnifs(const BiChipSnifs &image,char* NewNameRecipee,short newtype = 0,int copydata=0,IoMethod_t Method=kIoPlain,int MParam=0);
    ~BiChipSnifs();

    // Utilities
    void CheckNameRecipee(char* ChipNameRecipee);  
    ImageSnifs* Chip(int chip)  {return fChip[chip];}
    void SetChip(int chip, ImageSnifs* Image)  {fChip[chip]=Image;}
    void SetNLines(int NLines);
    void SetParanoMode(bool Mode);
    void SetAlgo(char* AlgoName);

    // Methods

    void SubstractOverscan();
    void OddEvenCorrect();
    void SubstractBias(BiChipSnifs* Bias);

    ImageSnifs* Assemble(char* ImageName,IoMethod_t Io=kIoPlain,int nLines=kIoAll);
    double GuessGainRatio(Section* S);

    void CreateVarianceFrame(char* VarianceNameRecipee="");
    void HandleSaturation();
    void AddOverscanVariance();
    // void PreprocessBias();
    //ImageSnifs* PreprocessAssemble(char* OutName,BiChipSnifs *bias);
    //ImageSnifs* PreprocessDark(char* OutName,BiChipSnifs *bias);
    //ImageSnifs* PreprocessFlat(char* OutName,BiChipSnifs *bias,ImageSnifs *dark);
    //ImageSnifs* Preprocess(char* OutName,BiChipSnifs *bias,ImageSnifs *dark,ImageSnifs *flat);

    // Hacks

    void HackFitsKeywords();
    void HackGainRatio();

  protected :
    ImageSnifs* fChip[2];
};


#endif
