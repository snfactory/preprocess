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

  There is a special case of bichips having 4 amplifiers : this corresponds to 
  the phtometric + guider images.
*/

/* ----- local includes -----*/
class ImageSnifs;
class Section;
class AlgoCams;
#include "image.hxx"

/* ===== Bi Chip Snifs ============================== */

class BiChipSnifs {
  public :

    BiChipSnifs(int NChips=2);
    BiChipSnifs(char* ChipNameRecipee,char* mode="Input", IoMethod_t Method=kIoPlain,int MParam=0);
    BiChipSnifs(const BiChipSnifs &image,char* NewNameRecipee,short newtype = 0,int copydata=0,IoMethod_t Method=kIoPlain,int MParam=0);
    ~BiChipSnifs();

    // Getters
    int NChips() const {return fNChips;}

    // Utilities
    void CheckNameRecipee(char* ChipNameRecipee);  
    ImageSnifs* Chip(int chip)  {return fChip[chip];}
    void SetChip(int chip, ImageSnifs* Image)  {fChip[chip]=Image;}
    void SetNLines(int NLines);
    void SetParanoMode(bool Mode);
    void SetAlgo(char* AlgoName);

    // Methods

    void UpdateFClass();
    //void SubstractOverscan();
    // void OddEvenCorrect();
    void SubstractBias(BiChipSnifs* Bias);
    void SubstractDark(BiChipSnifs* Dark);

    ImageSnifs* Assemble(char* ImageName,IoMethod_t Io=kIoPlain,int nLines=kIoAll);
    double GuessGainRatio(Section* S);

    void CreateVarianceFrame(char* VarianceNameRecipee="");
    void HandleSaturation();

    // Hacks

    void HackFitsKeywords( char * PrimaryName=0 );
    void HackGainRatio();

  protected :
    ImageSnifs** fChip;
    int fNChips;
  
};


#endif
