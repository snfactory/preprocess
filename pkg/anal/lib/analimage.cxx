/* === Doxygen Comment ======================================= */
/*! 
 * \file          analimage.cxx
 * \copyright     (c)  2003 SNIFS-Supernova Factory Experiment
 * \date          Wed Aug  6 18:32:01 2003
 * \author        Emmanuel Gangler <e.gangler@ipnl.in2p3.fr>
 * \version       0.0
 * \brief         
 *                
 * $Id$
 **/
/* =========================================================== */

/* ----- ROOT includes ------------------------------ */
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"

/* ----- image includes ------------------------------ */
#include "imagesnifs.hxx"

/* ----- local includes ------------------------------ */
#include "analimage.hxx"

ClassImp(ImageSignature)

/* ===== AnalImageSignature ============================== */

/* ----- constructor -----------------------------------------*/
AnalImageSignature::AnalImageSignature(TFile * Rfile)
{
  fImage=0;
  fSignature=new ImageSignature();
  fFile = Rfile;
  fTree=0;
}

/* ----- constructor -----------------------------------------*/
AnalImageSignature::~AnalImageSignature()
{
  delete fSignature;
}

/* ----- SetImage -----------------------------------------*/
void AnalImageSignature::SetImage(ImageSnifs * Image) 
{
  fImage = Image;
  if (fFile) {
    fTree = (TTree*) fFile->Get("Signature");
    if (!fTree) {
      fTree = new TTree("Signature","ImageSignatures");
    }
    TBranch * branch;
    char imageType[lg_name+1];
    ImageType(imageType);
    if ( (branch = fTree->GetBranch(imageType)) )
      branch->SetAddress(fSignature);
    else
      branch = fTree->Branch(imageType,"ImageSignature",&fSignature);
  }
}


/* ----- ImageType -----------------------------------------*/
void AnalImageSignature::ImageType(char* TypeString) 
  // Splits the image with its type
{
  strcpy(TypeString,"Normal");
}

/* ----- FillImageSignature -----------------------------------------*/
void AnalImageSignature::FillImageSignature() 
  // Splits the image with its type
{
  fSignature->fImageType=-1;
  fSignature->fExpTime=0;
  fSignature->fJulDate=1;
  fSignature->fRobustMean=2;
  fSignature->fRobustRMS=3;
  fSignature->fQuant99=4;
  fSignature->fQuant999=5;
  fSignature->fSatu=6;
  fSignature->fFocus=7;
  fSignature->fOverscanLevel=8;
  fSignature->fNoise=9;
  fSignature->fGain=10;
  fSignature->fColor=11;
  
}

/* ----- FillImageSignature -----------------------------------------*/
void AnalImageSignature::StoreImageSignature() 
{
  if (fTree)
    fTree->Fill();
}

