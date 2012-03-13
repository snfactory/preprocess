/* === Doxygen Comment ======================================= */
/*! 
 * \file          filter_imaged.cxx
 * \copyright     (c) 2004 SNIFS Collaboration
 * \date          Wed Aug  6 17:44:35 2003
 * \author        Emmanuel Gangler <e.gangler@ipnl.in2p3.fr>
 * \version       0.0
 * \brief         
 *                
 * $Id$
 **/
/* =========================================================== */

/* ----- std includes ----- */
#include <stdlib.h>
#include <vector>
using namespace std;

/* ----- local include ----- */
//#include <gsl/gsl_vector.h>
//#include <gsl/gsl_matrix.h>
//#include <gsl/gsl_linalg.h>


/* ----- local include ----- */
#include "imagesnifs.hxx"
#include "imagestacksnifs.hxx"
#include "filter.hxx"
#include "catorfile.hxx"
#include "kombinator.hxx"
//#include "utils.h"


/* ----- ScaleBy ---------------------------------------- */
void ScaleBy(vector<ImageSnifs*>* ins, vector<ImageSimple*>*refs, double cut) 
{
  vector<ImageSnifs*>::iterator iter;
  vector<ImageSimple*>::iterator iterRef;

  vector<gsl_vector*> coefs ;
  for (iter = ins->begin();
       iter != ins->end();++iter){
    // Guess the coefficients : get first the biggest vector
    gsl_vector * A=(*iter)->FitBy(*refs,cut);
    coefs.push_back(A);
    printf("image %s has a value of  ",(*iter)->Name());
    unsigned int niter;
    for (iterRef = refs->begin(),niter=0;niter<A->size-1;++iterRef,niter++) {
      (*iter)->Add((*iterRef),-gsl_vector_get(A,niter));
      printf("%f  ",gsl_vector_get(A,niter));
    }
    
    printf("%f\n",gsl_vector_get(A,niter));
    //(*iter)->Scale(1/gsl_vector_get(A,niter));
    
    //gsl_vector_free(A);
  }
  
  vector<gsl_vector*>::iterator itercoef ;
  double norm=0;
  for (itercoef = coefs.begin();
       itercoef != coefs.end();++itercoef){
    double elem=gsl_vector_get((*itercoef),(*itercoef)->size-1);
    norm+=elem*elem;
  }
  norm/=coefs.size();
  for (itercoef = coefs.begin(), iter = ins->begin();
       itercoef != coefs.end();++itercoef,++iter){
    (*iter)->Scale(norm/gsl_vector_get((*itercoef),(*itercoef)->size-1));
    gsl_vector_free(*itercoef);
  }
}


/* ----- main ---------------------------------------- */

int main(int argc, char **argv) {

  char **argval, **arglabel;
  double cut;
  char tmpname[lg_name+1],outName[lg_name+1];
  strcpy(tmpname,"tmp.fits");

  set_arglist("-intmp none -out none -cut 3.0 -darkmaps null -ndepth 1 -nlines 4096");
  // -intmp : list of real exposures (with cosmics and so on, but biasz removed)
  //    note that this input will be MODIFIED by the program
  // -out : the newly produced dark component
  // -cut : number of sigmas to consider an outlier
  // -darkmaps : list of dark maops already produced
  //    The last one may be a first guess for the newlycreated map
  // -ndepth : how much maps are produced. 
  //     -darkmaps should contain ndepth or ndepth-1 objects
  // - nlines : internal number to know how many lines of image 
  //     should be stored at once


  init_session(argv,argc,&arglabel,&argval);

  //char inName[lg_name+1],outName[lg_name+1],refName[lg_name+1];

  CatOrFile inCat(argval[0]);
  CatOrFile darkCat(argval[3]);

  vector<ImageSnifs *> *ins;
  vector<ImageSimple *> darks;
  vector<ImageSnifs*>::iterator iter;
  //while (refCat.NextFile(refName)) {
  //  Refs.push_back(new ImageSnifs(refName,"I"));
  //}
 
  strcpy(outName,argval[1]);
  cut = atof(argval[2]);
  unsigned int ndepth = atoi(argval[4]);
  int nlines = atoi(argval[5]);
  

  ImageStackSnifs* darkMaps;
  if is_set(argval[4]) {
    darkMaps=new ImageStackSnifs(&darkCat,"I",nlines);
    }
  else {
    darkMaps= new ImageStackSnifs();
  }
  for (iter=darkMaps->GetImages().begin();iter!=darkMaps->GetImages().end();++iter) {
    darks.push_back(*iter);
  }
  
  if (darks.size()!=ndepth) {
    print_error("test_dark : mismatch between -darkmaps size and ndepth");
    return 1;
  }
  
  // loading now the input stack
  ImageStackSnifs* inStack=new ImageStackSnifs(&inCat,"IO",nlines);
  ins=&inStack->GetImages();

  // Now taking the first coefficients
  ScaleBy(ins,&darks,cut);
  delete darks.back();
  darks.pop_back();
  // build now the map
  KGauss k(cut);
  ImageSnifs* outRef=inStack->Kombine(outName,&k);
  // and bootstrap once
  darks.push_back(outRef);
  ScaleBy(ins,&darks,cut);
  delete darks.back();
  darks.pop_back();
  outRef=inStack->Kombine(outName,&k);
  // and bootstrap 2ce (useless)
  darks.push_back(outRef);
  ScaleBy(ins,&darks,cut);
  delete darks.back();
  darks.pop_back();
  outRef=inStack->Kombine(outName,&k);

  // cleanup
  delete inStack;
  delete outRef;
  
   
  exit_session(0);
  
}
