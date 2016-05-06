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
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>


/* ----- local include ----- */
#include "imagesnifs.hxx"
//#include "filter.hxx"
#include "catorfile.hxx"
#include "utils.h"

/* ----- main ---------------------------------------- */

int main(int argc, char **argv) {

  char **argval, **arglabel;
  double cut;

  set_arglist("-in none -ref none -out none -cut 3.0");

  init_session(argv,argc,&arglabel,&argval);

  char inName[lg_name+1],outName[lg_name+1],refName[lg_name+1];

  CatOrFile inCat(argval[0]);
  CatOrFile outCat(argval[2]);

  CatOrFile refCat(argval[1]);
  vector<ImageSnifs *> Refs;
  vector<ImageSnifs*>::iterator iter, iter2;
  while (refCat.NextFile(refName)) {
    Refs.push_back(new ImageSnifs(refName,"I"));
  }
 
  cut = atof(argval[3]);
  
  gsl_matrix * XX = gsl_matrix_alloc(Refs.size(),Refs.size());
  gsl_vector * XY = gsl_vector_alloc(Refs.size());
  gsl_vector * A = gsl_vector_alloc(Refs.size());
  // Y=XA

  while (inCat.NextFile(inName) && outCat.NextFile(outName)) {
    print_msg("Fitting %s",inName);
    ImageSnifs *in = new ImageSnifs(inName,"I");
    ImageSnifs *out = new ImageSnifs(*in,outName,0,1);

    double num=0;
    double denom=0;
    int nbad;
      int niter,niter2;
    do {
      gsl_matrix_set_zero(XX);
      gsl_vector_set_zero(XY);

      for (iter = Refs.begin(),niter=0;iter != Refs.end();++iter,niter++){
        for (int j=0;j<in->Ny();j++) {
          for (int i=0;i<in->Nx();i++) {
            double var = in->Variance()->RdFrame(i,j);
            if (var>ut_big_value/10)
              var=ut_big_value;
            else
              var=1;
            (*gsl_vector_ptr(XY,niter))+=(*iter)->RdFrame(i,j)*in->RdFrame(i,j)/var;
          }
        }
        for (iter2 = iter,niter2=niter;iter2 != Refs.end();++iter2,niter2++){
          for (int j=0;j<in->Ny();j++) {
            for (int i=0;i<in->Nx();i++) {
              double var = in->Variance()->RdFrame(i,j);
              if (var>ut_big_value/10)
                var=ut_big_value;
              else
                var=1;
              
              (*gsl_matrix_ptr(XX,niter,niter2))+=(*iter)->RdFrame(i,j)*(*iter2)->RdFrame(i,j)/var;      
            }
            
          }
        }
      }
      for (unsigned int n1=0;n1<Refs.size();n1++) {
        for (unsigned int n2=n1+1;n2<Refs.size();n2++) {
          gsl_matrix_set(XX,n2,n1,gsl_matrix_get(XX,n1,n2));
        } 
      }
      
      gsl_linalg_cholesky_decomp (XX);
      gsl_linalg_cholesky_solve (XX,XY,A);

          //num+=in->RdFrame(i,j)*ref->RdFrame(i,j)/variance;
          //denom+=ref->RdFrame(i,j)*ref->RdFrame(i,j)/variance;
      for (unsigned int n=0;n<Refs.size();n++)
        print_msg("Coefficient %f",gsl_vector_get(A,n));
      // removing far
      nbad=0;
      for (int j=0;j<in->Ny();j++) {
        for (int i=0;i<in->Nx();i++) {
          double variance = in->Variance()->RdFrame(i,j);
          double delta=in->RdFrame(i,j);
          for (iter = Refs.begin(),niter=0;iter != Refs.end();++iter,niter++)
            delta -= gsl_vector_get(A,niter)*(*iter)->RdFrame(i,j);
          if (delta > cut*sqrt(variance)) {
            in->Variance()->WrFrame(i,j,ut_big_value);
            nbad++;
          }
        }
      }
      print_msg("removed %d pixels",nbad);
    }
    while(nbad>0);

    for (iter = Refs.begin(),niter=0;iter != Refs.end();++iter,niter++)
      out->Add((*iter),-gsl_vector_get(A,niter));


    delete in;
    delete out;
    
    }

  for (iter = Refs.begin();iter != Refs.end();++iter)
    delete *iter;
   
    
    
  exit_session(0);
  
}
