/* === Doxygen Comment ======================================= */
/*! 
 * \file          lrerror.cxx
 * \copyright     (c) 2004 SNIFS Collaboration
 * \date          Wed Aug  6 17:44:35 2003
 * \author        Emmanuel Gangler <e.gangler@ipnl.in2p3.fr>
 * \version       0.0
 * \brief         
 *                
 * $Id$
 **/
/* =========================================================== */

/* The purpose of this main is to compute the error from extrapolating left vs right pixel with respect to the given pixel*/

/* ----- local include ----- */
#include "imagesnifs.hxx"
#include "filter.hxx"
#include "catorfile.hxx"
#include "utils.h"
#include <gsl/gsl_fit.h>

/* ----- main ---------------------------------------- */

int main(int argc, char **argv) {

  char **argval, **arglabel;

  set_arglist("-in none -out none");
  init_session(argv,argc,&arglabel,&argval);

  char inName[lg_name+1],outName[lg_name+1];

  CatOrFile inCat(argval[0]);
  CatOrFile outCat(argval[1]);

  while (inCat.NextFile(inName) && outCat.NextFile(outName)) {
    print_msg("Filtering %s",inName);
    ImageSnifs *in = new ImageSnifs(inName,"I");
    ImageSimple *var = in->Variance();
    if (!in->Variance())
      print_error("input frame has no variance");
    ImageSnifs *out = new ImageSnifs(*in,outName);
    for (int i=1;i<in->Nx()-1;i++){
      for (int j=1;j<in->Ny()-1;j++){
        //        double error=(in->RdFrame(i+1,j)+in->RdFrame(i-1,j))/2 - in->RdFrame(i,j);
        //        double norm=sqrt((var->RdFrame(i+1,j)+var->RdFrame(i-1,j))/4 + var->RdFrame(i,j));

#ifdef A_COMPUTATION
        double val;
        if (in->RdFrame(i,j)!=0)
          val=(in->RdFrame(i+1,j)+in->RdFrame(i-1,j))/2 / in->RdFrame(i,j);
        else
          val=0;
        
        double variance=( (var->RdFrame(i+1,j)+var->RdFrame(i-1,j))/4 
                     + val*val*var->RdFrame(i,j) ) 
          / (var->RdFrame(i,j)* var->RdFrame(i,j));
#endif
#ifdef B_COMPUTATION
        double val;
        double variance;
        val = in->RdFrame(i,j)/in->RdFrame(i,2939);
#endif        
#ifdef C_COMPUTATION
        double val;
        double variance;
        val = in->RdFrame(i+1,j)/in->RdFrame(i-1,j);
#endif
        double val;
        int width=5;
        double profile[2][width];
        double x[width];
        double corr;
        for (int k=0;k<width;k++) { // 5 has to be tuned !
          // forget -1 as it may be contaminated
          profile[0][k]=in->RdFrame(i,j+1+k);
          profile[1][k]=in->RdFrame(i,j-1-k);
          x[width]=k+1;
        }
        //        double c00,c01,c1,cov00,cov01,cov11,sumsq;
        //gsl_fit_linear (x,1,profile[0],1,width, &c00,&c1,&cov00,&cov01, &cov11, &sumsq);
        //gsl_fit_linear (x,1,profile[1],1,width, &c01,&c1,&cov00,&cov01, &cov11, &sumsq);
        
        // compute the correction
        corr=ut_median(profile[0],width)-ut_median(profile[1],width);
        //val=corr;
        val=corr/in->RdFrame(i,j);

        out->WrFrame(i,j,val);
        //        out->Variance()->WrFrame(i,j,variance);
      }
    }
    delete in;
    delete out;
  }

  exit_session(0);
  
}
