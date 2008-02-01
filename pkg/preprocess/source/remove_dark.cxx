/* === Doxygen Comment ======================================= */
/*! 
 * \file          remove_dark.cxx
 * \copyright     (c) 2007 SNIFS Collaboration
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
#include "imagestacksnifs.hxx"
#include "filter.hxx"
#include "catorfile.hxx"
#include "kombinator.hxx"
#include "analyser.hxx"
#include "section.hxx"
#include "utils.h"


#ifdef REMEMBER
void removeSmart() {
    // data for B channel ;)
    // data for B channel
    double spread[42]={1.0,
                       0.69255299539170512,
                       0.16123502304147466,
                       0.021837173579109064,
                       0.0051981566820276494,
                       0.0028878648233486943,
                       0.0020706605222734257,
                       0.0020706605222734257,
                       0.0018740399385560676,
                       0.0015545314900153612,
                       0.0015545314900153612,
                       0.0012104454685099848,
                       0.0012104454685099848,
                       0.0012104454685099848,
                       0.0012104454685099848,
                       0.0010384024577572964,
                       0.0010384024577572964,
                       0.00082334869431643627,
                       0.00082334869431643627,
                       0.00058371735791090619,
                       0.00058371735791090619,
                       0.00058371735791090619,
                       0.00058371735791090619,
                       0.0004178187403993855,
                       0.00025806451612903221,
                       0.00025806451612903221,
                       0.00025806451612903221,
                       0.00025806451612903221,
                       0.00015975422427035328,
                       0.00015975422427035328,
                       0.00010445468509984635,
                       0.00010445468509984635,
                       0.00010445468509984635,
                       0.00010445468509984635,
                       0.00010445468509984635,
                       0.00010445468509984635,
                       0.00010445468509984635,
                       0.00010445468509984635,
                       0.00010445468509984635,
                       0.00010445468509984635,
                       0.00010445468509984635,
                       0.00010445468509984635};
    
    
    out->SetTo(1);
    double noData=0;
    // set the selection map
    //for (int j=0;j<in->Ny();j++)
    //  for (int i=0;i<in->Nx();i++)
    //    out->WrFrame(i,j,1);
    // detect signal on image
    for (int j=0;j<in->Ny();j++)
      for (int i=0;i<in->Nx();i++) {
        if (in->RdFrame(i,j)>signalCut) {
          double toomuch=(in->RdFrame(i,j)-Mean);
          int k;
          for (k=0;k<42;k++)
            if (toomuch*spread[k] < noise*0.01)
              break;
          for (int j2=MAX(j-k,0);j2<MIN(j+k+1,in->Ny());j2++)
            for (int i2=MAX(i-k,0);i2<MIN(i+k+1,in->Nx());i2++) {
              out->WrFrame(i2,j2,noData);
            }
        }
      }

}
#endif

gsl_vector* FitLine(ImageSimple* in, int j, double** model, unsigned int size, double cut) {

  // Y = XA
  if (!in->Variance())
    print_error("FitLine no variance in frame");

  // allocate some space
  gsl_matrix * XX = gsl_matrix_alloc(size,size);
  gsl_vector * XY = gsl_vector_alloc(size);
  gsl_vector * A = gsl_vector_alloc(size);
  // Y=XA
  vector<ImageSimple*>::iterator iter, iter2;

  int nbad;
  unsigned int niter,niter2;
  do {
    gsl_matrix_set_zero(XX);
    gsl_vector_set_zero(XY);
    
    for (niter = 0;niter<size;niter++){
      for (int i=0;i<in->Nx();i++) {
        double var = in->Variance()->RdFrame(i,j);
        if (var>sqrt(ut_big_value))
          var=ut_big_value;
        else
          var=1;
        (*gsl_vector_ptr(XY,niter))+=model[niter][i]*in->RdFrame(i,j)/var;
        }
      for (niter2=niter;niter2<size;niter2++){
        for (int i=0;i<in->Nx();i++) {
          double var = in->Variance()->RdFrame(i,j);
          if (var>sqrt(ut_big_value))
            var=ut_big_value;
          else
            var=1;
          
          (*gsl_matrix_ptr(XX,niter,niter2))+=model[niter][i]*model[niter2][i]/var;
          }  
      }
    }
    for (unsigned int n1=0;n1<size;n1++) {
      for (unsigned int n2=n1+1;n2<size;n2++) {
        gsl_matrix_set(XX,n2,n1,gsl_matrix_get(XX,n1,n2));
      } 
    }
    
    gsl_linalg_cholesky_decomp (XX);
    gsl_linalg_cholesky_solve (XX,XY,A);
    
    //for (unsigned int n=0;n<size;n++)
    //print_msg("Coefficient %f\n",gsl_vector_get(A,0));
    // removing far
    nbad=0;
    for (int i=0;i<in->Nx();i++) {
      double variance = in->Variance()->RdFrame(i,j);
      double delta=in->RdFrame(i,j);
      for (niter=0;niter<size;niter++)
        delta -= gsl_vector_get(A,niter)*model[niter][i];
      if (delta > cut*sqrt(variance)) {
        in->Variance()->WrFrame(i,j,ut_big_value);
        nbad++;
      }
    }
    //print_msg("removed %d pixels",nbad);
  }
  while(nbad>0);
  gsl_matrix_free(XX);
  gsl_vector_free(XY);
  
  return A;

}

gsl_vector* FitLines(ImageSimple* in, int win, double** model, double cut) {

  // Y = XA
  // in this model, the lines are free, but with pattern 1
  // only the pattern 0 is common to all lines
  if (!in->Variance())
    print_error("FitLine no variance in frame");

  int matsize=in->Ny()/win+1;
  // 0 is the special pattern and should be in

  // allocate some space
  gsl_matrix * XX = gsl_matrix_alloc(matsize,matsize);
  gsl_vector * XY = gsl_vector_alloc(matsize);
  gsl_vector * A = gsl_vector_alloc(matsize);
  // Y=XA
  vector<ImageSimple*>::iterator iter, iter2;

  int nbad;
  do {
    gsl_matrix_set_zero(XX);
    gsl_vector_set_zero(XY);

    for (int j=0;j<in->Ny();j+=win) {
      for (int i=0;i<in->Nx();i+=win) {
        double var = in->Variance()->RdFrame(i,j);
        if (var>sqrt(ut_big_value))
          var=ut_big_value;
        else
          var=1;
        (*gsl_vector_ptr(XY,0))+=model[0][i]*in->RdFrame(i,j)/var;
        (*gsl_vector_ptr(XY,j/win+1))+=model[1][i]*in->RdFrame(i,j)/var;
        
        // now the XX only (0/x correlations)
        (*gsl_matrix_ptr(XX,0,0))+=model[0][i]*model[0][i]/var;
        (*gsl_matrix_ptr(XX,0,j/win+1))+=model[0][i]*model[1][i]/var;
        (*gsl_matrix_ptr(XX,j/win+1,j/win+1))+=model[1][i]*model[1][i]/var;
      }
    }
    for (int n1=0;n1<matsize;n1++) {
      for (int n2=n1+1;n2<matsize;n2++) {
        gsl_matrix_set(XX,n2,n1,gsl_matrix_get(XX,n1,n2));
      } 
    }
    
    gsl_linalg_cholesky_decomp (XX);
    gsl_linalg_cholesky_solve (XX,XY,A);
    
    //for (unsigned int n=0;n<size;n++)
    if (DEBUG)
      print_msg("Coefficient %f",gsl_vector_get(A,0));
    // removing far
    nbad=0;
    int ibad=-1,jbad=-1;
    double maxbad=-ut_big_value;
    for (int j=0;j<in->Ny();j+=win)
      for (int i=0;i<in->Nx();i+=win) {
        double variance = in->Variance()->RdFrame(i,j);
        double delta=in->RdFrame(i,j);
        delta -= gsl_vector_get(A,0)*model[0][i];
        delta -= gsl_vector_get(A,j/win+1)*model[1][i];
        // on purpose we forget vector 1 contibution ...
        if (delta > cut*sqrt(variance) && in->RdFrame(i,j)-gsl_vector_get(A,0)*model[0][i]>maxbad) {
          maxbad=in->RdFrame(i,j)-gsl_vector_get(A,0)*model[0][i];
          ibad=i;
          jbad=j;
          nbad++;
        }
    }
    if (nbad>0) {
      in->Variance()->WrFrame(ibad,jbad,ut_big_value);
    }

    //print_msg("removed %d pixels",nbad);
  }
  while(nbad>0);
  gsl_matrix_free(XX);
  gsl_vector_free(XY);
  
  return A;

}


gsl_vector* FitImage(ImageSimple* in, int win, ImageSimple *model, double cut) {

  // Y = XA
  // in this model, the lines are free, but with pattern 1
  // only the pattern 0 is common to all lines
  if (!in->Variance())
    print_error("FitLine no variance in frame");

  int matsize=1+in->Ny()/win;
  // 0 is the model coefficient, 1 is the constant coefficient

  // allocate some space
  gsl_matrix * XX = gsl_matrix_alloc(matsize,matsize);
  gsl_vector * XY = gsl_vector_alloc(matsize);
  gsl_vector * A = gsl_vector_alloc(matsize);
  // Y=XA
  vector<ImageSimple*>::iterator iter, iter2;

  int nbad,lcount=0;
  do {
    gsl_matrix_set_zero(XX);
    gsl_vector_set_zero(XY);

    for (int j=0;j<in->Ny();j+=win) {
      for (int i=0;i<in->Nx();i+=win) {
        double var = in->Variance()->RdFrame(i,j);
        if (var>sqrt(ut_big_value))
          var=ut_big_value;
        else
          var=1;
        (*gsl_vector_ptr(XY,0))+=model->RdFrame(i,j)*in->RdFrame(i,j)/var;
        (*gsl_vector_ptr(XY,1+j/win))+=1.0*in->RdFrame(i,j)/var;
        
        // now the XX only (0/x correlations)
        (*gsl_matrix_ptr(XX,0,0))+=model->RdFrame(i,j)*model->RdFrame(i,j)/var;
        (*gsl_matrix_ptr(XX,0,1+j/win))+=model->RdFrame(i,j)*1.0/var;
        (*gsl_matrix_ptr(XX,1+j/win,1+j/win))+=1.0/var;
      }
    }
    for (int n1=0;n1<matsize;n1++) {
      for (int n2=n1+1;n2<matsize;n2++) {
        gsl_matrix_set(XX,n2,n1,gsl_matrix_get(XX,n1,n2));
      } 
    }
    
    gsl_linalg_cholesky_decomp (XX);
    gsl_linalg_cholesky_solve (XX,XY,A);
    
    //for (unsigned int n=0;n<size;n++)
    if (DEBUG) {
      if (!lcount%10)
        print_msg("Coefficient %f",gsl_vector_get(A,0));
      lcount++;
    }
    
    
    // removing far
    nbad=0;
    int ibad=-1,jbad=-1;
    double maxbad=-ut_big_value;
    for (int j=0;j<in->Ny();j+=win)
      for (int i=0;i<in->Nx();i+=win) {
        double variance = in->Variance()->RdFrame(i,j);
        double delta=in->RdFrame(i,j);
        delta -= gsl_vector_get(A,0)*model->RdFrame(i,j);
        delta -= gsl_vector_get(A,1+j/win)*1.0;
        if (delta > cut*sqrt(variance) && delta>maxbad) {
          maxbad=delta;
          ibad=i;
          jbad=j;
          nbad++;
        }
    }
    if (nbad>0) {
      in->Variance()->WrFrame(ibad,jbad,ut_big_value);
    }

    //print_msg("removed %d pixels",nbad);
  }
  while(nbad>0);
  gsl_matrix_free(XX);
  gsl_vector_free(XY);
  
  return A;

}



/* ----- main ---------------------------------------- */

int main(int argc, char **argv) {

  char **argval, **arglabel;
  char inName[lg_name+1],outName[lg_name+1],tmpName[lg_name+1];
  sprintf(tmpName,"tmp.fits");

  set_arglist("-in none -out none -win 128 -dark null");
  init_session(argv,argc,&arglabel,&argval);

  //char inName[lg_name+1],outName[lg_name+1],refName[lg_name+1];

  ImageSnifs* darkRef=0,*filtredDark=0;
  CatOrFile inCat(argval[0]);
  CatOrFile outCat(argval[1]);
  int win=atoi(argval[2]);
  if (is_set(argval[3])) {
    darkRef=new ImageSnifs(argval[3],"I");
    filtredDark=new ImageSnifs(*darkRef,tmpName,0);
  }
  

  //ImageSnifs* mask=new ImageSnifs(argval[2],"I");

  while (inCat.NextFile(inName) && outCat.NextFile(outName)) {
    ImageSnifs *in = new ImageSnifs(inName,"I");
    ImageSnifs *out = new ImageSnifs(*in,outName,0,1);
    // first guess zeBackground

    /*
    Section S(1,in->Nx(),1,in->Ny());
    ImageAnalyser ana(out,&S);
    double Mean,Rms;
    ana.SigmaClippedInfoVarKnown(3,&Mean,&Rms);
    printf("Image info %f +/- %f\n",Mean,Rms);
    


    double rdNoise=in->RdNoise();
    double sigmaCut=3.0;
    double noise=sqrt(rdNoise*rdNoise + Mean);
    double signalCut=Mean+sigmaCut*noise;
    */
  
    // set the selection map
    //for (int j=0;j<in->Ny();j++)
    //  for (int i=0;i<in->Nx();i++)
    //    out->WrFrame(i,j,1);
    // detect signal on image

    // oulah !!!
    /*
    in->Variance()->Mask(in,signalCut,3,ut_big_value);
    */
    ImageFilterSigmaClip * F = new ImageFilterSigmaClip(win,win,ImageFilter::kPixelize);
    F->SetSigma(3.0);
    F->SetInputImage(in);
    F->SetOutputImage(out);
    ((ImageFilter*)F)->Filter();
    

    double *line[2];
    line[0]=new double[2048];
    line[1]=new double[2048];
    for (int i=0;i<2048;i++) {
      line[0][i]=(i/1024)-0.5;
      line[1][i]=1;
    }
    

    gsl_vector* A;
    /*
    for (int j=0;j<out->Ny();j+=win) {
      A=FitLine(out,j,line,2,3.0);
      delete A;
    }
    */

    if (darkRef) {
      F->SetInputImage(darkRef);
      F->SetOutputImage(filtredDark);
      ((ImageFilter*)F)->Filter();
      A=FitImage(out,win,filtredDark,3.0);
      if (DEBUG) {
        printf("params %f ",gsl_vector_get(A,0));
        for (int j=0;j<out->Ny()/win;j++) 
          printf(" %f ",gsl_vector_get(A,j+1));
        printf("\n");
      }      

      for (int j=0;j<out->Ny();j++) {
        for (int i=0;i<out->Nx();i++) {
          //out->Variance()->WrFrame(i,j,out->Variance()->RdFrame(i-i%win,j-j%win));
          out->WrFrame(i,j,in->RdFrame(i,j)-gsl_vector_get(A,0)*darkRef->RdFrame(i,j));
          out->Variance()->WrFrame(i,j,in->Variance()->RdFrame(i,j));
        }
      }
    } else {

      A=FitLines(out,win,line,3.0);
      if (DEBUG) {
        printf("params %f ",gsl_vector_get(A,0));
        for (int j=0;j<out->Ny()/win;j++) 
          printf(" %f ",gsl_vector_get(A,j+1));
        printf("\n");
      }
      

      for (int j=0;j<out->Ny();j++) {
        for (int i=0;i<out->Nx();i++) {
          //        out->Variance()->WrFrame(i,j,out->Variance()->RdFrame(i-i%win,j-j%win));
          out->WrFrame(i,j,in->RdFrame(i,j)-gsl_vector_get(A,0)*line[0][i]);
          out->Variance()->WrFrame(i,j,in->Variance()->RdFrame(i,j));
        }
      }
    }
    


    out->WrDesc("DARKSTEP",DOUBLE,1,gsl_vector_ptr(A,0));
    //delete A;    

    /*
    ana.SetImage(out);
    ana.SigmaClippedInfoVarKnown(3,&Mean,&Rms);
    printf("Image info %f +/- %f\n",Mean,Rms);
    noise=sqrt(rdNoise*rdNoise + Mean)/win;
    signalCut=Mean+sigmaCut*noise;
    in->Variance()->Mask(out,signalCut,0,ut_big_value);
    */

    //out->SetTo(1);
    //out->Mask(in->Variance(),ut_big_value/10,0,0);
    //out->Mask(in->Variance(),ut_big_value/10,0,ut_big_value);
    

    /*
    ImageFilterSigmaClip * F2 = new ImageFilterSigmaClip(64,64,ImageFilter::kPixelize);
    F2->SetSigma(3.0);
    F2->SetInputImage(in);
    F2->SetOutputImage(out);
    ((ImageFilter*)F2)->Filter();

    
    ana.SetImage(out);
    ana.SigmaClippedInfoVarKnown(3,&Mean,&Rms);
    printf("Image info %f +/- %f\n",Mean,Rms);
    noise=sqrt(rdNoise*rdNoise + Mean);
    signalCut=Mean+sigmaCut*noise;
    in->Variance()->Mask(in,signalCut,3,ut_big_value);
    printf("Image info %f +/- %f\n",Mean,Rms);    
    */

    //delete F;
    //delete F2;
    delete in;
    delete out;


  }
  if (darkRef) {
    delete darkRef;
    delete filtredDark;
  }
  
  
   
  exit_session(0);
  
}
