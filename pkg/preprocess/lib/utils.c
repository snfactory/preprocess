/*
This files contains some generic utilities
 */


#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "utils.h"
#include <gsl/gsl_statistics.h>
#include "IFU_io.h"

static int ascending(const void *a, const void*b) {
  if (*(double*)a==*(double*)b) return 0;  
  return (*(double*)a > *(double*)b ) ? 1 : -1;
}

/* ----- ut_median ---------------------------------------- */
double ut_median(double* values, int n)
     /* median of unsorted values */
{
  qsort(values,n,sizeof(double),&ascending);
  /* alternatively, I could have used  gsl_stats_median_from_sorted_data 
     but it is fucking the flies ! */
  /* i saw also something usable in mathlib : the names might conflict !*/
  if (n%2) 
    return (values[n/2-1]+values[n/2])/2;
  else 
    return values[n/2];
  
}

/* ----- ut_min ---------------------------------------- */
double ut_min(double* values, int n)
     /* min of array */
{
  double min = values[0];
  int i;
  
  for (i=1;i<n;i++) {
    if (values[i]<min)
      min=values[i];
  }
  return min;
}

//---------------------------------- ut_trunc_sigma_known ------------------
void ut_trunc_sigma_known(double** values, int * nitems, double sigma, double sigcut){
     /* trunc the data set, removing any value outside +/- sigcut sigma
       if sigcut<0, does noting .
       it works iteratively by removing too far away values. => slow if
       the sigcut is too low !
       
       values : pointer to array of values - modified
       nitems : number of values -> modified
       sigcut : if <0, does nothing

       Beware : the number of remaining items might be 0 !

     */
  int ok=0;

  if (sigcut<0) 
    return;
  
  qsort(*values,*nitems,sizeof(double),ascending);

  while (!ok && *nitems>1) {
    
    double mean=ut_mean(*values,*nitems);
    double valinf = (mean-(*values)[0])/sigma;
    double valsup = ((*values)[*nitems-1]-mean)/sigma;
    if (valsup>sigcut && valsup>=valinf) {
      double bound = MAX(sigcut,valinf);
      /* remove more than one for efficiency */
      while(valsup>=bound && *nitems>1) {
        (*nitems)--;
        valsup = ((*values)[*nitems-1]-mean)/sigma;
      }
    } else if (valinf>sigcut && valinf>=valsup) {
      double bound = MAX(sigcut,valinf);
      while(valinf>=bound && *nitems>1 ){
        (*nitems)--;
        (*values)++;
        valinf = (mean-(*values)[0])/sigma;
      }
    } else {
      ok=1;
    }
  }
}


//---------------------------------- ut_trunc_sigma_unknown ------------------
void ut_trunc_sigma_unknown(double** values, int * nitems, double sigcut){
     /* trunc the data set, removing any value outside +/- sigcut sigma
       if sigcut<0, does noting .
       it works iteratively by removing too far away values. => slow if
       the sigcut is too low !
       Beware : sigma is computed from the distribution itself.
       Points far away will push the sigma towards greater values.
       
       values : pointer to array of values - modified
       nitems : number of values -> modified
       sigcut : if <0, does nothing

       Beware : the number of remaining items might be 0 !

     */
  int ok=0;

  if (sigcut<0) 
    return;
  
  qsort(*values,*nitems,sizeof(double),ascending);

  while (!ok && *nitems>1) {
    double mean,sigma,valinf,valsup;

    mean = ut_mean((*values),*nitems);
    sigma = gsl_stats_sd_m(*values,1,*nitems,mean);
    valinf = (mean-(*values)[0])/sigma;
    valsup = ((*values)[*nitems-1]-mean)/sigma;
    if (valsup>sigcut && valsup>=valinf) {
      (*nitems)--;
    } else if (valinf>sigcut && valinf>=valsup) {
      (*nitems)--;
      (*values)++;
    } else {
      ok=1;
    }
  }
}

/*-------------------- ut_fraction_sigcut ----------------------*/
double ut_fraction_sigcut(double fraction) {
  /* given a distribution with fraction outliers, at infinite 
   distance with respect to normal values, 
  returns the standard deviation of the outliers, when the
  RMS is computed with it included.*/

  /* the computation is performed with fraction values at 1 and 
   1-fraction at 0 */
  double sigcut,rms;
  /* mean = fraction */
  rms = sqrt(fraction*fraction*(1-fraction) + (1-fraction)*(1-fraction)*fraction);
  
  sigcut = (1-fraction)/rms;
  return sigcut;
}

/*-------------------- ut_mean ----------------------*/
double ut_mean(double* val,int n){

  int i;
  double sum = 0;
  for (i=0;i<n;i++)
    sum += val[i];
  return sum/n;
}

/*-------------------- ut_autocovariance ----------------------*/
void ut_autocovariance(double* data, double* autocorr, int n) {
/* data is of dim n, autocorr of dim n */
  int decal,i;
  double mean;
  //  double var;
  double *mdata = malloc(sizeof(double)*n);

  /* compute the mean */
  mean = ut_mean(data,n);

  /* substract it (n computations) */
  for (i=0;i<n;i++)
    mdata[i] = data[i]-mean;

  /* compute variance and what is needed for decal=0*/
  decal=0;
  autocorr[0]=0;
  for (i=0;i<n;i++) {
    autocorr[0]+=mdata[i]*mdata[i];
  }
  autocorr[0] /= n;
  //var = autocorr[0];
  //autocorr[0] = 1.0;

  /* n^2 part */
  for (decal=1; decal<n; decal ++) {
    autocorr[decal]=0;    
    for (i=0;i<n-decal;i++) {
      autocorr[decal]+=mdata[i]*mdata[i+decal];
    }
    autocorr[decal] /= (n-decal);
    //    autocorr[decal] /= var;
  }
  free(mdata);
  
}


/*-------------------- ut_varname_from_imname ----------------------*/
void ut_varname_from_imname(char* imname,char* varname){
  /* returns the name of the variance image from the name of 
   the image */

  /* check standard name */
  char extname[lg_name+1];
  extname[0]='\0';
  varname[0]='\0';
  if (strstr(imname,"[chip00]"))
    strcpy(extname,"[var00]");
  if (strstr(imname,"[chip01]"))
    strcpy(extname,"[var01]");
  if (strstr(imname,"[image]"))
    strcpy(extname,"[variance]");

  if (extname[0]) {
    /* builds the variance name */
    char *ptName;
    strcpy(varname,imname);
    ptName = strchr(varname,'[');
    strcpy(ptName,extname);
  }
}

/*-------------------- ut_open_check_name ----------------------*/
char* ut_open_check_name(char* name){
  /* if name is a name without extension, tries to see if a 
     [image] exists */
  
  if (!strchr(name,'[')) {
    char new_name[lg_name+1];
    sprintf(new_name,"%s[image]",name);
    if (exist(new_name))
      strcpy(name,new_name);
  }
  return name;
}

/*-------------------- ut_open_check_name ----------------------*/
char* ut_create_check_name(char* name){
  /* if name is a name without extension, appends an [image] 
     by default */
  
  if (!strchr(name,'[')) {
    char new_name[lg_name+1];
    sprintf(new_name,"%s[image]",name);
    strcpy(name,new_name);
  }
  return name;
}




