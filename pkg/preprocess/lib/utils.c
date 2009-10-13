/*
This files contains some generic utilities
 */

/* ----- System includes -----*/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/* ----- IFU and GSL includes -----*/
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_sort.h>
#include "IFU_io.h"
#include "IFU_math.h"
/* ----- Specific IFU includes ----- */
#include "data_io.h"
#ifdef FITS
#include "fitsio.h"
#endif

/* Removing MIDAS by hand - not really clean, but should allow to compile with an IFU with --disable-midas enabled*/
#undef MIDAS

#ifdef MIDAS
#include "midas_defs.h"
#include "tbldef.h"
#include "fctext.h"
#include "ldbext.h"
#include "proto_st.h"

#endif
extern IO_Format InputIO, OutputIO;

/* ----- local includes ----- */
#include "utils.h"


/* ----- static functions ---------------------------------------- */
static int ascending(const void *a, const void*b) {
  if (*(double*)a==*(double*)b) return 0;  
  return (*(double*)a > *(double*)b ) ? 1 : -1;
}

/* ----- ut_median ---------------------------------------- */
void ut_sort_ascending(double* values, int n)
     /* median of unsorted values */
{
  qsort(values,n,sizeof(double),&ascending);
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

/* ----------------------------- ut_trunc_sigma_known -------------- */
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


/* ------------------------------- ut_trunc_sigma_unknown ------------------ */
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

/* ------------------------------- ut_trunc_sigma_unknown ------------------ */
void ut_trunc_sigma_unknown_fast_sorted(double** values, int * nitems, double sigcut){
     /* trunc the data set, removing any value outside +/- sigcut sigma
       if sigcut<0, does noting .
       it works grossly, removing each time all what is outside the window, 
       which is not accurate

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
  
  /*  qsort(*values,*nitems,sizeof(double),ascending); */

  while (!ok && *nitems>1) {

    double mean,sigma,valinf,valsup;
    int nitems_old = *nitems;

    mean = ut_mean((*values),*nitems);
    sigma = gsl_stats_sd_m(*values,1,*nitems,mean);
    valinf = mean-sigma*sigcut;
    valsup = mean + sigma*sigcut;
    while ((*values)[0] < valinf && *nitems>1) {
      (*nitems)--;
      (*values)++;
    }
    while ((*values)[*nitems-1] > valsup && *nitems>1) {
      (*nitems)--;
    }
    if (nitems_old == (*nitems))
      ok=1;
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

/*-------------------- ut_mean ----------------------*/
double ut_mode(double* val,int n){
  /* This function returns the highest probability density a set of
  datas. This only works reliabely if there is only 1 maximum, 
  and in the limit where the data is sufficiently sampled. 
  No error computation is performed ...*/
  /* The search is dyadic and iterative, each step looking for the part
  where the median is the narrower. The time is O(NlnN) (+ the time to run
  indexx)*/ 
  /* This is a fast and inaccurate algorithm. For a general algorithm, see the
  Parzenwindow method - the drawback is that it depends on a window size
  a priori*/  


  size_t * index;
  double med;
  int nmin=0;
  int nmax=n; /*  out of the table */

  index = malloc(n*sizeof(size_t));
  
  /* test */

  gsl_sort_index(index,val,1,n);

  
  do { 

    /* median computation */
    if ((nmax-nmin) % 2 == 0)
      med = (val[index[ (nmin + nmax )/2 ]] + val[index[ (nmin+nmax)/2 -1 ]])/2;
    else
      med = val[index[ (nmin + nmax)/2 ]];
    
    /* The density is higher for the upper part */
    if ( fabs(med - val[index[nmin]])  > fabs (med - val[index[nmax-1]]) )
      nmin =(nmin + nmax - 1 )/2;
    else if ( fabs(med - val[index[nmin]])  < fabs (med - val[index[nmax-1]]) )
      /* rounding in the good direction  in order to include the median place */
      nmax =  ( nmin + nmax)/2 + 1;
    else /* stop in case of equality */
      nmax = nmin;
  } while (nmax-nmin > 2);
  
  return med;
}

/*-------------------- ut_mean ----------------------*/
double ut_mode2(double* val,int n){
  /* This function returns the highest probability density a set of
  datas. This only works reliabely if there is only 1 maximum, 
  and in the limit where the data is sufficiently sampled. 
  No error computation is performed ...*/
  /* The search is dyadic and iterative, each step looking for the part
  which contains most of the data with respect to the highest value. 
  The time is O(NlnN) (+ the time to run indexx)*/ 
  /* This is a fast and inaccurate algorithm. For a general algorithm, see the
  Parzenwindow method - the drawback is that it depends on a window size
  a priori*/  


  size_t * index;
  double half;
  int nmin=0;
  int nmax=n; /*  out of the table */

  index = malloc(n*sizeof(size_t));
  
  /* test */

  gsl_sort_index(index,val,1,n);

  
  do { 
    int i;
    /* half window computation */
    half=(val[index[nmin]]+val[index[nmax-1]])/2;
    for (i=nmin;i<nmax;i++) {
      if (val[index[i]]>half) 
        break;
    }
    if (i==nmax)
      return half;
    else if (i<nmin+(nmax-nmin)/2) 
      nmin=i;
    else if (i>nmin+(nmax-nmin)/2 || ((nmax-nmin)%2) )
      nmax=i;
    else {
      nmin++;
      nmax--;
    }
  } while (nmax-nmin>2);
  return half;
}

/*-------------------- ut_autocovariance ----------------------*/
void ut_autocovariance(double* data, double* autocorr, int n) {
/* data is of dim n, autocorr of dim n */
  int decal,i;
  double mean;
  /*  double var; */
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
  /* var = autocorr[0]; */
  /* autocorr[0] = 1.0; */

  /* n^2 part */
  for (decal=1; decal<n; decal ++) {
    autocorr[decal]=0;    
    for (i=0;i<n-decal;i++) {
      autocorr[decal]+=mdata[i]*mdata[i+decal];
    }
    autocorr[decal] /= (n-decal);
    /*    autocorr[decal] /= var; */
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
  if (strstr(imname,"[chip02]"))
    strcpy(extname,"[var02]");
  if (strstr(imname,"[chip03]"))
    strcpy(extname,"[var03]");
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

/*-------------------- ut_varname_from_imname ----------------------*/
int ut_is_bichip_detcom(char* filename){
  /* returns 1 if the name corresponds to a detcom image */

  char testName[lg_name+1];
  if (strstr(filename,"%d"))
    return 1;
  if (strrchr(filename,'['))
    return 0;
  // .fits file extension default
  if (strstr(filename,".fits"))
    sprintf(testName,"%s[chip00]",filename);
  else
    sprintf(testName,"%s.fits[chip00]",filename);
  if (exist(testName))
    return 1;
  else {
    if (strstr(filename,".fits"))
      sprintf(testName,"%s[CHIP00]",filename);
    else
      sprintf(testName,"%s.fits[CHIP00]",filename);
    if (exist(testName))
      return 1;
    else 
      return 0;
  }
}


/*-------------------- ut_build_tmp_name ----------------------*/
void ut_build_tmp_name(char* filename, char* tmp_prefix){
  /* builds a temporary memory file name from the image name and a
  prefix */
  char imname[lg_name+1], *pt_name;

  if (!(pt_name = strrchr(filename,'/')))
    pt_name = filename;
  else 
    pt_name = pt_name+1;
  sprintf(imname,"mem://%s_%s",tmp_prefix,pt_name);
  //  sprintf(imname,"%s_%s",tmp_prefix,pt_name);
  strcpy(filename,imname);
}

/*-------------------- ut_primary_header_name ----------------------*/
void ut_primary_header_name(char* full_name, char* primary_name){
  /* the primary header is the image name without any extension*/

  char* pt_name;
  strcpy(primary_name,full_name);

  if ((pt_name = strchr(primary_name,'[')))
    pt_name[0]='\0';

}

/* ===== IFU lib disagreement ======================================== */

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!.func                        open_frame_fast()
!
!.purp          opens a 2D frame and updates the image structure items
!.desc
! open_frame(frame,name,mode)	
!
! IMAGE2D *frame;       image structure
! char *name;           frame name
! char *mode;           open mode (Input,Ouput,IO)
!.ed
-------------------------------------------------------------------- */
int 
open_frame_fast(IMAGE2D *frame, char *name, char *mode)		
{
  char errtext[132], filename[lg_name+1];
  int status, nbaxes, iomode, int_datatype;
  float cuts[4];
  int info[5];
#ifdef IRAF
  int two_dim=2;
  int len;
#endif
#ifdef FITS
  fitsfile *fptr;
  int nbread;
  int npix;
  int group = 0;
  double pixref;
#endif

  memset(frame->ident,' ',lg_ident);
  frame->ident[lg_ident] = '\0';
  memset(frame->cunit,' ',lg_unit);
  frame->cunit[lg_unit] = '\0';
  memset(frame->history,' ',lg_hist);
  frame->history[lg_hist] = '\0';
  frame->external_info = NULL;
  frame->file_type = T_IMA2D;
  frame->data_format = InputIO.basic_io;

  strcpy(filename,name);
  first_blk(filename); 
  strcpy(frame->name,filename);
  append_ima_extension(frame->name,InputIO.basic_io);

  strcpy(filename,frame->name);

  if (!exist(filename)) { /* check if fil exists */
    status = ERR_OPEN;
    sprintf(errtext,"open_frame: frame %s",filename);
    Handle_Error(errtext,status);
    return(status);
  }

  switch(mode[0]) {
  case 'I' : 
    if (mode[1] == 'O')
      frame->iomode = (int)IO_MODE;
    else
      frame->iomode = (int)I_MODE;
    break;
  case 'O' : frame->iomode = (int)O_MODE;
    break;
  default  : frame->iomode = (int)I_MODE;
    break;
  }
	
  iomode = get_iomode_code(InputIO.basic_io,frame->iomode);

  switch (InputIO.basic_io) {

#ifdef MIDAS
  case MIDAS_FORMAT :
    status = SCFINF(filename,2,info);  
    if (status == 0) {
      status = SCIGET(filename, info[1], iomode, F_IMA_TYPE, 2, 
                      &nbaxes, &(frame->nx), &(frame->startx), &(frame->stepx),
                      frame->ident, frame->cunit, (char **)(&(frame->data)), 
                      &(frame->imno));
      frame->data_type = info[1];
      frame->data_type = decode_datatype(InputIO.basic_io,frame->data_type);

      if (nbaxes!=2) /* We open a spectrum like an image, and that's not good */
        status = ERR_OPEN; 

    }
    break;
#endif
#ifdef IRAF
  case IRAF_FORMAT :
  case STSDAS_FORMAT :
    len = strlen(filename);
    uimopn(filename,&iomode,&(frame->imno),&status,len);
    if (status != 0) 
      break;
    uimgid(&(frame->imno),&int_datatype,&two_dim,&(frame->nx),&status);
    frame->data_type = decode_datatype(InputIO.basic_io,(short)(int_datatype));
    if (status != 0)
      break;
    alloc_frame_mem(frame, datatype);
    switch(frame->data_type) {
    case SHORT :
      uigs2s(&(frame->imno),&one,&(frame->nx),&one,&(frame->ny),
             frame->data.s_data,&status);
      break;
    case INT :
    case LONG :
      uigs2l(&(frame->imno),&one,&(frame->nx),&one,&(frame->ny),
             frame->data.l_data,&status);
      break;
    case FLOAT :
      uigs2r(&(frame->imno),&one,&(frame->nx),&one,&(frame->ny),
             frame->data.f_data,&status);
      break;
    case DOUBLE :
      uigs2d(&(frame->imno),&one,&(frame->nx),&one,&(frame->ny),
             frame->data.d_data,&status);
      break;
    }
    disable_user_warnings();
    RD_desc(frame,"IDENT",CHAR,lg_ident,frame->ident);
    restore_user_warnings();
    break;
#endif
#ifdef FITS
  case FITS_A_FORMAT :
  case FITS_B_FORMAT :
    status =0;
    if (fits_open_file(&fptr,filename,iomode,&status)) {
      status = ERR_ACCESS; break;
    }
    frame->external_info = (void *)fptr;
    if (fits_read_key(fptr, TINT,"NAXIS", &nbaxes,NULL, &status)) {
      status = ERR_READ; break;
    }
    if (nbaxes != 2) {
      status = ERR_IMA_HEAD; break;
    }
    if (fits_read_key(fptr, TINT, "NAXIS1",
                      &(frame->nx), NULL, &status)) {
      status = ERR_READ; break;
    }
    if (fits_read_key(fptr, TINT, "NAXIS2",
                      &(frame->ny), NULL, &status)) {
      status = ERR_READ; break;
    }
    if (status == 0) {
      pixref = 1.0;
      fits_read_key(fptr, TDOUBLE, "CRPIX1", &pixref, NULL, &status);
      if (status) { status = 0; pixref = 1; }
      fits_read_key(fptr, TDOUBLE, "CRVAL1", &(frame->startx), NULL, &status);
      if (status) { status = 0; frame->startx = (double)1; }
      fits_read_key(fptr, TDOUBLE, "CDELT1", &(frame->stepx), NULL, &status);
      if (status) { status = 0; frame->stepx = (double)1; }
      frame->startx -= (pixref-1)*frame->stepx;
      pixref = 1.0;
      fits_read_key(fptr, TDOUBLE, "CRPIX2", &pixref, NULL, &status);
      if (status) { status = 0; pixref = 1; }
      fits_read_key(fptr, TDOUBLE, "CRVAL2", &(frame->starty), NULL, &status);
      if (status) { status = 0; frame->starty = (double)1; }
      fits_read_key(fptr, TDOUBLE, "CDELT2", &(frame->stepy), NULL, &status);
      if (status) { status = 0; frame->stepy = (double)1; }
      frame->starty -= (pixref-1)*frame->stepy;
    }
    else
      break;

    int_datatype = (fptr->Fptr->tableptr)->tdatatype;
    frame->data_type = decode_datatype(InputIO.basic_io,(short)int_datatype);
    if (frame->data_type == SHORT) {
      if (fptr->Fptr->tableptr[1].tscale == 1 && fptr->Fptr->tableptr[1].tzero == 32768)
        /* unsigned short !!! */
        frame->data_type = LONG;
    }

    if (alloc_frame_mem(frame, frame->data_type) < 0) {
      fits_close_file(fptr,&status);
      status = ERR_ALLOC;
      break;
    }

    npix = frame->nx*frame->ny;
    switch (frame->data_type) {
    case SHORT :
      if (fits_read_img_sht(fptr,group,1L,npix,(short)0,
                            frame->data.s_data,&nbread,&status)) {
        status = ERR_READ; break;
      }
      break;
    case LONG :
    case INT :
      if (fits_read_img_lng(fptr,group,1L,npix,(int)0,
                            frame->data.l_data,&nbread,&status)) {
        status = ERR_READ; break;
      }
      break;
    case FLOAT :
      if (fits_read_img_flt(fptr,group,1L,npix,(float)0,
                            frame->data.f_data,&nbread,&status)) {
        status = ERR_READ; break;
      }
      break;
    case DOUBLE :
      if (fits_read_img_dbl(fptr,group,1L,npix,(double)0,
                            frame->data.d_data,&nbread,&status)) {
        status = ERR_READ; break;
      }
      break;
    }
    break;
#endif
  }

  if (status) {
    sprintf(errtext,"open_frame: frame %s",filename);
    status = get_tiger_errcode(frame->data_format,status);
    Handle_Error(errtext,status);
  }
  else {
    disable_user_warnings();
    status = RD_desc(frame,"LHCUTS",FLOAT,4,cuts);
    RD_desc(frame,"HISTORY",CHAR,lg_hist,frame->history);
    restore_user_warnings();

    frame->endx = frame->startx + (frame->nx -1)*frame->stepx;
    frame->endy = frame->starty + (frame->ny -1)*frame->stepy;
    if (status <= 0) {
      /* image_minmax is a really slow routine, and most of the time useless */
      frame->min = -ut_big_value;
      frame->max = +ut_big_value;
      /*			image_minmax(frame); */
    }
    else {
      frame->min = cuts[2];
      frame->max = cuts[3];
    }
    status = 0;
    /* parse wcs if contained in file */
    status = parse_wcs(frame);
  }
  return(status);
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!.func                       close_frame_fast()
!
!.purp             closes a currently active 2D frame 
!.desc
! close_frame(frame)	
!
! IMAGE2D *frame;       image structure
!.ed
-------------------------------------------------------------------- */
int 
close_frame_fast(IMAGE2D *frame)			/* close active frame */
{
  char  errtext[132], filename[lg_name+1];
  int   stat, int_datatype;
  float cuts[4];
#ifdef IRAF
  int one=1;
#endif
#ifdef FITS
  fitsfile *fptr;
  int npix;
#endif

  strcpy(filename,frame->name);

  if (frame->iomode == (int)I_MODE) {
    switch (frame->data_format) {
#ifdef MIDAS
    case MIDAS_FORMAT :
      stat = SCFCLO(frame->imno);
      break;
#endif
#ifdef IRAF
    case IRAF_FORMAT :
    case STSDAS_FORMAT :
      uimclo(&(frame->imno),&stat);
      break;
#endif
#ifdef FITS
    case FITS_A_FORMAT :
    case FITS_B_FORMAT :
      stat =0;
      fptr = (fitsfile *)frame->external_info;
      fits_close_file(fptr,&stat);
      free_frame_mem(frame);
      frame->external_info = NULL;
      break;
#endif
    }
    if (stat) {
      sprintf(errtext,"close_frame: frame %s",filename);
      stat = get_tiger_errcode(frame->data_format,stat);
      Handle_Error(errtext,stat);
    }
    return(stat);
  }

  /*
  if (frame->data.d_data != NULL) {
    image_minmax(frame);

    cuts[0]=(float)frame->min; cuts[2]=(float)frame->min;
    cuts[1]=(float)frame->max; cuts[3]=(float)frame->max;
    stat = WR_desc(frame,"LHCUTS",FLOAT,4,cuts);
  }
  */

  WR_history(frame, (Anyfile *)0);

  switch (frame->data_format) {
#ifdef MIDAS
  case MIDAS_FORMAT :
    stat = SCFCLO(frame->imno);
    break;
#endif
#ifdef IRAF
  case IRAF_FORMAT :
  case STSDAS_FORMAT :
    switch(frame->data_type) {
    case SHORT :
      uips2s(&(frame->imno),&one,&(frame->nx),&one,&(frame->ny),
             frame->data.s_data,&stat);
      break;
    case INT :
    case LONG :
      uips2l(&(frame->imno),&one,&(frame->nx),&one,&(frame->ny),
             frame->data.l_data,&stat);
      break;
    case FLOAT :
      uips2r(&(frame->imno),&one,&(frame->nx),&one,&(frame->ny),
             frame->data.f_data,&stat);
      break;
    case DOUBLE :
      uips2d(&(frame->imno),&one,&(frame->nx),&one,&(frame->ny),
             frame->data.d_data,&stat);
      break;
    }
    if (stat == 0)  
      uimclo(&(frame->imno),&stat);
    free_frame_mem(frame);
    break;
#endif
#ifdef FITS
  case FITS_A_FORMAT :
  case FITS_B_FORMAT :
    stat = 0;
    fptr = (fitsfile *)frame->external_info;
    if (frame->iomode != (int)I_MODE) {
      if (frame->data.d_data != NULL) {
        int_datatype = get_datatype_code(OutputIO.basic_io,frame->data_type);
        npix = frame->nx*frame->ny;
        if (fits_write_img(fptr,int_datatype,1L,npix,
                           frame->data.s_data,&stat)) {
          stat = ERR_WRIT;
        }
      }
    }
    if (! stat) {
      fits_close_file(fptr,&stat);
      stat = wcs_free(frame);
    }
    free_frame_mem(frame);
    frame->external_info = NULL;
    break;
#endif
  }
  if (stat) {
    sprintf(errtext,"close_frame: frame %s",filename);
    stat = get_tiger_errcode(frame->data_format,stat);
    Handle_Error(errtext,stat);
  } else {
    if (TK && (frame->iomode == O_MODE || frame->iomode == IO_MODE))
      {
        printf("@ N {%s}\n",filename);
      }
  }

  return(stat);
}


/* ----- juldat --------------------------------------------------------- */
/*C
C       FUNCTION JULDAT RETURNS THE JULIAN DATE AS A DOUBLE
C       PRECISION NUMBER.
C       ALGORITHM FROM "ALAMANAC FOR COMPUTERS" p.B2
C*/
double juldat( int year, int month, int day, double ut)
     /* year, month, day, ut are the UTC ones
      ut is the fractional hour */
{
        double tmp;
        tmp = 367.0*year + 0.5 + ut/24.0;
        tmp = tmp - ((7*(year+((month+9)/12)))/4) + ((275*month)/9) + day + 1721013;
        return tmp;
}

/* =====  string utilities ===== */

int ut_parse_line(char* line, char** tokens) {
  /* the token strings will contain the line elements separated by spaces */
  /* the original line will be modified (some '/0' char added)*/
  char* pt_line = line;
  int itoken=0;

  while (pt_line && pt_line[0]!='\0') {
    tokens[itoken++] = pt_line;
    pt_line = strchr(pt_line,' ');
    if (pt_line)
      while((pt_line)[0]==' ') {
        pt_line[0]='\0';
        pt_line++;
      }
  }
  return itoken;
}
