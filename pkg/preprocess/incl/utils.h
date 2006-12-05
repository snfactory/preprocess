#ifndef UTILS_H 
#define UTILS_H

#ifdef __cplusplus
extern "C" {
#endif

const double ut_big_value = 1e31;
const int lg_message = 256;

#include "IFU_io.h" 

/* ----- statistical utilities ----- */
void ut_sort_ascending(double* values, int n);

double ut_median(double* values, int n);

double ut_min(double* values, int n);
  
double ut_mean(double* values, int n);

double ut_mean_var_sigmas_known(double* values, int n, double* mean, double* var);

double ut_mode(double* values, int n);

void ut_autocovariance(double* data, double* autocorr, int n);

void ut_trunc_sigma_known(double** values, int* n_init, double sigma, double sigcut);

void ut_trunc_sigma_unknown(double** values, int* n_init, double sigcut);

void ut_trunc_sigma_unknown_fast_sorted(double** values, int* n_init, double sigcut);

double ut_fraction_sigcut(double fraction);
  
/* ----- name utilities ----- */
void ut_varname_from_imname(char* imname,char* varname);

char* ut_open_check_name(char* name);

char* ut_create_check_name(char* name);
  
int ut_is_bichip_detcom(char* filename);

void ut_build_tmp_name(char* filename, char* tmp_prefix);

void ut_primary_header_name(char* full_name, char* primary_name);
  

/* ----- library patches ---------- */

int close_frame_fast(IMAGE2D *frame);

int open_frame_fast(IMAGE2D *frame, char *name, char *mode);

double juldat( int year, int month, int day, double ut);

/* -----  string utilities ----- */

int ut_parse_line(char* line, char** tokens);
  

#ifdef __cplusplus
}
#endif

#endif
