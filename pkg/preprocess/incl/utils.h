#ifndef UTILS_H 
#define UTILS_H

#ifdef __cplusplus
extern "C" {
#endif

#define MIN(x,y) (x<y ? x : y)
#define MAX(x,y) (x>y ? x : y)

/* ----- median of the values ----- */
double ut_median(double* values, int n);

double ut_min(double* values, int n);
  
double ut_mean(double* values, int n);

void ut_trunc_sigma_known(double** values, int* n_init, double sigma, double sigcut);

void ut_trunc_sigma_unknown(double** values, int* n_init, double sigcut);

double ut_fraction_sigcut(double fraction);
  
/* ----- name utilities ----- */
void ut_varname_from_imname(char* imname,char* varname);

char* ut_open_check_name(char* name);
char* ut_create_check_name(char* name);
  

#ifdef __cplusplus
}
#endif

#endif
