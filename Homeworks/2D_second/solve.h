#include <stdio.h>
#include <math.h>

double f_0 (double /*x*/, double /*y*/);
double f_1 (double x, double /*y*/);
double f_2 (double /*x*/, double y);
double f_3 (double x, double y);
double f_4 (double x, double y);
double f_5 (double x, double y);
double f_6 (double x, double y);
double f_7 (double x, double y);

void matrix_mult_vector_msr (int n, double *A, int *I, double *x, double *y, int p, int k);
  
int minimal_errors_msr_matrix (int n, double *A, int *I, double *b,
                               double *x, /* in-out */
                               double *r, double *u, double *v,
                               double eps, int max_it, double omega, int p, int k);

int minimal_errors_msr_matrix_full (int n, double *A, int *I, double *b,
                                   double *x, /* in-out */
                                   double *r, double *u, double *v, double eps, 
                                   int max_it, int max_step, double omega, 
                                   int p, int k);

void thread_rows (int n, int p, int k, int &i1, int &i2);

void apply_preconditional_msr_matrix (int n, double *A, int *I, double *v,
                                      double *r, double omega, int p, int k);

double scalar_product (int n, double *x, double *y, int p, int k);
void mult_sub_vector (int n, double *x, double *y, double alpha, int p, int k);

//static double *results = nullptr;
double max_f_ab (double a, double b, double c, double d, int func_id, double (*f) (double, double));

int init_reduce_sum (int p);
int delete_reduce_sum ();

double reduce_sum_det (int p, int k, double s);


void init_b_fast (int n_x, int n_y, double a, double b, double c, double d, double (*f) (double, double y), double *B);


double F_ij (int nx, int ny, double hx, double hy, double b, double c, int i, int j, double (*f) (double x, double y));

// ###############################################################
// ###############################################################
// ###############################################################

void solve_01 (int n, double x[], double f_x[], double alpha[], double z[]);

double Pf_01 (double x, double a, double b, int n, double x_mas[], double alpha[]);

void solve_02 (int n, double x[], double f_x[], double mas_4n[], double dd0, double ddn);

int bin_search_to_add (double x, double * a, int n);
double Pf_02 (double x, double a, double b, int n, double x_mas[], double mas_4n[]);


// ###############################################################
// ###############################################################
// ###############################################################

void init_b (int n_x, int n_y, double a, double b, double c, double d, double (*f) (double, double), double *B, int p, int k, int disturbance, double max_abs_f);

double f_dist (int nx, int ny, double a, double b, double c, double d, double i, double j, double (*f) (double x, double y), int disturbance, double max_abs_f);

double F_ij (int nx, int ny, double a, double b, double c, double d, int i, int j, double (*f) (double x, double y), int disturbance, double max_abs_f);

double F_ij_fast (int nx, int ny, double a, double b, double c, double d, int i, int j, double (*f) (double x, double y));

void print_B (int nx, int ny, double *B, int pr);

// ###############################################################
// ###############################################################
// ###############################################################

int get_len_msr (int n_x, int n_y);

int IA_ij (int n_x, int n_y, double h_x, double h_y, int i, int j, int is, int js, int s, int *I, double *A);

// ###############################################################
// ###############################################################
// ###############################################################

int get_len_msr (int n_x, int n_y);
void ij2l (int n_x, int /*n_y*/, int i, int j, int &l);
void l2ij (int n_x, int /*n_y*/, int &i, int &j, int l);

int get_off_diag (int nx, int ny, double hx, double hy, int i, int j, int *I = nullptr, double *A = nullptr);

//возвращает количество внедиагональных точек (i, j)
int get_len_nss_off_diag (int nx, int ny);

void get_diag (int nx, int ny, double hx, double hy, int i, int j, int /*I*/, double *A);

void fill_I (int nx, int ny, double hx, double hy, int *I);

int fill_IA (int nx, int ny, double hx, double hy, int *I, double *A, int p, int k);

void print_MSR_mat (int nx, int ny, int *I, double *A);
void print_MSR (int nx, int ny, int *I, double *A, int pr);

void print_mas (int n, double *x, int pr);
//void IA_ij (int n_x, int n_y, double h_x, double h_y, int i, int j, int is, int js, int s, int *I = nullptr, double *A = nullptr);

// ###############################################################
// ###############################################################
// ###############################################################

double p_f (int nx, int ny, double a, double b, double c, double d, double *x_mas, double x, double y);

double residual_C_norm_error (int nx, int ny, double a, double b, double c, double d, double *x_mas, double (*f) (double x, double y), int p, int k);

double residual_L1_norm_error (int nx, int ny, double a, double b, double c, double d, double *x_mas, double (*f) (double x, double y), int p, int k);

double residual_C_norm_diff (int nx, int ny, double a, double b, double c, double d, double *x_mas, double (*f) (double x, double y), int p, int k);

double residual_L1_norm_diff (int nx, int ny, double a, double b, double c, double d, double *x_mas, double (*f) (double x, double y), int p, int k);
