#include <iostream>

#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define ABS(a) ((a) > 0 ? (a) : (-a))

// work with blocks
void get_block(double *a, double *block, int n, int cols, int m, int k, int l,
               int i, int j, int p, int pi, int use_local_j = 1,
               int bl_rows = -1, int bl_cols = -1);
void put_block(double *a, double *block, int n, int cols, int m, int k, int l,
               int i, int j, int p, int pi);
void rows_permutation_p(double *A, double *block1, double *block2, int n, int cols,
                        int m, int k, int l, int i1, int i2,
                        int begin, int p, int pi);

// reading
int read_matrix(double *a, int n, int m, int p, int pi, const char *name,
                double *buf, MPI_Comm com);
int read_array(FILE *fp, double *a, int len);

// printing
void print_matrix(double *matrix, int n, int r);
void print_matrix_l_x_n(double *matrix, int l, int n);
void print_matrix_mpi(double *a, int n, int m, int p, int pi,
    double *buf, int max_print, MPI_Comm com);
void print_array(double *a, int n, int m, int l, int max_print, int & printed_rows, int p, int rows);

// filling
void unit_matrix(double *matrix, int n);
void unit_matrix_mpi(double *matrix, int n, int m, int p, int pi);
void zero_matrix(double *matrix, int n, int m);
void zero_matrix_mpi(double *matrix, int n, int m, int p, int pi);
int fill_matrix_p(double *matrix, int n, int s, int m, int p, int pi);
double function(int s, int n, int i, int j);
void init_matrix(double *a, int n, int m, int p, int pi, int s);


// matrix operations
void matrix_multiply(const double *A, const double *B, double *C, int p,
                     int q, int r);
void matrix_sum(double *A, double *B, double *C, int n, int m);
void matrix_subtr(double *A, double *B, int n, int m);
int get_inverse_matrix(double *A, double *B, int m);

// norm calculation
double get_norm(double *matrix, int m);
double get_norm_p(double *matrix, int n, int m, int p, int pi);
double get_norm_pi(double *matrix, int n, int cols);

// time calculation
double get_time();

// index calculation
int l2g(int n, int m, int p, int pi, int i_loc);
int g2l(int n, int m, int p, int pi, int i_glob);
int get_max_cols(int n, int m, int p);
int get_bl_cols(int n, int m, int p, int pi);
int get_loc_cols(int n, int m, int p, int pi);
int get_pi(int n, int m, int p, int i_glob);
void get_column(double *matrix, double *buffer,
    int n, int m, int width,
    int p, int pi, int i);

// main functions
void init_matrix(double *a, int n, int m, int p, int pi, int s);
int mpi_calculate(
    double *matrix,          // n x (m * bl_cols)
    double *inversed_matrix, // n x (m * bl_cols)
    double *buffer,
    int n,
    int m,
    int p,
    int pi,
    MPI_Comm com);

double residual_calculate_mpi(
    double *matrix,
    double *inversed,
    int n,
    int m,
    int p,
    int pi,
    MPI_Comm com
);