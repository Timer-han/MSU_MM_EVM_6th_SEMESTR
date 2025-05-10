#include <iostream>

#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define ABS(a) ((a) > 0 ? (a) : (-a))

void get_block(double *a, double *block, int n, int cols, int m, int k, int l,
               int i, int j, int p, int pi);
void put_block(double *a, double *block, int n, int cols, int m, int k, int l,
               int i, int j, int p, int pi);


void print_matrix(double *matrix, int n, int r);
void print_matrix_l_x_n(double *matrix, int l, int n);


void unit_matrix(double *matrix, int n);
void unit_matrix_p(double *matrix, int n, int m, int p, int pi);
void zero_matrix(double *matrix, int n, int m);
void zero_matrix_p(double *matrix, int n, int m, int p, int pi);
int fill_matrix_p(double *matrix, int n, int s, int m, int p, int pi);


void matrix_multiply(const double *A, const double *B, double *C, int p,
                     int q, int r);
void matrix_sum(double *A, double *B, double *C, int n, int m);
void matrix_subtr(double *A, double *B, int n, int m);


double get_norm(double *matrix, int m);
double get_norm_p(double *matrix, int n, int m, int p, int pi);
int get_inverse_matrix(double *A, double *B, int m);


void rows_permutation_p(double *A, double *block1, double *block2, int n, int cols,
    int m, int k, int l, int i1, int i2,
    int begin, int p, int pi);

double get_time();