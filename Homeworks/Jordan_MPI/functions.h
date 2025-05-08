#include <iostream>

#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define ABS(a) ((a) > 0 ? (a) : (-a))

static double EPS = 1e-16;



// class Args
// {
// public:
//     // Данные для отладки
// 	int p = 0; // количество потоков
// 	int pi = 0; // номер потока
//     MPI_Comm com = MPI_COMM_WORLD;
// 	const char * name = nullptr;
// 	double error_flag = 0;

//     // Данные для вычислений
// 	int n = 0;  // размер матрицы
// 	int m = 0;  // размер блока
// 	int l = 0;  // размер последнего блока
// 	int k = 0;  // количество блоков размера m в строке
// 	int s = 0;  // номер формулы
// 	int r = 0;  // размер выводимой
//     double r1 = 0; // невязка 1
//     double r2 = 0; // невязка 2
//     double t1 = 0; // time 1
//     double t2 = 0; // time 2
//     FILE *file;
//     double *block = nullptr;
//     double *matrix = nullptr;
//     double *inversed_matrix = nullptr;
//     double *thread_part = nullptr;
//     double *norm = nullptr;
// };

void *thread_func(void *args);
int find_diff(double *matrix, double *inversed_matrix, double *block, double *norm, FILE*file, int n, int m, int s, double &r1, double &r2, int p, int pi);
int fill_matrix(double *matrix, int n, int s);
int read_matrix_from_file(double *matrix, int n, FILE *file);
void unit_matrix(double *matrix, int n);
void rows_permutation(double *A, double *block1, double *block2, int n,
                      int m, int k, int l, int i1, int i2,
                      int begin);
void print_matrix(double *matrix, int n, int r);
void print_matrix_l_x_n(double *matrix, int l, int n);
void get_block(double *a, double *block, int n, int m, int k, int l,
               int i, int j);
void put_block(double *a, double *block, int n, int m, int k, int l,
               int i, int j);
void matrix_sum(double *A, double *B, double *C, int n, int m);
void matrix_subtr(double *A, double *B, int n, int m);
double get_norm(double *matrix, int m);
int get_inverse_matrix(double *A, double *B, int m);
// void mult(double *a, double *b, double *c, int n, int m);
double mult_sub_norm_p(double *a, double *b, double *pc, double *norm, int n, int m, int p, int pi);
void matrix_multiply(const double *A, const double *B, double *C, int p,
                     int q, int r);
void zero_matrix(double *matrix, int n, int m);

void zero_matrix_p(double *matrix, int n, int m, int p, int pi);
void unit_matrix_p(double *matrix, int n, int m, int p, int pi);
int fill_matrix_p(double *matrix, int n, int s, int m, int p, int pi);
double get_norm_p(double *matrix, int n, int m, int p, int pi);
int process_args(Args *a);
void rows_permutation_p(double *A, double *block1, double *block2, int n,
                      int m, int k, int l, int i1, int i2,
                      int begin, int p, int pi);
void *thread_func(void *args);
double get_time();

int l2g (int n, int m, int p, int pi, int i_loc);
int g2l (int n, int m, int p, int pi, int i_glob);
int get_max_cols(int n, int m, int p);
int get_cols(int n, int m, int p, int pi);
int get_loc_cols(int n, int m, int p, int pi);
int get_k( int n, int m, int p, int i_glob);
int function(int s, int n, int i, int j);
void initmatrix( double *a, int n, int m, int p, int pi, int s);
int read_matrix( double *a, int n, int m, int p, int pi, const char *name, double *buf,	MPI_Comm com);
int read_array(FILE *fp, double *a, int len);