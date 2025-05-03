#include <iostream>

#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define ABS(a) ((a) > 0 ? (a) : (-a))

static double EPS = 1e-16;


enum class io_status
{
    bad_allocation,
	error_open,
	error_read,
    irreversible,
	success,
	undef,
	unknown_formula
};


enum class reduce
{
    abs_max,
    abs_min,
    abs_max_first,
    abs_min_first,
    max,
    min,
    mult,
	sum,
    synchronize
};


class Args
{
public:
    // Данные для отладки
	size_t p = 0; // количество потоков
	size_t pi = 0; // номер потока
	pthread_t tid = 1;
	const char * name = nullptr;
	io_status error_type = io_status::undef;
	double error_flag = 0;

    // Данные для вычислений
	size_t n = 0;  // размер матрицы
	size_t m = 0;  // размер блока
	size_t l = 0;  // размер последнего блока
	size_t k = 0;  // количество блоков размера m в строке
	size_t s = 0;  // номер формулы
	size_t r = 0;  // размер выводимой
    double r1 = 0; // невязка 1
    double r2 = 0; // невязка 2
    double t1 = 0; // time 1
    double t2 = 0; // time 2
    FILE *file;
    double *block = nullptr;
    double *matrix = nullptr;
    double *inversed_matrix = nullptr;
    double *thread_part = nullptr;
    double *norm = nullptr;
};

void synchronize(int p, double *a, int n, reduce reduce_type);
void *thread_func(void *args);
int find_diff(double *matrix, double *inversed_matrix, double *block, double *norm, FILE*file, int n, int m, int s, double &r1, double &r2, size_t p, size_t pi);
int fill_matrix(double *matrix, size_t n, size_t s);
int read_matrix_from_file(double *matrix, size_t n, FILE *file);
void unit_matrix(double *matrix, size_t n);
void rows_permutation(double *A, double *block1, double *block2, size_t n,
                      size_t m, size_t k, size_t l, size_t i1, size_t i2,
                      size_t begin);
void print_matrix(double *matrix, size_t n, size_t r);
void print_matrix_l_x_n(double *matrix, size_t l, size_t n);
void get_block(double *a, double *block, size_t n, size_t m, size_t k, size_t l,
               size_t i, size_t j);
void put_block(double *a, double *block, size_t n, size_t m, size_t k, size_t l,
               size_t i, size_t j);
void matrix_sum(double *A, double *B, double *C, size_t n, size_t m);
void matrix_subtr(double *A, double *B, size_t n, size_t m);
double get_norm(double *matrix, size_t m);
int get_inverse_matrix(double *A, double *B, size_t m);
// void mult(double *a, double *b, double *c, size_t n, size_t m);
double mult_sub_norm_p(double *a, double *b, double *pc, double *norm, size_t n, size_t m, size_t p, size_t pi);
void matrix_multiply(const double *A, const double *B, double *C, size_t p,
                     size_t q, size_t r);
void zero_matrix(double *matrix, size_t n, size_t m);

void zero_matrix_p(double *matrix, size_t n, size_t m, size_t p, size_t pi);
void unit_matrix_p(double *matrix, size_t n, size_t m, size_t p, size_t pi);
int fill_matrix_p(double *matrix, size_t n, size_t s, size_t m, size_t p, size_t pi);
double get_norm_p(double *matrix, size_t n, size_t m, size_t p, size_t pi);
int process_args(Args *a);
void rows_permutation_p(double *A, double *block1, double *block2, size_t n,
                      size_t m, size_t k, size_t l, size_t i1, size_t i2,
                      size_t begin, size_t p, size_t pi);
void *thread_func(void *args);
double get_time();
int l2g (int n, int m, int p, int pi, int i_loc);
int g2l ( int n, int m, int p, int pi, int i_glob
);
int get_max_cols( int n, int m, int p
);
int get_cols( int n, int m, int p, int pi
);
int get_loc_cols(
    int n,
    int m,
    int p,
    int pi
);
int get_k( int n, int m, int p, int i_glob
);
int function(
    int s,
    int n,
    int i,
    int j
);
void initmatrix( double *a, int n, int m, int p, int pi, int s
);
int read_matrix( double *a, int n, int m, int p, int pi, const char *name, double *buf, // буффер - блочная строка n * m
	MPI_Comm com
);
int read_array(FILE *fp, double *a, int len);