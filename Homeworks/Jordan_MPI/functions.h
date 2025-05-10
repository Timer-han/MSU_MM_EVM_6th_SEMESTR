#include <iostream>

#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define ABS(a) ((a) > 0 ? (a) : (-a))




double get_time(void);

/* Базовые операции с элементами блоков */
void get_block(double *a, double *block,
               int n, int cols, int m, int k, int l,
               int i, int j, int p, int pi);
void put_block(double *a, double *block,
               int n, int cols, int m, int k, int l,
               int i, int j, int p, int pi);

/* Вывод матриц */
void print_matrix(double *matrix, int n, int r);
void print_matrix_l_x_n(double *matrix, int l, int n);
void print_matrix_mpi(double *a,
                      int n, int m, int p, int pi,
                      double *buf, int max_print,
                      MPI_Comm com);

/* Инициализация и чтение матриц */
int function(int s, int n, int i, int j);
void initmatrix(double *a,
                int n, int m, int p, int pi,
                int s);
int read_matrix(double *a,
                int n, int m, int p, int pi,
                const char *name,
                double *buf,
                MPI_Comm com);
int read_matrix_from_file(double *matrix,
                          int n,
                          FILE *file);
int fill_matrix(double *matrix,
                int n, int s);

/* Простые преобразования матриц */
void unit_matrix(double *matrix, int n);
void zero_matrix(double *matrix, int n, int m);
void zero_matrix_p(double *matrix,
                   int n, int m, int p, int pi);
void get_column(double *matrix, double *buffer,
                int n, int m, int p, int pi,
                int i);

/* Блочное переставление строк */
void rows_permutation(double *A,
                      double *block1, double *block2,
                      int n, int m, int k, int l,
                      int i1, int i2,
                      int begin);
void rows_permutation_p(double *A,
                        double *block1, double *block2,
                        int n, int cols,
                        int m, int k, int l,
                        int i1, int i2,
                        int begin,
                        int p, int pi);

/* Алгебраические операции */
void matrix_multiply(const double *A,
                     const double *B,
                     double *C,
                     int p, int q, int r);
void matrix_sum(double *A,
                double *B,
                double *C,
                int n, int m);
void matrix_subtr(double *A,
                  double *B,
                  int n, int m);

/* Норма и обратная матрица */
double get_norm(double *matrix, int m);
int get_inverse_matrix(double *A,
                       double *B,
                       int m);

/* MPI–вычисления */
int mpi_calculate(double *matrix,          /* n x (m * bl_cols) */
                  double *inversed_matrix, /* n x (m * bl_cols) */
                  int n,
                  int m,
                  int p,
                  int pi,
                  MPI_Comm com);