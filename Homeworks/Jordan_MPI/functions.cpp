#include "mpi.h"
#include "functions.h"

#include <cfloat>
#include <stdlib.h>
#include <string.h>
#include <sys/sysinfo.h>
#include <sys/time.h>
#include <cmath>
#include <time.h>

static double EPS = 1e-16;

void get_block(double *a, double *block, int n, int cols, int m, int k, int l,
               int i, int j, int p, int pi)
{
    int row_count, col_count;
    int row_offset, col_offset;
    int block_row, block_col;
    int local_j = g2l(n, m, p, pi, j);

    // Определяем размеры блока
    if (i < k)
        row_count = m;
    else
        row_count = l;

    if (j < k)
        col_count = m;
    else
        col_count = l;

    // Вычисляем смещение в исходной матрице
    row_offset =       i * m;
    col_offset = local_j * m;

    // Копируем элементы блока
    for (block_row = 0; block_row < row_count; block_row++)
    {
        for (block_col = 0; block_col < col_count; block_col++)
        {
            block[block_row * col_count + block_col] =
                a[(row_offset + block_row) * cols + (col_offset + block_col)];
        }
    }
}

void put_block(double *a, double *block, int n, int cols, int m, int k, int l,
               int i, int j, int p, int pi)
{
    int row_size, col_size;
    int row_offset, col_offset;
    int block_row, block_col;
    int local_j = g2l(n, m, p, pi, j);

    // Определяем размеры блока
    if (i < k)
        row_size = m;
    else
        row_size = l;

    if (j < k)
        col_size = m;
    else
        col_size = l;

    // Вычисляем смещение в матрице
    row_offset =       i * m;
    col_offset = local_j * m;

    // Копируем элементы из блока обратно в матрицу
    for (block_row = 0; block_row < row_size; block_row++)
    {
        for (block_col = 0; block_col < col_size; block_col++)
        {
            a[(row_offset + block_row) * cols + (col_offset + block_col)] =
                block[block_row * col_size + block_col];
        }
    }
}

void print_matrix(double *matrix, int n, int r)
{
    int rows = (r > n) ? n : r;
    int cols = (r > n) ? n : r;
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            printf(" %10.3e", matrix[i * n + j]);
        }
        printf("\n");
    }
}



void unit_matrix(double *matrix, int n)
{
    memset(matrix, 0, n * n * sizeof(double));
    for (int i = 0; i < n; i++) {
        matrix[i * n + i] = 1;
    }
}

void zero_matrix(double *matrix, int n, int m)
{
    memset(matrix, 0, n * m * sizeof(double));
}

void matrix_multiply(const double *A, const double *B, double *C, int p,
                     int q, int r) {
    int i, j, k;
    for (i = 0; i < p; ++i) {
        for (k = 0; k < r; ++k) {
            C[i * r + k] = 0.0;
        }
    }


    for (i = 0; i < p; ++i) {
        for (j = 0; j < q; ++j) {

            double a_ij = A[i * q + j];
            k = 0;
            for (; k + 3 < r; k += 4) {
                C[i * r + k]     += a_ij * B[j * r + k];
                C[i * r + k + 1] += a_ij * B[j * r + k + 1];
                C[i * r + k + 2] += a_ij * B[j * r + k + 2];
                C[i * r + k + 3] += a_ij * B[j * r + k + 3];
            }

            for (; k < r; ++k) {
                C[i * r + k] += a_ij * B[j * r + k];
            }
        }
    }
}

void matrix_sum(double *A, double *B, double *C, int n,
                int m) // C = A + B
{
    int i;
    for (i = 0; i < n * m; i++) {
        C[i] = A[i] + B[i];
    }
}

void matrix_subtr(double *A, double *B, int n, int m) // A -= B
{
    int i;
    for (i = 0; i < n * m; i++) {
        A[i] -= B[i];
    }
}

double get_norm(double *matrix, int m)
{
    int i, j;
    double max = 0, sum;
    for (j = 0; j < m; j++) {
        sum = 0;
        for (i = 0; i < m; i++) {
            sum += std::abs(matrix[j + i * m]);
        }
        // printf("sum is %8.3e\n", sum);
        max = std::max(max, sum);
    }
    return max;
}

int get_inverse_matrix(double *A, double *B, int m)
{
    int i, j, k, max_ind;
    double max, buf;
    unit_matrix(B, m);
    for (j = 0; j < m; j++) {
        // print_matrix(A, m, m);
        // printf("\n");
        max = std::abs(A[j + j * m]);
        max_ind = j;
        for (i = j + 1; i < m; i++) {
            if (std::abs(A[j + i * m]) > max) {
                max = std::abs(A[j + i * m]);
                max_ind = i;
            }
        }
        // printf("max: %lf, max_ind: %d\n", max, max_ind);

        if (max < EPS) {
            return -1;
        }
        // print_matrix(A, m, m);
        // printf("^^^A0\n");
        // print_matrix(B, m, m);
        // printf("^^^B0\n");
        // printf("max_ind is %d\n", max_ind);

        // Меняем местами строки, если нашли нужную
        if (max_ind != j) {
            for (i = j; i < m; i++) {
                buf = A[i + j * m];
                A[i + j * m] = A[i + max_ind * m];
                A[i + max_ind * m] = buf;
            }

            for (i = 0; i < m; i++) {
                buf = B[i + j * m];
                B[i + j * m] = B[i + max_ind * m];
                B[i + max_ind * m] = buf;
            }
        }
        
        max = A[j + j * m];
        A[j + j * m] = 1;
        for (i = j + 1; i < m; i++) {
            A[i + j * m] /= max;
        }

        for (i = 0; i < m; i++) {
            B[i + j * m] /= max;
        }

        // print_matrix(A, m, m);
        // printf("^^^A\n");
        // print_matrix(B, m, m);
        // printf("^^^B\n");
        for (i = 0; i < m; i++) {
            if (i != j) {
                max = A[j + i * m];
                A[j + i * m] = 0;
                for (k = j + 1; k < m; k++) {
                    A[k + i * m] -= max * A[k + j * m];
                }
                for (k = 0; k < m; k++) {
                    B[k + i * m] -= max * B[k + j * m];
                }
                // if (std::abs(max) > EPS) {
                // }
            }
        }
        // print_matrix(A, m, m);
        // printf("^^^A2\n\n");
        // print_matrix(B, m, m);
        // printf("^^^B2\n\n");

    }
    return 0;
}

// double mult_sub_norm_p(double *a, double *b, double *pc, double *norm, int n, int m, int p, int pi)
// {
//     int i, j, s, r, t, q;
//     int k = n / m;
//     int l = n - k * m; // n = k * m + l
//     int bl = (l != 0 ? k + 1 : k); // Общее количество блоков
//     int v, h, ah;
//     double max_norm = 0;

//     for (i = 0; i < n; i++) norm[i] = 0;

//     // Проходим по всем блокам в матрице C
//     for (i = 0; i < bl; i++) {
//         for (j = pi; j < bl; j += p) {
//             // Определяем размер текущего блока C[v x h]
//             v = (i < k ? m : l); // вертикальный размер блока
//             h = (j < k ? m : l); // горизонтальный размер блока

//             // Указатель на начало текущего блока C
//             // double *pc = c + (i * m) * n + j * m;

//             // Инициализируем блок C нулями
//             for (r = 0; r < v; r++) {
//                 for (t = 0; t < h; t++) {
//                     pc[r * h + t] = 0.0;
//                 }
//             }

//             // Перемножаем соответствующие блоки A и B и добавляем к C
//             for (s = 0; s < bl; s++) {
//                 // Определяем размер внутреннего блока
//                 ah = (s < k ? m : l);

//                 // Указатели на текущие блоки A и B
//                 double *pa = a + (i * m) * n + s * m; // блок A[i][s]
//                 double *pb = b + (s * m) * n + j * m; // блок B[s][j]

//                 // Основные циклы с разверткой для блоков 3x3
//                 int r_end = (v / 3) * 3;
//                 int t_end = (h / 3) * 3;

//                 // Обработка блоков размером 3x3
//                 for (r = 0; r < r_end; r += 3) {
//                     for (t = 0; t < t_end; t += 3) {
//                         double s00 = 0.0, s01 = 0.0, s02 = 0.0;
//                         double s10 = 0.0, s11 = 0.0, s12 = 0.0;
//                         double s20 = 0.0, s21 = 0.0, s22 = 0.0;

//                         for (q = 0; q < ah; q++) {
//                             double a0q = pa[(r + 0) * n + q];
//                             double a1q = pa[(r + 1) * n + q];
//                             double a2q = pa[(r + 2) * n + q];
//                             double bq0 = pb[q * n + (t + 0)];
//                             double bq1 = pb[q * n + (t + 1)];
//                             double bq2 = pb[q * n + (t + 2)];

//                             s00 += a0q * bq0;
//                             s01 += a0q * bq1;
//                             s02 += a0q * bq2;

//                             s10 += a1q * bq0;
//                             s11 += a1q * bq1;
//                             s12 += a1q * bq2;

//                             s20 += a2q * bq0;
//                             s21 += a2q * bq1;
//                             s22 += a2q * bq2;
//                         }

//                         pc[(r + 0) * h + (t + 0)] += s00;
//                         pc[(r + 0) * h + (t + 1)] += s01;
//                         pc[(r + 0) * h + (t + 2)] += s02;

//                         pc[(r + 1) * h + (t + 0)] += s10;
//                         pc[(r + 1) * h + (t + 1)] += s11;
//                         pc[(r + 1) * h + (t + 2)] += s12;

//                         pc[(r + 2) * h + (t + 0)] += s20;
//                         pc[(r + 2) * h + (t + 1)] += s21;
//                         pc[(r + 2) * h + (t + 2)] += s22;

//                         // if (r == t && i == j) {
//                         //     max_norm = MAX(std::abs(1. - s00), max_norm);
//                         //     max_norm = MAX(std::abs(1. - s11), max_norm);
//                         //     max_norm = MAX(std::abs(1. - s22), max_norm);
//                         // } else {
//                         //     max_norm = MAX(std::abs(s00), max_norm);
//                         //     max_norm = MAX(std::abs(s11), max_norm);
//                         //     max_norm = MAX(std::abs(s22), max_norm);
//                         // }

//                         // max_norm = MAX(std::abs(s01), max_norm);
//                         // max_norm = MAX(std::abs(s02), max_norm);
//                         // max_norm = MAX(std::abs(s10), max_norm);
//                         // max_norm = MAX(std::abs(s12), max_norm);
//                         // max_norm = MAX(std::abs(s20), max_norm);
//                         // max_norm = MAX(std::abs(s21), max_norm);

//                     }

//                     // Обработка оставшихся столбцов в блоке
//                     for (t = t_end; t < h; t++) {
//                         double s0 = 0.0, s1 = 0.0, s2 = 0.0;

//                         for (q = 0; q < ah; q++) {
//                             double a0q = pa[(r + 0) * n + q];
//                             double a1q = pa[(r + 1) * n + q];
//                             double a2q = pa[(r + 2) * n + q];
//                             double bqt = pb[q * n + t];

//                             s0 += a0q * bqt;
//                             s1 += a1q * bqt;
//                             s2 += a2q * bqt;
//                         }

//                         pc[(r + 0) * h + t] += s0;
//                         pc[(r + 1) * h + t] += s1;
//                         pc[(r + 2) * h + t] += s2;
//                         // max_norm = MAX(std::abs(s0), max_norm);
//                         // max_norm = MAX(std::abs(s1), max_norm);
//                         // max_norm = MAX(std::abs(s2), max_norm);
                        
//                     }
//                 }

//                 // Обработка оставшихся строк в блоке
//                 for (r = r_end; r < v; r++) {
//                     for (t = 0; t < t_end; t += 3) {
//                         double s0 = 0.0, s1 = 0.0, s2 = 0.0;

//                         for (q = 0; q < ah; q++) {
//                             double a0q = pa[r * n + q];
//                             double bq0 = pb[q * n + (t + 0)];
//                             double bq1 = pb[q * n + (t + 1)];
//                             double bq2 = pb[q * n + (t + 2)];

//                             s0 += a0q * bq0;
//                             s1 += a0q * bq1;
//                             s2 += a0q * bq2;
//                         }

//                         pc[r * h + (t + 0)] += s0;
//                         pc[r * h + (t + 1)] += s1;
//                         pc[r * h + (t + 2)] += s2;

//                         // max_norm = MAX(std::abs(s0), max_norm);
//                         // max_norm = MAX(std::abs(s1), max_norm);
//                         // max_norm = MAX(std::abs(s2), max_norm);

//                     }

//                     // Обработка оставшихся элементов
//                     for (t = t_end; t < h; t++) {
//                         double sum = 0.0;

//                         for (q = 0; q < ah; q++) {
//                             sum += pa[r * n + q] * pb[q * n + t];
//                         }

//                         pc[r * h + t] += sum;
//                         // if (r == t && i == j) max_norm = MAX(std::abs(1 - sum), max_norm);
//                         // else max_norm = MAX(std::abs(sum), max_norm);
//                     }
//                 }
//             }
//             for (r = 0; r < v; r++) {
//                 for (t = 0; t < h; t++) {
//                     if (i == j && r == t) {
//                         norm[t + m * j] += std::abs(1 - pc[r * h + t]);
//                     } else {
//                         norm[t + m * j] += std::abs(pc[r * h + t]);
//                     }
//                 }
//             }
//         }
//     }

//     // synchronize(p, &max_norm, 1, reduce::abs_max);
//     synchronize(p);
//     for (i = 0; i < n; i++) max_norm = MAX(max_norm, norm[i]);
//     synchronize(p);


//     return max_norm;
// }

void print_matrix_l_x_n(double *matrix, int l, int n)
{
    for (int i = 0; i < l; i++) {
        for (int j = 0; j < n; j++) {
            printf(" %10.3e", matrix[i * n + j]);
        }
        printf("\n");
    }
}

// int find_diff(double *matrix, double *inversed_matrix, double* block, double *norm, FILE *file, int n, int m, int s, double &r1, double &r2, int p, int pi)
// {
//     if (n < 11000) {
//         if (s == 0 && pi == 0) {
//             if (read_matrix_from_file(matrix, n, file) != 0) {
//                 return 2;
//             }
//         } else if (s != 0) {
//             if (fill_matrix(matrix, n, s) != 0) {
//                 return 2;
//             }
//         }
//         synchronize(p);

//         r1 = mult_sub_norm_p(matrix, inversed_matrix, block, norm, n, m, p, pi);
//         r2 = mult_sub_norm_p(inversed_matrix, matrix, block, norm, n, m, p, pi);
//         // printf("F pi: %d r1 = %8.3e, r2 = %8.3e\n", pi, r1, r2);

//     } else {
//         r1 = 0;
//         r2 = 0;
//     }
//     return 0;
// }


void zero_matrix_mpi(double *matrix, int n, int m, int p, int pi)
{
    int cols = get_loc_cols(n, m, p, pi);
    // std::cout << "count: " << cols * n << std::endl;
    memset(matrix, 0, cols * n * sizeof(double));
}

void unit_matrix_mpi(double *matrix, int n, int m, int p, int pi)
{
    int cols = get_loc_cols(n, m, p, pi);
    // std::cout << "cols: " << cols << std::endl;
    zero_matrix_mpi(matrix, n, m, p, pi);

    for (int i = pi * m; i < n; i += p * m) {
        int j_loc = g2l(n, m, p, pi, i);
        // std::cout << "i, j_loc: " << i << ", " << j_loc << std::endl;
        for (int j = 0; j < m && j_loc + j < cols; j++) {
            // std::cout << "pos: " << (i + j) * cols + (j_loc + j) << std::endl;
            matrix[(i + j) * cols + (j_loc + j)] = 1;
        }
    }
}


int fill_matrix_p(double *matrix, int n, int s, int m, int p, int pi)
{
    int i, j, k;
    switch (s) {
    case 1:
        for (i = 0; i < n; i++)
            for (j = pi * m; j < n; j+= p * m)
                for (k = j; k < j + m && k < n; k++)
                    matrix[i * n + k] = n - MAX(i, k);
        break;
    case 2:
        for (i = 0; i < n; i++)
            for (j = pi * m; j < n; j+= p * m)
                for (k = j; k < j + m && k < n; k++)
                    matrix[i * n + k] = MAX(i, k) + 1;
        break;
    case 3:
        for (i = 0; i < n; i++)
            for (j = pi * m; j < n; j+= p * m)
                for (k = j; k < j + m && k < n; k++)
                    matrix[i * n + k] = (i > k ? i - k : k - i);
        break;
    case 4:
        for (i = 0; i < n; i++)
            for (j = pi * m; j < n; j+= p * m)
                for (k = j; k < j + m && k < n; k++)
                    matrix[i * n + k] = 1. / (i + k + 1);
        break;
    default:
        fprintf(stderr, "[-] Unknown formula %d\n", s);
        return -1;
    }
    return 0;
}


double get_norm_p(double *matrix, int n, int m, int p, int pi)
{
    int i, j, k;
    double max = 0, sum;
    for (j = pi * m; j < n; j += p * m) {
        for (k = j; k < j + m && k < n; k++) {
            sum = 0;
            for (i = 0; i < n; i++) {
                sum += std::abs(matrix[k + i * n]);
            }
            max = std::max(max, sum);
            // printf("sum is %8.3e\n", sum);
        }
    }
    return max;
}



void rows_permutation_p(double *A, double *block1, double *block2, int n, int cols,
                      int m, int k, int l, int i1, int i2,
                      int begin, int p, int pi)
{
    int i;
    if (begin % p <= pi) begin += pi - begin % p;
    else begin += p + pi - begin % p;



    for (i = begin; i < k + 1; i += p) {
        get_block(A, block1, n, cols, m, k, l, i1, i, p, pi);
        get_block(A, block2, n, cols, m, k, l, i2, i, p, pi);
        put_block(A, block2, n, cols, m, k, l, i1, i, p, pi);
        put_block(A, block1, n, cols, m, k, l, i2, i, p, pi);
    }
}

double get_time() {
    struct timeval buf;
    gettimeofday(&buf, 0);
    return buf.tv_sec + buf.tv_usec / 1e6;
}


// void *thread_func(void *args)
// {
// 	Args *a = (Args*) args;
	
// 	int p = a->p, pi = a->pi, n = a->n, m = a->m, s = a->s, k = a->k, l = a->l, r = a->r;
//     int diag, i, j, bl = (l == 0) ? k : k + 1, min_norm_ind, row, x, y, z, begin;
//     double *matrix = a->matrix, *inversed_matrix = a->inversed_matrix, *block = a->block, *norm_arr = a->norm;
//     // double r1 = a->r1, r2 = a->r2;
//     double norm = 0, min_norm, buf_array[2];
//     FILE *file = a->file;
//     pthread_t tid = a->tid;
//     cpu_set_t cpu;
// 	int n_cpus = get_nprocs(); // число процессоров
// 	int cpu_id = n_cpus - 1 - (pi % n_cpus);

// 	CPU_ZERO(&cpu); // вместо конструктора
//     CPU_SET(cpu_id, &cpu);

//     pthread_setaffinity_np(tid, sizeof(cpu), &cpu);

//     double *block_A = new double[m*m];
//     double *block_B = new double[m*m];
//     double *block_C = new double[m*m];

//     if (!block_A || !block_B || !block_C) {
//         if (block_A) delete[] block_A;
//         if (block_B) delete[] block_B;
//         if (block_C) delete[] block_C;
//         a->error_flag = 1;
//         a->error_type = io_status::bad_allocation;
//     }
// 	synchronize(p, & a->error_flag, 1);
// 	if (a -> error_flag > 0) {
//         delete[] block_A;
//         delete[] block_B;
//         delete[] block_C;
// 		return nullptr;
// 	}

//     a->t1 = get_time();

//     zero_matrix_p(matrix, n, m, p, pi);
//     unit_matrix_p(inversed_matrix, n, m, p, pi);
//     synchronize(p);


//     if (s == 0) {
//         if (pi == 0) {
//             if (read_matrix_from_file(matrix, n, file) != 0) {
//                 a->error_flag = 1;
//                 a->error_type = io_status::error_read;
//             }
//         }
//     } else {
//         if (fill_matrix_p(matrix, n, s, m, p, pi) != 0) {
//             a->error_flag = 1;
//             a->error_type = io_status::unknown_formula;
//     
//     }

//     if (pi == 0) {
//         printf("[+] Given matrix:\n");
//         print_matrix(matrix, n, 4);
//     }
// 	// printf("Yep!\n");
	
// 	synchronize(a->p, & a->error_flag, 1);
// 	if (a -> error_flag > 0) {
//         delete[] block_A;
//         delete[] block_B;
//         delete[] block_C;
// 		return nullptr;
// 	}

//     norm = get_norm_p(matrix, n, m, p, pi);
// 	synchronize(a->p, &norm, 1);
// 	if (pi == 0) EPS *= norm;
//     // printf("norm = %8.3e, eps = %8.3e\n", norm, EPS);
// 	synchronize(a->p);

//     // Пробегаюсь по диагональным элементам
//     for (diag = 0; diag < bl; diag++) {
//         // print_matrix(matrix, n, n);
//         // printf("\n");
//         min_norm = -1;
//         min_norm_ind = -1;
//         // Ищу в столбце матрицу у кооторой обратная имеет наименьшую норму
//         if (diag < k) {
//             for (row = diag + pi; row < k; row += p) {
//                 get_block(matrix, block_A, n, m, k, l, diag, row);
//                 if (get_inverse_matrix(block_A, block_B, m) == 0) {
//                     // print_matrix(block_B, m, print_size);
//                     norm = get_norm(block_B, m);
//                     // printf("norm is %8.3e\n", norm);
//                     if (min_norm < 0 || norm < min_norm) {
//                         min_norm = norm;
//                         min_norm_ind = row;
//                     }
//                 }
//             }
//         } else {
//             if (pi == 0) {
//                 get_block(matrix, block_A, n, m, k, l, k, k);
//                 if (get_inverse_matrix(block_A, block_B, l) == 0) {
//                     norm = get_norm(block_B, l);
//                     min_norm = norm;
//                     min_norm_ind = k;
//                 }
//             }
//         }
//         buf_array[0] = min_norm;
//         buf_array[1] = min_norm_ind * 1. > 1e+13 ? -1. : min_norm_ind;
//         // buf_array[1] = buf_array[1] * 1. > 1e+13 ? -1. : buf_array[1];

        
//         // printf("%d min_norm: %8.3e, ind: %.0e ---buf\n", pi, buf_array[0], buf_array[1]);
//         // printf("%d min_norm: %8.3e, ind: %lu\n", pi, min_norm, min_norm_ind);
//         synchronize(p, buf_array, 2, reduce::abs_min_first);
//         // if (pi == 0) printf("%d min_norm: %8.3e, ind: %.0e\n", pi, buf_array[0], buf_array[1]);
//         // printf("Yep!\n");
//         if (buf_array[1] < 0) {

//             fprintf(stderr, "[-] Matrix is irreversible, place: 923!\n");
//             a->error_flag = 1;
//             a->error_type = io_status::irreversible;
//             delete[] block_A;
//             delete[] block_B;
//             delete[] block_C;
//             return nullptr;
//         }
//         // min_norm = buf_array[0];
//         min_norm_ind = buf_array[1];
//         // printf("Step: %d, thread: %d\n", diag, pi);

//         // Переставляю строки так, чтобы на текущем диагональном элементе была
//         // матрица с наименьшей нормой её обратной
//         rows_permutation_p(
//             matrix, block_A, block_B, n, m, k, l, min_norm_ind, diag, diag, p, pi
//         );
//         rows_permutation_p(
//             inversed_matrix, block_A, block_B, n, m, k, l, min_norm_ind, diag, 0, p, pi
//         );
        
//         // Нахожу обратную в потоке pi
//         if (diag % p == pi) {
//             get_block(matrix, block_A, n, m, k, l, diag, diag);
//             if (diag != k) {
//                 if (get_inverse_matrix(block_A, block, m) != 0) {
//                     fprintf(stderr, "[-] Matrix is irreversible, place 949!\n");
//                     a->error_flag = 1;
//                     a->error_type = io_status::irreversible;
//                     delete[] block_A;
//                     delete[] block_B;
//                     delete[] block_C;
//                     return nullptr;
//                 }
//             } else {
//                 if (get_inverse_matrix(block_A, block, l) != 0) {
//                     fprintf(stderr, "[-] Matrix is irreversible!\n");
//                     a->error_flag = 1;
//                     a->error_type = io_status::irreversible;
//                     delete[] block_A;
//                     delete[] block_B;
//                     delete[] block_C;
//                     return nullptr;
//                 }
//             }
//         }
//         synchronize(p);

//         memcpy(block_B, block, m * m * sizeof(double));
        
//         // Вставляю единичную на место текущего диагонального элемента
//         // put_block(matrix, block_A, n, m, k, l, diag, diag);

//         // Начиная со следующего каждый элемент в строке домножаю на обратную к диагональному
//         // block_B - обратная к диагональному
//         begin = diag + 1;
//         if (begin % p <= pi) begin += pi - begin % p;
//         else begin += pi + p - begin % p;

//         for (i = begin; i < bl; i+=p) {
//             get_block(matrix, block_A, n, m, k, l, diag, i);
//             x = (diag == k ? l : m);
//             y = (i == k    ? l : m);
//             matrix_multiply(block_B, block_A, block_C, x, x, y);
//             put_block(matrix, block_C, n, m, k, l, diag, i);
//         }

//         // Каждый элемент той же строки матрицы В домножаю на эту же обратную
//         for (i = pi; i < bl; i+=p) {
//             get_block(inversed_matrix, block_A, n, m, k, l, diag, i);
//             // printf("------matrix------\n");
//             // printf("Diag = %d, i = %d\n", diag, i);
//             // print_matrix(block_A, (i < k ? m : l), 4);
//             x = (diag == k ? l : m);
//             y = (i == k    ? l : m);
//             // printf("x = %d, y = %d\n", x, y);
//             // printf("------------------\n");
//             matrix_multiply(block_B, block_A, block_C, x, x, y);
//             put_block(inversed_matrix, block_C, n, m, k, l, diag, i);
//         }

//         synchronize(p);

//         // Для каждой строки, кроме той, на которой есть текущий диагональный
//         for (i = 0; i < bl; i++) {
//             if (i != diag) {
//                 // Запоминаем блок диагонального столбца и текущей строки

//                 get_block(matrix, block_A, n, m, k, l, i, diag);
//                 // Для каждого следующего элемента текущей строки, вычитаем из него А х С
//                 for (j = begin; j < bl; j+=p) { // j - номер столбца

//                     get_block(matrix, block_C, n, m, k, l, diag, j);

//                     x = (i == k    ? l : m);
//                     y = (diag == k ? l : m);
//                     z = (j == k    ? l : m);
//                     matrix_multiply(block_A, block_C, block_B, x, y, z);

//                     get_block(matrix, block_C, n, m, k, l, i, j);
//                     matrix_subtr(block_C, block_B, x, z);
//                     put_block(matrix, block_C, n, m, k, l, i, j);
//                 }

//                 // Аналогично для матрицы В
//                 for (j = pi; j < bl; j+=p) {
//                     get_block(inversed_matrix, block_C, n, m, k, l, diag, j);
//                     x = (i == k    ? l : m);
//                     y = (diag == k ? l : m);
//                     z = (j == k    ? l : m);
//                     matrix_multiply(block_A, block_C, block_B, x, y, z);

//                     get_block(inversed_matrix, block_C, n, m, k, l, i, j);
//                     matrix_subtr(block_C, block_B, x, z);
//                     put_block(inversed_matrix, block_C, n, m, k, l, i, j);
//                 }
//             }
//         }
//         synchronize(p);
//     }

// 	synchronize(p, & a-> error_flag, 1);
// 	if (a -> error_flag > 0) {
// 		delete[] block_A;
//         delete[] block_B;
//         delete[] block_C;
// 		return nullptr;
// 	}

//     a->t1 -= get_time();

//     if (pi == 0) {
//         printf("\n[+] Inversed matrix:\n");
//         print_matrix(inversed_matrix, n, r);
//     }


//     a->t2 = get_time();

//     if (find_diff(matrix, inversed_matrix, block_A, norm_arr, file, n, m, s, a->r1, a->r2, p, pi) != 0) {
//         a->error_type = io_status::error_read;
//         a->error_flag = 1;
//     }
// 	synchronize(p, & a-> error_flag, 1);
//     if (a -> error_flag > 0) {
// 		delete[] block_A;
//         delete[] block_B;
//         delete[] block_C;
// 		return nullptr;
// 	}

//     a->t2 -= get_time();


// 	a -> error_type = io_status::success;

//     delete[] block_A;
//     delete[] block_B;
//     delete[] block_C;
// 	return nullptr;
// }


// local to global
int l2g (
	int /*n*/,
	int m,
	int p,
	int pi,
	int i_loc
) { 
	// Номер блочного локального столбца
	int i_loc_m = i_loc / m;
	// Номер блочного глобального столбца
	int i_glob_m = i_loc_m * p + pi;
	return i_glob_m * m + i_loc % m;
}

// global to local
int g2l (
	int /*n*/,
	int m,
	int p,
	int /*pi*/,
	int i_glob
) {
	// номер блочного глобального столбца
	int i_glob_m = i_glob / m;
	// номер блочного локального столбца
	int i_loc_m = i_glob_m / p;
	return i_loc_m * m + i_glob % m;
}

// максимальное число столбцов локальной матрицы
int get_max_cols(
	int n,
	int m,
	int p
) {
	int b = (n + m - 1) / m;
	return (b + p - 1) / p * m;
}

// число блочных столбцов на процесс
int get_bl_cols(
	int n,
	int m,
	int p,
	int pi
) {
	int b = (n + m - 1) / m;
	return b % p <= pi ? b / p : b / p + 1;
}

// число столбцов в локальной матрице
int get_loc_cols(
    int n,
    int m,
    int p,
    int pi
) {
    // b - число блочных столбцов в глобальной матрице
    int b = get_bl_cols(n, m, p, pi);
    int b_loc = b * m - m;
    int l = n % m;
    int k = (n + m - 1) / m - 1;

    if (k % p == pi && l > 0) {
        b_loc += l;
    } else {
        b_loc += m;
    }
    
    return b_loc;
}

// в каком процессе лежит столбец?
int get_pi(
	int /*n*/,
	int m,
	int p,
	int i_glob
) {
	int i_glob_m = i_glob / m;
	return i_glob_m % p;
}

int function(
    int s,
    int n,
    int i,
    int j
) {
    // функция для инициализации матрицы
    // a_ij = f(s, n, i, j)
    if (s == 1) return n - MAX(i, j);
    if (s == 2) return MAX(i, j) + 1;
    if (s == 3) return (i > j ? i - j : j - i);
    if (s == 4) return 1. / (i + j + 1);
    if (s == 5) return i * n + j;
    return -1;
}

//a_ij = f(s, n, i, j)
void init_matrix(
	double *a,
	int n,
	int m,
	int p,
	int pi,
	int s
) {
	int i_loc, j_loc, i_glob, j_glob, cols;
	// сколько строк в процессе

    cols = get_loc_cols(n, m, p, pi);
    for (j_loc = 0; j_loc < cols; j_loc++) {
        j_glob = l2g(n, m, p, pi, j_loc);
        for (i_loc = 0; i_loc < n; i_loc++) {
            i_glob = i_loc;
            a[i_loc * cols + j_loc] = function(s, n, i_glob, j_glob);
        }
    }
}

int read_matrix(
	double *a,
	int n,
	int m,
	int p,
	int pi,
	const char *name,
	double *buf, // буффер - блочная строка n * m
	MPI_Comm com
) {
	int main_pi = 0; // кто читает файл
	FILE *fp = nullptr;
	int err = 0;
	if (pi == main_pi) {
		fp = fopen("r", name);
		if (fp == nullptr) err = 1;
	}
	MPI_Bcast(&err, 1, MPI_INT, main_pi, com);
	if (err) return err; // во всех процессах
	
	memset(buf, 0, n * m * sizeof(double));
	
	// число блочных строк
	int b, max_b = (n + m - 1) / m;
	for (b = 0; b < max_b; b++) {
		// владелец строки
		// int owner = b % p;
		int rows = b * m + m <= n ? m : n - b * m;
		// rows = min (m, n - b * m)
		
		// лок номер столбца
		// int b_loc = b / p;
		if (pi == main_pi) {
			err += read_array(fp, buf, n * rows);
		}

        MPI_Bcast(buf, n * rows, MPI_DOUBLE, main_pi, com);

        int cols = get_bl_cols(n, m, p, pi);
        int cols_loc = get_loc_cols(n, m, p, pi);
        int pos = 0;
        for (int i = 0; i < rows; i++) {
            for (int j = pi; j < cols; j += p) {
                for (int k = j * m; k < (j + 1) * m && k < n; k++) {
                    a[(i + b * m) * cols_loc + pos++] = buf[i * n + k];
                }
            }
        }
	}
	
	if (pi == main_pi) {
		fclose(fp);
		fp = nullptr;
	}
	MPI_Bcast(&err, 1, MPI_INT, main_pi, com);
	if (err) return err;
	return 0;
}

int read_array(FILE *fp, double *a, int len)
{
	for (int i = 0; i < len; i++) {
		if (fscanf(fp, "%lf", a + i) != 1) return -2;
	}
	return 0;
}

double get_norm_pi(double *matrix, int n, int cols)
{
    int i, j;
    double max = 0, sum;
    for (j = 0; j < cols; j++) {
        sum = 0;
        for (i = 0; i < n; i++) {
            sum += std::abs(matrix[j + i * cols]);
        }
        max = std::max(max, sum);
    }
    return max;
}


void get_column(
    double *matrix,
    double *buffer,
    int n,
    int m,
    int p,
    int pi,
    int i
)
{
    // int bl_cols = get_bl_cols(n, m, p, pi);
    int cols = get_loc_cols(n, m, p, pi);
    for (int j = 0; j < m; j++) {
        for (int k = 0; k < n; k++) {
            buffer[k * m + j] = matrix[i * m + k * cols + j];
        }
    }
}


// печать матрицы
void print_matrix_mpi(
	double * a,
	int n,
	int m,
	int p,
	int pi,
	double *buf, // n * m блочная строка
	int max_print,
	MPI_Comm com
) {
	int main_pi = 0; // только 0 в большинстве систем
	int printed_rows = 0;
    MPI_Status st;
    int cols = get_loc_cols(n, m, p, pi);
    int k = n / m;
    
    // // printf("[+] Process %d, cols: %d\n", pi, cols);
    // MPI_Barrier(com);
    // printf("----------------------------------------------------\n");
    // for (int i = 0; i < p; i++) {
    //     MPI_Barrier(com);
    //     if (i == pi) {
    //         printf("\n[+] Process %d:\n", pi);
    //         for (int j = 0; j < n; j++) {
    //             for (int k = 0; k < cols; k++) {
    //                 printf(" %10.3e", a[j * cols + k]);
    //             }
    //             printf("\n");
    //         }
    //     }
    // }
    // printf("----------------------------------------------------\n");
    // MPI_Barrier(com);


    
    // Отправка всех строк толщиной m
    for (int i = 0; i < k; i++) {
        MPI_Barrier(com);
        if (pi == main_pi) {
            memcpy(buf, a + i * m * cols, cols * m * sizeof(double));

            int p_shift = cols * m;
            for (int pk = 1; pk < p; pk++) {
                int pk_cols = get_loc_cols(n, m, p, pk);
                // printf("[+] pk, pk_cols: %d, %d\n", pk, pk_cols);
                MPI_Recv(buf + p_shift, pk_cols * m, MPI_DOUBLE, pk, 0, com, &st);
                p_shift += pk_cols * m;
            }

            // // печать всей блочной строки
            // printf("-------------------------------------------------\n");
            // printf("\n[+] Step %d:\n", i);
            // for (int j = 0; j < m; j++) {
            //     for (int k = 0; k < n; k++) {
            //         printf(" %10.3e", buf[j * n + k]);
            //     }
            //     printf("\n");
            // }
            // printf("-------------------------------------------------\n");


            print_array(buf, n, m, m, max_print, printed_rows, p, m);
        }
        else {
            MPI_Send(a + i * m * cols, cols * m, MPI_DOUBLE, main_pi, 0, com);
        }
        MPI_Bcast(&printed_rows, 1, MPI_INT, main_pi, com);
        if (printed_rows >= max_print) {
            return;
        }
    }

    // printf("[+] printed_rows: %d\n", printed_rows);

    int l = n % m;
    if (l == 0) {
        return;
    }

    if (pi == main_pi) {
        memcpy(buf, a + k * m * cols, cols * l * sizeof(double));

        int p_shift = cols * l;
        for (int pk = 1; pk < p; pk++) {
            int pk_cols = get_loc_cols(n, m, p, pk);
            MPI_Recv(buf + p_shift, pk_cols * l, MPI_DOUBLE, pk, 0, com, &st);
            p_shift += pk_cols * l;
        }

        // // печать всей блочной строки
        // printf("\n-------------------------------------------------\n");
        // printf("[+] Step %d:\n", k);
        // for (int j = 0; j < l; j++) {
        //     for (int k = 0; k < n; k++) {
        //         printf(" %10.3e", buf[j * n + k]);
        //     }
        //     printf("\n");
        // }
        // printf("-------------------------------------------------\n");
        print_array(buf, n, m, l, max_print, printed_rows, p, l);
    }
    else {
        MPI_Send(a + k * m * cols, cols * l, MPI_DOUBLE, main_pi, 0, com);
    }
}

// печать прямоугольной матрицы с адресом a, длиной строки n, числом строк rows
// число напечатаных строк <= max_print
// возвращает число напечатаных строк
void print_array(
	double *a,
	int n,
    int m,
    int ,
	int max_print,
    int & printed_rows,
    int p,
    int rows
) {
    int bl_cols = get_bl_cols(n, m, p, 0);
    max_print = std::min(max_print, n);

    for (int row = 0; row < rows && printed_rows < max_print; row++) {
        int printed = 0;
        for (int bl_col = 0; bl_col < bl_cols; bl_col++) {
            int skip = 0;
            for (int pi = 0; pi < p; pi++) {
                int cols = get_loc_cols(n, m, p, pi);
                // printf(" |pi: %d, cols: %d|", pi, cols);
                for (int i = 0; i < std::min(m, cols - bl_col * m); i++) {
                    printf(" %10.3e", a[skip + bl_col * m + i + row * cols]);
                    printed++;
                    if (printed >= max_print) {
                        printed_rows++;
                        printf("\n");
                        break;
                    }
                }
                skip += cols * rows;
                if (printed >= max_print) break;
            }
            if (printed >= max_print) break;
        }
        if (printed < max_print) printf("\n");
    }
    // printf("\n");
}

int mpi_calculate(
    double * matrix,            // n x (m * bl_cols)
    double * inversed_matrix,   // n x (m * bl_cols)
    double * buffer,
    int n,
    int m,
    int p,
    int pi,
    MPI_Comm com
)
{
    // int bl_cols = get_bl_cols(n, m, p, pi);
    int cols = get_loc_cols(n, m, p, pi);
    int err = 0, sum = 0;

    
    int bl = (n + m - 1) / m;
    int k = n / m;
    int l = n % m;
    int min_norm_ind;
    double min_norm;

    double *block_A = new double[m * m];
    double *block_B = new double[m * m];
    double *block_C = new double[m * m];
    double *buf_array = new double[2 * p];
    double min_elem[2];

    if (!block_A || !block_B || !block_C) {
        err = 1;
    }

    MPI_Allreduce(&err, &sum, 1, MPI_INT, MPI_SUM, com);
    if (sum) {
        if (pi == 0)
            fprintf(stderr, "[-] Error in allocation: %d\n", __LINE__);
        if (block_A)
            delete[] block_A;
        if (block_B)
            delete[] block_B;
        if (block_C)
            delete[] block_C;
        if (buf_array)
            delete[] buf_array;
        return -1;
    }
    
    double norm = get_norm_pi(matrix, n, cols);
    MPI_Allreduce(&norm, buffer, 1, MPI_DOUBLE, MPI_MAX, com);
    if (pi == 0) EPS *= norm;
    MPI_Barrier(com);
    
    // Пробегаюсь по диагональным элементам
    for (int diag = 0; diag < bl; diag++) {
        min_norm = -1.;
        min_norm_ind = -1;
        get_column(matrix, buffer, n, m, p, pi, diag);
        MPI_Bcast(buffer, n * m, MPI_DOUBLE, diag % p, com);
        // Ищу в столбце матрицу у кооторой обратная имеет наименьшую норму
        if (diag < k) {
            for (int row = diag + pi; row < k; row += p) {
                get_block(matrix, block_A, n, cols, m, k, l, diag, row, p, pi);
                if (get_inverse_matrix(block_A, block_B, m) == 0) {
                    // print_matrix(block_B, m, print_size);
                    norm = get_norm(block_B, m);
                    // printf("norm is %8.3e\n", norm);
                    if (min_norm < 0 || norm < min_norm) {
                        min_norm = norm;
                        min_norm_ind = row;
                    }
                }
            }
        } else {
            if (pi == k % p) {
                get_block(matrix, block_A, n, cols, m, k, l, k, k, p, pi);
                if (get_inverse_matrix(block_A, block_B, l) == 0) {
                    norm = get_norm(block_B, l);
                    min_norm = norm;
                    min_norm_ind = k;
                }
            }
        }
        min_elem[0] = min_norm;
        min_elem[1] = min_norm_ind;
        // buf_array[1] = buf_array[1] * 1. > 1e+13 ? -1. : buf_array[1];

        
        // printf("%d min_norm: %8.3e, ind: %.0e ---buf\n", pi, buf_array[0], buf_array[1]);
        // printf("%d min_norm: %8.3e, ind: %lu\n", pi, min_norm, min_norm_ind);
        MPI_Allgather(min_elem, 2, MPI_DOUBLE, buf_array, 2, MPI_DOUBLE, com);

        min_norm = -1;
        min_norm_ind = -1;
        for (int i = 0; i < p; i++) {
            if (buf_array[i * 2] < min_norm || (min_norm < EPS && buf_array[i * 2] >= EPS)) {
                min_norm = buf_array[i * 2];
                min_norm_ind = round(buf_array[i * 2 + 1]);
            }
        }

        if (min_norm_ind < 0) {
            if (pi == 0)
                fprintf(stderr, "[-] Matrix is irreversible, place: %d!\n", __LINE__);
            delete[] block_A;
            delete[] block_B;
            delete[] block_C;
            delete[] buf_array;
            return -2;
        }

        // Переставляю строки так, чтобы на текущем диагональном элементе была
        // матрица с наименьшей нормой её обратной
        rows_permutation_p(
            matrix, block_A, block_B, n, cols, m, k, l, min_norm_ind, diag, diag, p, pi
        );
        rows_permutation_p(
            inversed_matrix, block_A, block_B, n, cols, m, k, l, min_norm_ind, diag, 0, p, pi
        );
        
        // Нахожу обратную
        get_block(buffer, block_A, n, m, m, k, l, diag, 0, p, pi);
        if (diag != k) {
            if (get_inverse_matrix(block_A, block_B, m) != 0) {
                if (pi == 0)
                    fprintf(stderr, "[-] Matrix is irreversible, place %d!\n", __LINE__);
                delete[] block_A;
                delete[] block_B;
                delete[] block_C;
                delete[] buf_array;
                return -2;
            }
        } else {
            if (get_inverse_matrix(block_A, block_B, l) != 0) {
                if (pi == 0)
                    fprintf(stderr, "[-] Matrix is irreversible, place %s:%d!\n", __FILE__, __LINE__);
                delete[] block_A;
                delete[] block_B;
                delete[] block_C;
                delete[] buf_array;
                return -2;
            }
        }

        MPI_Barrier(com);

        
        // Вставляю единичную на место текущего диагонального элемента
        // put_block(matrix, block_A, n, m, k, l, diag, diag);

        // Начиная со следующего каждый элемент в строке домножаю на обратную к диагональному
        // block_B - обратная к диагональному
        int begin = diag + 1;
        if (begin % p <= pi) begin += pi - begin % p;
        else begin += pi + p - begin % p;

        int x, y, z;
        for (int i = begin; i < bl; i+=p) {
            get_block(matrix, block_A, n, cols, m, k, l, diag, i, p, pi);
            x = (diag == k ? l : m);
            y = (i == k    ? l : m);
            matrix_multiply(block_B, block_A, block_C, x, x, y);
            put_block(matrix, block_C, n, cols, m, k, l, diag, i, p, pi);
        }

        // Каждый элемент той же строки матрицы В домножаю на эту же обратную
        for (int i = pi; i < bl; i+=p) {
            get_block(inversed_matrix, block_A, n, cols, m, k, l, diag, i, p, pi);
            // printf("------matrix------\n");
            // printf("Diag = %d, i = %d\n", diag, i);
            // print_matrix(block_A, (i < k ? m : l), 4);
            x = (diag == k ? l : m);
            y = (i == k    ? l : m);
            // printf("x = %d, y = %d\n", x, y);
            // printf("------------------\n");
            matrix_multiply(block_B, block_A, block_C, x, x, y);
            put_block(inversed_matrix, block_C, n, cols, m, k, l, diag, i, p, pi);
        }

        MPI_Barrier(com);

        // Для каждой строки, кроме той, на которой есть текущий диагональный
        for (int i = 0; i < bl; i++) {
            if (i != diag) {
                // Запоминаем блок диагонального столбца и текущей строки

                get_block(matrix, block_A, n, cols, m, k, l, i, diag, p, pi);
                // Для каждого следующего элемента текущей строки, вычитаем из него А х С
                for (int j = begin; j < bl; j+=p) { // j - номер столбца

                    get_block(matrix, block_C, n, cols, m, k, l, diag, j, p, pi);

                    x = (i == k    ? l : m);
                    y = (diag == k ? l : m);
                    z = (j == k    ? l : m);
                    matrix_multiply(block_A, block_C, block_B, x, y, z);

                    get_block(matrix, block_C, n, cols, m, k, l, i, j, p, pi);
                    matrix_subtr(block_C, block_B, x, z);
                    put_block(matrix, block_C, n, cols, m, k, l, i, j, p, pi);
                }

                // Аналогично для матрицы В
                for (int j = pi; j < bl; j+=p) {
                    get_block(inversed_matrix, block_C, n, cols, m, k, l, diag, j, p, pi);
                    x = (i == k    ? l : m);
                    y = (diag == k ? l : m);
                    z = (j == k    ? l : m);
                    matrix_multiply(block_A, block_C, block_B, x, y, z);

                    get_block(inversed_matrix, block_C, n, cols, m, k, l, i, j, p, pi);
                    matrix_subtr(block_C, block_B, x, z);
                    put_block(inversed_matrix, block_C, n, cols, m, k, l, i, j, p, pi);
                }
            }
        }
        MPI_Barrier(com);
    }

	MPI_Barrier(com);

    // a->t1 -= get_time();

    if (pi == 0) {
        printf("\n[+] Inversed matrix:\n");
        print_matrix_mpi(
            inversed_matrix, n, m, p, pi, buffer, 4, com
        );
        // print_matrix(inversed_matrix, n, r);
    }


    // a->t2 = get_time();

    // if (find_diff(matrix, inversed_matrix, block_A, norm_arr, file, n, m, s, a->r1, a->r2, p, pi) != 0) {
    //     err = 1;
    // }
    // MPI_Allreduce(&err, &err, 1, MPI_INT, MPI_SUM, com);
    // if (err > 0) {
	// 	delete[] block_A;
    //     delete[] block_B;
    //     delete[] block_C;
    //     delete[] buf_array;
	// 	return -3;
	// }

    // a->t2 -= get_time();



    delete[] block_A;
    delete[] block_B;
    delete[] block_C;
    delete[] buf_array;
	return 0;
}

