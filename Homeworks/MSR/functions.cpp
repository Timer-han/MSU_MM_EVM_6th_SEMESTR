#include <cmath>
#include <sys/resource.h>
#include <sys/sysinfo.h>
#include <sys/time.h>
#include <pthread.h>  

#include "functions.h"
#include "reduce.h"

double f_0(double /* x */, double /* y */)
{
	return 1;
}

double f_1(double x, double /* y */)
{
	return x;
}

double f_2(double /* x */, double y)
{
	return y;
}

double f_3(double x, double y)
{
	return x + y;
}

double f_4(double x, double y)
{
	return sqrt(x * x + y * y);
}

double f_5(double x, double y)
{
	return x * x + y * y;
}

double f_6(double x, double y)
{
	return exp(x * x - y * y);
}

double f_7(double x, double y)
{
	return 1. / (25 *(x * x + y * y) + 1);
}

void select_func(int k, double (*&f)(double, double)) {
	switch (k) {
		case 0:
			f = f_0;
			break;
		case 1:
			f = f_1;
			break;
		case 2:
			f = f_2;
			break;
		case 3:
			f = f_3;
			break;
		case 4:
			f = f_4;
			break;
		case 5:
			f = f_5;
			break;
		case 6:
			f = f_6;
			break;
		case 7:
			f = f_6;
			break;
    }
}

int minimal_residual_msr_matrix(
    int n,
    double *A,
    int *I,
    double *b,
    double *x,
    double *r,
    double *u,
    double *v,
    double eps,
    int max_it,
    int p,
    int pi)
{
    double prec, b_norm2, tau;
    double c1, c2;
    int it = 0;
    b_norm2 = scalar_product(n, b, b, p, pi); // (b, b)
    prec = b_norm2 * eps * eps;
    // r = Ax
    mult_msr_matrix_vector(n, A, I, x, r, p, pi);
    // 1 точка синхронизации
	// r = Ax - b, r-=b, r-=1*b
    mult_sub_vector(n, r, b, 1., p, pi);
	// r -= 1.*b, 1 точка синхронизации
    for (; it < max_it; ++it)
    {
		// Mr = r - решить систему
        apply_preconditioner_msr_matrix(
            n, A, I, v, r, p, pi);
        // 1 точка синхронизациии
		// u = Av, u = AM^(-1)r
        mult_msr_matrix_vector(n, A, I, v, u, p, pi); 
		// 1 точка синхронизации
		// c_1 = (u, r)
		// c_2 = (u, u)
        c1 = scalar_product(n, u, r, p, pi);
        c2 = scalar_product(n, u, u, p, pi);

        if (c1 < prec || c2 < prec)
            break;

        tau = c1 / c2;
        mult_sub_vector(n, x, v, tau, p, pi); 
		// 1 точка синхронизации
		// r -= tau * u
        mult_sub_vector(n, r, u, tau, p, pi);
		// 1 точка синхронизации
    }
    if (it > max_it)
    {
        return -1;
    }
    return it;
}

int minimal_residual_msr_matrix_full(
    int n,
    double *A,
    int *I,
    double *b,
    double *x, /*Начальное, а в конце будет ответ*/
    double *r,
    double *u,
    double *v,
    double eps,
    int max_it,
    int max_step,
    int p,
    int pi)
{
    int step, ret, its = 0;
    for (step = 0; step < max_step; ++step) {
        ret = minimal_residual_msr_matrix(
            n, A, I, b, x, r, u, v, eps, max_it, p, pi
        );
        if (ret >= 0) {
            its += ret;
            break;
        }
        its += max_it;
    }
    if (step >= max_step)
        return -1;
    return its;
}

void thread_rows(
    int n,
    int p,
    int pi,
    int &i1,
    int &i2)
{
    i1 = n * pi,
    i1 /= p;
    i2 = n * (pi + 1);
    i2 /= p;
}

double scalar_product(
    int n,
    double *x,
    double *y,
    int p,
    int pi)
{
    int i1;
    int i2;
    int i;
    double s = 0;
    thread_rows(n, p, pi, i1, i2);
    for (i = i1; i < i2; i++)
        s += x[i] * y[i];

    s = reduce_sum_det(p, pi, s);
    // ответ в каждом потоке
    return s;
}

void mult_sub_vector(
    int n,
    double *x,
    double *y,
    double t,
    int p,
    int pi)
{
    int i, i1, i2;
    thread_rows(n, p, pi, i1, i2);
    for (i = i1; i < i2; i++)
        x[i] -= t * y[i];
    synchronize(p);
}

void apply_preconditioner_msr_matrix(
    int n,
    double *A,
    int * /*I*/,
    double *v,
    double *r,
    int p,
    int pi)
{
    int i, i1, i2;
    thread_rows(n, p, pi, i1, i2);
    for (i = i1; i < i2; i++)
        v[i] = r[i] / A[i];
    synchronize(p);
}

void mult_msr_matrix_vector(
    int n,
    double *A,
    int *I,
    double *x,
    double *y,
    int p,
    int pi)
{
    int i, i1, i2;
    int l, J;
    double s;
    thread_rows(n, p, pi, i1, i2);
    for (i = i1; i < i2; i++)
    {
		// диагональный элемент
        s = A[i] * x[i]; // A_ii * x_i
		// число != 0 вне диагональных элементов
        l = I[i + 1] - I[i];
		// начало строки i
        J = I[i];
        for (int j = 0; j < l; j++)
            s += A[J + j] * x[I[J + j]];
        y[i] = s;
    }
    synchronize(p);
}

void ij2l(int nx, int /*ny*/, int i, int j, int &l)
{
    l = i + j * (nx + 1);
}

void l2ij(int nx, int /*ny*/, int &i, int &j, int l)
{
    j = l / (nx + 1);
    i = l - j * (nx + 1);
}

int get_len_msr(int nx, int ny)
{
    return (nx + 1) * (ny + 1) + 6 * (nx - 1) * (ny - 1) + 4 * (2 * (nx - 1) + 2 * (ny - 1)) + 2 * 3 + 2 * 2;
}

#define F(I, J)                \
    do                         \
    {                          \
        ij2l(nx, ny, I, J, k); \
        if (I_ij)              \
        {                      \
            I_ij[m] = k;       \
        }                      \
        m++;                   \
    } while (0)

int get_all_diag(int nx, int ny, int i, int j, int *I_ij)
{
    int k, m = 0;
    if (i < nx)
    {
        F(i + 1, j);
    }
    if (j > 0)
    {
        F(i, j - 1);
    }
    if (i > 0 && j > 0)
    {
        F(i - 1, j - 1);
    }
    if (i > 0)
    {
        F(i - 1, j);
    }
    if (j < ny)
    {
        F(i, j + 1);
    }
    if (i < nx && j < ny)
    {
        F(i + 1, j + 1);
    }
    return m;
}

int get_len_msr_all_diag(int nx, int ny)
{
    int i, j, m = 0;

    for (i = 0; i <= nx; i++)
        for (j = 0; j <= ny; j++)
            m += get_all_diag(nx, ny, i, j);

    return m;
}

int allocate_msr_matrix(
    int nx,
    int ny,
    double **p_A,
    int **p_I)
{
    int diag_len = (nx + 1) * (ny + 1);
    int all_diag = get_len_msr_all_diag(nx, ny);
    int len = diag_len + all_diag + 1;

    double *A = nullptr;
    int *I = nullptr;

    A = new double[len];
    I = new int[len];
    if (!A || !I)
    {
        if (A) delete[] A;
        if (I) delete[] I;
        return 1;
    }

    *p_A = A;
    *p_I = I;
    return 0;
}

void fill_I(int nx, int ny, int *I)
{
    int n = (nx + 1) * (ny + 1);
    int l, i, j;
    int m, r = n + 1;
    for (l = 0; l < n; ++l)
    {
        l2ij(nx, ny, i, j, l);
        I[l] = r;
        m = get_all_diag(nx, ny, i, j, I + r);
        r += m;
    }
    I[l] = r;
}

void fill_A_ij(
    int nx,
    int ny,
    double hx,
    double hy,
    int i,
    int j,
    double *A_diag,
    double *A_all_diag)
{
    double s = hx * hy;
    int l;

    if (i > 0 && i < nx && j > 0 && j < ny)
    {
        *A_diag = 6 * s / 12;
        for (l = 0; l < 6; ++l)
            A_all_diag[l] = s / 12;
    }
    else if (j == 0 && i > 0 && i < nx)
    {
        *A_diag = 3 * s / 12;
        A_all_diag[0] = 1 * s / 24;
        A_all_diag[1] = 1 * s / 24;
        A_all_diag[2] = 2 * s / 24;
        A_all_diag[3] = 2 * s / 24;
    }
    else if (j == ny && i > 0 && i < nx)
    {
        *A_diag = 3 * s / 12;
        A_all_diag[0] = 1 * s / 24;
        A_all_diag[1] = 2 * s / 24;
        A_all_diag[2] = 2 * s / 24;
        A_all_diag[3] = 1 * s / 24;
    }
    else if (i == 0 && j > 0 && j < ny)
    {
        *A_diag = 3 * s / 12;
        A_all_diag[0] = 2 * s / 24;
        A_all_diag[1] = 1 * s / 24;
        A_all_diag[2] = 1 * s / 24;
        A_all_diag[3] = 2 * s / 24;
    }
    else if (i == nx && j > 0 && j < ny)
    {
        *A_diag = 3 * s / 12;
        A_all_diag[0] = 1 * s / 24;
        A_all_diag[1] = 2 * s / 24;
        A_all_diag[2] = 2 * s / 24;
        A_all_diag[3] = 1 * s / 24;
    }
    else if (i == 0 && j == 0)
    {
        *A_diag = 2 * s / 12;
        A_all_diag[0] = 1 * s / 24;
        A_all_diag[1] = 1 * s / 24;
        A_all_diag[2] = 2 * s / 24;
    }
    else if (i == nx && j == ny)
    {
        *A_diag = 2 * s / 12;
        A_all_diag[0] = 1 * s / 24;
        A_all_diag[1] = 2 * s / 24;
        A_all_diag[2] = 1 * s / 24;
    }
    else if (i == 0 && j == ny)
    {
        *A_diag = 1 * s / 12;
        A_all_diag[0] = 1 * s / 24;
        A_all_diag[1] = 1 * s / 24;
    }
    else if (i == nx && j == 0)
    {
        *A_diag = 1 * s / 12;
        A_all_diag[0] = 1 * s / 24;
        A_all_diag[1] = 1 * s / 24;
    }
}

void fill_A(
    int nx,
    int ny,
    double hx,
    double hy,
    int *I,
    double *A,
    int p,
    int pi)
{
    int i, j, l, i1, i2;
    int n = (nx + 1) * (ny + 1);

    thread_rows(n, p, pi, i1, i2);

    for (l = i1; l < i2; ++l)
    {
        l2ij(nx, ny, i, j, l);
        double *A_diag = A + l;
        double *A_all_diag = A + I[l];
        fill_A_ij(nx, ny, hx, hy, i, j, A_diag, A_all_diag);
    }

    synchronize(p);
}

#define FF(I, J) (f(x0 + (I) * hx, y0 + (J) * hy))

double F_IJ(
    int nx,
    int ny,
    double hx,
    double hy,
    double x0,
    double y0,
    int i,
    int j,
    double (*f)(double, double))
{
    double w = hx * hy / 192; // ввели вес

    if (i > 0 && i < nx && j > 0 && j < ny) // внутренние
    {
        return w * (
            36 * FF(i, j)
            + 20 * (
                FF(i + 0.5, j)
                + FF(i, j - 0.5)
                + FF(i - 0.5, j - 0.5)
                + FF(i - 0.5, j)
                + FF(i, j + 0.5)
                + FF(i + 0.5, j + 0.5)
            ) + 4 * (
                FF(i + 0.5, j - 0.5)
                + FF(i - 0.5, j - 1)
                + FF(i - 1, j - 0.5)
                + FF(i - 0.5, j + 0.5)
                + FF(i + 0.5, j + 1)
                + FF(i + 1, j + 0.5)
            ) + 2 * (
                FF(i + 1, j)
                + FF(i, j - 1)
                + FF(i - 1, j - 1)
                + FF(i - 1, j)
                + FF(i, j + 1)
                + FF(i + 1, j + 1)
            )
        );
    }
    if (i > 0 && i < nx && j == 0) // нижняя сторона
    {
        return w * (
            18 * FF(i, j)
            + 10 * (FF(i + 0.5, j) + FF(i - 0.5, j))
            + 20 * (FF(i, j + 0.5) + FF(i + 0.5, j + 0.5))
            + 4 * (
                FF(i - 0.5, j + 0.5)
                + FF(i + 0.5, j + 1)
                + FF(i + 1, j + 0.5)
            )
            + 1 * (FF(i - 1, j) + FF(i + 1, j))
            + 2 * (FF(i, j + 1) + FF(i + 1, j + 1))
        );
    }
    if (i > 0 && i < nx && j == ny)  // верхняя сторона
    {
        return w * (
            18 * FF(i, j)
            + 10 * (FF(i - 0.5, j) + FF(i + 0.5, j))
            + 20 * (FF(i, j - 0.5) + FF(i - 0.5, j - 0.5))
            + 4 * (
                FF(i + 0.5, j - 0.5)
                + FF(i - 0.5, j - 1)
                + FF(i - 1, j - 0.5)
            ) 
            + 1 * (FF(i - 1, j) + FF(i + 1, j))
            + 2 * (FF(i, j - 1) + FF(i - 1, j - 1))
        );
    }
    if (i == 0 && j > 0 && j < ny) // левая сторона без углов
    {
        return w * (
            18 * FF(i, j)
            + 10 * (FF(i, j - 0.5) + FF(i, j + 0.5))
            + 20 * (FF(i + 0.5, j) + FF(i + 0.5, j + 0.5))
            + 4 * (
                FF(i + 0.5, j - 0.5)
                + FF(i + 0.5, j + 1)
                + FF(i + 1, j + 0.5)
            )
            + 1 * (FF(i, j - 1) + FF(i, j + 1))
            + 2 * (FF(i + 1, j) + FF(i + 1, j + 1))
        );
    }
    if (i == nx && j > 0 && j < ny) // правая сторона без углов
    {
        return w * (
            18 * FF(i, j)
            + 10 * (FF(i, j - 0.5) + FF(i, j + 0.5))
            + 20 * (FF(i - 0.5, j) + FF(i - 0.5, j - 0.5))
            + 4 * (
                FF(i - 0.5, j - 1)
                + FF(i - 1, j - 0.5)
                + FF(i - 0.5, j + 0.5)
            )
            + 1 * (FF(i, j - 1) + FF(i, j + 1))
            + 2 * (FF(i - 1, j) + FF(i - 1, j - 1))
        );
    }
    if (i == 0 && j == 0) // левый нижний угол
    {
        return w * (
            12 * FF(i, j)
            + 10 * (FF(i + 0.5, j) + FF(i, j + 0.5))
            + 20 * (FF(i + 0.5, j + 0.5))
            + 4 * (FF(i + 1, j + 0.5) + FF(i + 0.5, j + 1))
            + 1 * (FF(i + 1, j) + FF(i, j + 1))
            + 2 * (FF(i + 1, j + 1))
        );
    }
    if (i == nx && j == ny) // парвый верхний угол
    {
        return w * (
            12 * FF(i, j)
            + 10 * (FF(i - 0.5, j)+ FF(i, j - 0.5))
            + 20 * (FF(i - 0.5, j - 0.5))
            + 4 * (FF(i - 0.5, j - 1) + FF(i - 1, j - 0.5))
            + 1 * (FF(i, j - 1) + FF(i - 1, j))
            + 2 * (FF(i - 1, j - 1))
        );
    }
    if (i == 0 && j == ny) // левый верхний угол
    {
        return w * (
            6 * FF(i, j)
            + 10 * (FF(i + 0.5, j) + FF(i, j - 0.5))
            + 4 * (FF(i + 0.5, j - 0.5))
            + 1 * (FF(i + 1, j) + FF(i, j - 1))
        );
    }
    if (i == nx && j == 0) // правый нижний угол
    {
        return w * (
            6 * FF(i, j)
            + 10 * (FF(i - 0.5, j) + FF(i, j + 0.5))
            + 4 * (FF(i - 0.5, j + 0.5))
            + 1 * (FF(i - 1, j) + FF(i, j + 1))
        );
    }
    return 1e308; // типа ошибка. Сюда не должно быть попаданий
}

void fill_B(
    int nx,
    int ny,
    double hx,
    double hy,
    double x0,
    double y0,
    double *b,
    double (*f)(double, double),
    int p,
    int pi)
{
    int n = (nx + 1) * (ny + 1);
    int i1, i2, i, j, l;

    thread_rows(n, p, pi, i1, i2);

    for (l = i1; l < i2; ++l)
    {
        l2ij(nx, ny, i, j, l);
        b[l] = F_IJ(nx, ny, hx, hy, x0, y0, i, j, f);
    }

    synchronize(p);
}

double r1(
    int nx,
    int ny,
    double hx,
    double hy,
    double x0,
    double y0,
    double *x,
    double (*f)(double, double),
    int p,
    int pi)
{
    int n = (nx + 1) * (ny + 1);
    int i1, i2, i, j, l;
    double res = -1.;

    thread_rows(n, p, pi, i1, i2);

    for (l = i1; l < i2; ++l)
    {
        l2ij(nx, ny, i, j, l);

        if (i == nx || j == ny)
        {
            continue;
        }

        res = std::max(
            res,
            std::abs(
                f(x0 + (i + 2. / 3) * hx, y0 + (j + 1. / 3) * hy)
                - (x[l] + x[l + 1] + x[l + 1 + nx + 1]) / 3
            )
        );
        res = std::max(
            res,
            std::abs(
                f(x0 + (i + 1. / 3) * hx, y0 + (j + 2. / 3) * hy)
                - (x[l] + x[l + nx + 1] + x[l + 1 + nx + 1]) / 3
            )
        );
    }

    sync_max(p, &res, 1);
    return res;
}

double r2(
    int nx,
    int ny,
    double hx,
    double hy,
    double x0,
    double y0,
    double *x,
    double (*f)(double, double),
    int p,
    int pi)
{
    int n = (nx + 1) * (ny + 1);
    int i1, i2, i, j, l;
    double res = 0;

    thread_rows(n, p, pi, i1, i2);

    for (l = i1; l < i2; ++l)
    {
        l2ij(nx, ny, i, j, l);

        if (i == nx || j == ny)
        {
            continue;
        }

        res += std::abs(
            f(x0 + (i + 2. / 3) * hx, y0 + (j + 1. / 3) * hy)
            - (x[l] + x[l + 1] + x[l + 1 + nx + 1]) / 3
        );
        res += std::abs(
            f(x0 + (i + 1. / 3) * hx, y0 + (j + 2. / 3) * hy)
            - (x[l] + x[l + nx + 1] + x[l + 1 + nx + 1]) / 3
        );
    }

    res = reduce_sum_det(p, pi, res);
    return res * hx * hy / 2;
}

double r3(
    int nx,
    int ny,
    double hx,
    double hy,
    double x0,
    double y0,
    double *x,
    double (*f)(double, double),
    int p,
    int pi)
{
    int n = (nx + 1) * (ny + 1);
    int i1, i2, i, j, l;
    double res = -1;

    thread_rows(n, p, pi, i1, i2);

    for (l = i1; l < i2; ++l)
    {
        l2ij(nx, ny, i, j, l);

        res = std::max(res, std::abs(f(x0 + i * hx, y0 + j * hy) - x[l]));
    }

    sync_max(p, &res, 1);
    return res;
}

double r4(
    int nx,
    int ny,
    double hx,
    double hy,
    double x0,
    double y0,
    double *x,
    double (*f)(double, double),
    int p,
    int pi)
{
    int n = (nx + 1) * (ny + 1);
    int i1, i2, i, j, l;
    double res = 0;

    thread_rows(n, p, pi, i1, i2);

    for (l = i1; l < i2; ++l)
    {
        l2ij(nx, ny, i, j, l);
        res += std::abs(f(x0 + i * hx, y0 + j * hy) - x[l]);
    }

    res = reduce_sum_det(p, pi, res);
    return res * hx * hy;
}


double GetCpuTime() {
	struct rusage buf;
	getrusage(RUSAGE_THREAD, &buf);
	return buf.ru_utime.tv_sec + buf.ru_utime.tv_usec / 1.e6;
}

void* thread_func(void* args) {
    thread_data* arg = (thread_data*) args;
    
    double a = arg->a;
    double b = arg->b;
    double c = arg->c;
    double d = arg->d;
    int nx = arg->nx;
    int ny = arg->ny;
	double (*f)(double, double) = arg->f;
	double eps = arg->eps;
	int maxit = arg->maxit;
	int p = arg->p;
	int pi = arg->pi;
	
	double* A = arg->A;
	double* B = arg->B;
	double* x = arg->x;
	int* I = arg->I;
	
	double* r = arg->r;
	double* u = arg->u;
	double* v = arg->v;
    
    cpu_set_t cpu;
    CPU_ZERO(&cpu);
    int n_cpus = get_nprocs();
    int cpu_id = n_cpus - 1 - (pi % n_cpus);
    CPU_SET(cpu_id, &cpu);
    pthread_t tid = pthread_self();

    pthread_setaffinity_np(tid, sizeof(cpu), &cpu);
    
    int n = (nx + 1) * (ny + 1);
	double hx = (b - a) / nx;
    double hy = (d - c) / ny;
	
	fill_A(nx, ny, hx, hy, I, A, p, pi);
    fill_B(nx, ny, hx, hy, a, c, B, f, p, pi); 
	
	int maxsteps = 300;

	arg->t1 = GetCpuTime();
	arg->it = minimal_residual_msr_matrix_full(n, A, I, B, x, r, u, v, eps, maxit, maxsteps, p, pi);
	arg->t1 = GetCpuTime() - arg->t1;
	
	arg->t2 = GetCpuTime();
	arg->r1 = r1(nx, ny, hx, hy, a, c, x, f, p, pi);
	arg->r2 = r2(nx, ny, hx, hy, a, c, x, f, p, pi);
	arg->r3 = r3(nx, ny, hx, hy, a, c, x, f, p, pi);
	arg->r4 = r4(nx, ny, hx, hy, a, c, x, f, p, pi);
	arg->t2 = GetCpuTime() - arg->t2;

	return nullptr;
}
