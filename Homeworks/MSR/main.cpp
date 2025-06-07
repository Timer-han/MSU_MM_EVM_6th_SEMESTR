#include <fstream>
#include <cmath>
#include <pthread.h>  
#include <sys/resource.h>
#include <sys/sysinfo.h>
#include <sys/time.h>

#include "reduce.h"
#include "functions.h"


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

int main(int argc, char *argv[])
{
    double a, b, c, d, eps;
    int nx, ny, pi, maxit, p;

    if (argc < 11)
    {
        printf("[-] Not enough arguments!\n");
        printf("[+] Usage: %s a b c d nx ny pi eps maxit p\n", argv[0]);
        return 1;
    }

    if (
        sscanf(argv[1], "%lf", &a) != 1 ||
        sscanf(argv[2], "%lf", &b) != 1 ||
        sscanf(argv[3], "%lf", &c) != 1 ||
        sscanf(argv[4], "%lf", &d) != 1 ||
        sscanf(argv[5], "%d", &nx) != 1 ||
        sscanf(argv[6], "%d", &ny) != 1 ||
        sscanf(argv[7], "%d", &pi) != 1 ||
        sscanf(argv[8], "%lf", &eps) != 1 ||
        sscanf(argv[9], "%d", &maxit) != 1 ||
        sscanf(argv[10], "%d", &p) != 1)
    {
        printf("[-] Invalid arguments!\n");
        return 2;
    }

    init_reduce_sum(p);

    double (*f)(double, double);
    select_func(pi, f);

    int n = (nx + 1) * (ny + 1);

    double *A;
    int *I;
    if (allocate_msr_matrix(nx, ny, &A, &I) != 0)
    {
        printf("[-] Failed to allocate MSR matrix!\n");
        return -1;
    }
    fill_I(nx, ny, I);
    double *B = new double[n];
    double *x = new double[n];
    double *r = new double[n];
    double *u = new double[n];
    double *v = new double[n];

    if (!A || !I || !B || !x || !r || !u || !v)
    {
        printf("[-] Memory allocation failed!\n");
        delete_reduce_sum();
        if (A)
            delete[] A;
        if (I)
            delete[] I;
        if (B)
            delete[] B;
        if (x)
            delete[] x;
        if (r)
            delete[] r;
        if (u)
            delete[] u;
        if (v)
            delete[] v;
        return 3;
    }

    for (int i = 0; i < n; ++i)
    {
        x[i] = 0;
        r[i] = 0;
        u[i] = 0;
        v[i] = 0;
    }

    thread_data *data = new thread_data[p];
    pthread_t *tid = new pthread_t[p];

    for (int i = 0; i < p; ++i)
    {
        data[i].a = a;
        data[i].b = b;
        data[i].c = c;
        data[i].d = d;
        data[i].nx = nx;
        data[i].ny = ny;
        data[i].f = f;
        data[i].eps = eps;
        data[i].maxit = maxit;
        data[i].p = p;
        data[i].pi = i;

        data[i].A = A;
        data[i].I = I;
        data[i].x = x;
        data[i].B = B;

        data[i].r = r;
        data[i].u = u;
        data[i].v = v;
    }

    for (int i = 1; i < p; ++i)
    {
        pthread_create(&tid[i], 0, thread_func, data + i);
    }
    thread_func(data + 0);

    for (int i = 1; i < p; ++i)
    {
        pthread_join(tid[i], 0);
    }

    int task = 1;
    double t1 = 0;
    double t2 = 0;

    for (int i = 0; i < p; ++i)
    {
        t1 = std::max(t1, data[i].t1);
        t2 = std::max(t2, data[i].t2);
    }

    printf("%s : Task = %d R1 = %e R2 = %e R3 = %e R4 = %e T1 = %.2f T2 = %.2f It = %d E = %e K = %d Nx = %d Ny = %d P = %d\n",
           argv[0], task, data->r1, data->r2, data->r3, data->r4, t1, t2, data->it, eps, pi, nx, ny, p);

    delete_reduce_sum();
    delete[] A;
    delete[] I;
    delete[] B;
    delete[] x;
    delete[] r;
    delete[] u;
    delete[] v;
    delete[] data;
    delete[] tid;

    return 0;
}
