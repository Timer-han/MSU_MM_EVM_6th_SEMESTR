#include <fstream>
#include <cmath>
#include <pthread.h>  
#include <sys/resource.h>
#include <sys/sysinfo.h>
#include <sys/time.h>

#include "reduce.h"
#include "functions.h"


double f_0([[maybe_unused]] double x, [[maybe_unused]] double y) {
	return 1;
}

double f_1(double x, [[maybe_unused]] double y) {
	return x;
}

double f_2([[maybe_unused]] double x, double y) {
	return y;
}

double f_3(double x, double y) {
	return x + y;
}

double f_4(double x, double y) {
	return sqrt(x * x + y * y);
}

double f_5(double x, double y) {
	return x * x + y * y;
}

double f_6(double x, double y) {
	return exp(x * x - y * y);
}

double f_7(double x, double y) {
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

double GetCpuTime() {
	struct rusage buf;
	getrusage(RUSAGE_THREAD, &buf);
	return buf.ru_utime.tv_sec + buf.ru_utime.tv_usec / 1.e6;
}

class thread_data {
public:
    double a;
    double b;
    double c;
    double d;
    int nx;
    int ny;
	double (*f)(double, double);
	double eps;
	int maxit;
	int p;
	int k;
	
	double* A;
	int* I;
	double* x;
	double* B;
	
	double* r;
	double* u;
	double* v;
	
	double r1;
	double r2;
	double r3;
	double r4;
	
	double t1;
	double t2;
	
	int it;
};

void* thread_func(void* arg) {
    thread_data* data = (thread_data*) arg;
    
    double a = data->a;
    double b = data->b;
    double c = data->c;
    double d = data->d;
    int nx = data->nx;
    int ny = data->ny;
	double (*f)(double, double) = data->f;
	double eps = data->eps;
	int maxit = data->maxit;
	int p = data->p;
	int k = data->k;
	
	double* A = data->A;
	int* I = data->I;
	double* x = data->x;
	double* B = data->B;
	
	double* r = data->r;
	double* u = data->u;
	double* v = data->v;
    
    cpu_set_t cpu;
    CPU_ZERO(&cpu);
    int n_cpus = get_nprocs();
    int cpu_id = n_cpus - 1 - (k % n_cpus);
    CPU_SET(cpu_id, &cpu);
    pthread_t tid = pthread_self();

    pthread_setaffinity_np(tid, sizeof(cpu), &cpu);
    
    int n = (nx + 1) * (ny + 1);
	double hx = (b - a) / nx;
    double hy = (d - c) / ny;
	int maxsteps = 300;
	
	fill_A(nx, ny, hx, hy, I, A, p, k);
    fill_B(nx, ny, hx, hy, a, c, B, f, p, k); 
	
	data->t1 = GetCpuTime();
	data->it = minimal_residual_msr_matrix_full(n, A, I, B, x, r, u, v, eps, maxit, maxsteps, p, k);
	data->t1 = GetCpuTime() - data->t1;
	
	data->t2 = GetCpuTime();
	data->r1 = r1(nx, ny, hx, hy, a, c, x, f, p, k);
	data->r2 = r2(nx, ny, hx, hy, a, c, x, f, p, k);
	data->r3 = r3(nx, ny, hx, hy, a, c, x, f, p, k);
	data->r4 = r4(nx, ny, hx, hy, a, c, x, f, p, k);
	data->t2 = GetCpuTime() - data->t2;

	return nullptr;
}

int main(int argc, char* argv[]) {
    argc = argc;
    double a = std::stod(argv[1]);
    double b = std::stod(argv[2]);
    double c = std::stod(argv[3]);
    double d = std::stod(argv[4]);
    int nx = std::stoi(argv[5]);
    int ny = std::stoi(argv[6]);
	int k = std::stoi(argv[7]);
	double eps = std::stod(argv[8]);
	int maxit = std::stoi(argv[9]);
	int p = std::stoi(argv[10]);
	
	InitReduceSum(p);
	
	double (*f)(double, double);
	select_func(k, f);
	
    int n = (nx + 1) * (ny + 1);
	
	double* A;
	int* I;
	if (allocate_msr_matrix(nx, ny, &A,  &I) != 0) {
		return -1;
	}
    fill_I(nx, ny, I);
	double* B = new double[n];
	double* x = new double[n];
	double* r = new double[n];
	double* u = new double[n];
	double* v = new double[n];
	
	for (int i = 0; i < n; ++i) {
		x[i] = 0;
		r[i] = 0;
		u[i] = 0;
		v[i] = 0;
	}
    
    thread_data* data_arr = new thread_data[p];
    pthread_t* tid = new pthread_t[p];
    for (int i = 0; i < p; ++i) {
        data_arr[i].a = a;
        data_arr[i].b = b;
        data_arr[i].c = c;
        data_arr[i].d = d;
		
		data_arr[i].nx = nx;
		data_arr[i].ny = ny;
		
		data_arr[i].f = f;
		data_arr[i].eps = eps;
		data_arr[i].maxit = maxit;
		data_arr[i].p = p;
		data_arr[i].k = i;
		
		data_arr[i].A = A;
		data_arr[i].I = I;
		data_arr[i].x = x;
		data_arr[i].B = B;
		
		data_arr[i].r = r;
		data_arr[i].u = u;
		data_arr[i].v = v;
    }
    
    for (int i = 1; i < p; ++i) {
        pthread_create(&tid[i], 0, thread_func, data_arr + i);
    }
    thread_func(data_arr + 0);
	
    for (int i = 1; i < p; ++i) {
        pthread_join(tid[i], 0);
    }
    
    int task = 1;
	double t1 = 0;
	double t2 = 0;
	
	for (int i = 0; i < p; ++i) {
        t1 = std::max(t1, data_arr[i].t1);
        t2 = std::max(t2, data_arr[i].t2);
    }
    
    printf ("%s : Task = %d R1 = %e R2 = %e R3 = %e R4 = %e T1 = %.2f T2 = %.2f It = %d E = %e K = %d Nx = %d Ny = %d P = %d\n", 
			argv[0], task, data_arr->r1, data_arr->r2, data_arr->r3, data_arr->r4, t1, t2, data_arr->it, eps, k, nx, ny, p);
    
	FreeReduceSum();
	delete[] A;
	delete[] I;
	delete[] B;
	delete[] x;
	delete[] r;
	delete[] u;
	delete[] v;
	delete[] data_arr;
    delete[] tid;
	
    return 0;
}
