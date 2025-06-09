#include <cmath>
#include <pthread.h>
#include <sys/resource.h>
#include <sys/sysinfo.h>
#include <sys/time.h>
#include "functions.h"
#include "reduce.h"

int minimal_residual_msr_matrix(int n, double* A, int* I, double* b, double* x, double* r, double* u, double* v, double eps, int maxit, int p, int k) {
    double prec;
    double b_norm2;
    double tau;
    double c1;
    double c2;
    int it;
    b_norm2 = scalar_product(n, b, b, p, k); // (b, b)
    prec = b_norm2 * eps * eps;
    multiply_msr_matrix_vector(n, A, I, x, r, p, k); // r=Ax
    mult_sub_vector(n, r, b, 1., p, k); // r-=1*b
    for (it = 0; it < maxit; ++it) {
        apply_preconditioner_msr_matrix(n, A, I, v, r, p, k);
        multiply_msr_matrix_vector(n, A, I, v, u, p, k); // u=Av
        c1 = scalar_product(n, u, r, p, k); // c1=(u,r)
        c2 = scalar_product(n, u, u, p, k); // c2=(u,u)
        
        if (c1 < prec || c2 < prec) {
            break;
        }
        
        tau = c1 / c2;
        mult_sub_vector(n, x, v, tau, p, k); // x-=tau*v
        mult_sub_vector(n, r, u, tau, p, k); // r-=tau*u
    }
    if (it > maxit) {
        return -1;
    }
    return it;
}

int minimal_residual_msr_matrix_full(int n, double* A, int* I, double* b, double* x, double* r, double* u, double* v, double eps, int maxit, int maxsteps, int p, int k) {
    int step;
    int ret;
    int its = 0;
    
    for (step = 0; step < maxsteps; ++step) {
        ret = minimal_residual_msr_matrix(n, A, I, b, x, r, u, v, eps, maxit, p, k);
        if (ret >= 0) {
            its += ret;
            break;
        }
        its += maxit;
    }
    if (step >= maxsteps) {
        return -1;
    }
    return its;
}

void thread_rows(int n , int p, int k, int &i1, int &i2) {
    i1 = n * k, 
    i1 /= p;
    i2 = n * (k + 1);
    i2 /= p;
}

double scalar_product(int n, double* x, double* y, int p, int k) {
	int i1;
	int i2;
	int i;
	double s = 0;
	thread_rows(n, p, k, i1, i2);
	for (i = i1; i < i2; ++i) {
		s += x[i] * y[i];
	}
	s = ReduceSumDet(p, k, s);
	return s;
}

void mult_sub_vector(int n, double* x, double* y, double t, int p, int k) {
	int i;
	int i1;
	int i2;
	thread_rows(n, p, k, i1, i2);
	for (i = i1; i < i2; ++i) {
		x[i] -= t * y[i];
	}
	ReduceSum(p);
}

void apply_preconditioner_msr_matrix(int n, double* A, int* /*I*/, double* v, double* r, int p, int k) {
	int i;
	int i1;
	int i2;
	thread_rows(n, p, k, i1, i2);
	for (i = i1; i < i2; ++i) {
		v[i] = r[i] / A[i];
	}
	ReduceSum(p);
}

void multiply_msr_matrix_vector(int n, double* A, int* I, double* x, double* y, int p, int k) {
    int i;
	int i1;
	int i2;
	int l;
	int J;
	double s;
	thread_rows(n, p, k, i1, i2);
	for (i = i1; i < i2; ++i) {
		s = A[i] * x[i]; //диагональ
		l = I[i + 1] - I[i]; // число элементов
		J = I[i];
		for (int j = 0; j < l; ++j) {
			s += A[J + j] * x[I[J + j]];
		}
		y[i] = s;
	}
	ReduceSum(p);
}

void ij2l(int nx, int /*ny*/, int i, int j, int &l) {
	l = i + j * (nx + 1);
}

void l2ij(int nx, int /*ny*/, int &i, int &j, int l) {
	j = l / (nx + 1);
	i = l - j * (nx + 1);
}

int get_len_msr(int nx, int ny) {
	return (nx + 1) * (ny + 1) + 6 * (nx - 1) * (ny - 1) + 4 * (2 * (nx - 1) + 2 * (ny - 1)) + 2 * 3 + 2 * 2;
}

#define F(I, J) do {ij2l(nx, ny, I, J, k); if (I_ij) {I_ij[m] = k;} ++m;} while(0)

int get_off_diag(int nx, int ny, int i, int j, int* I_ij) {
	int m = 0;
	int k;
	if (i < nx) {
		F(i + 1, j);
	}
	if (j > 0) {
		F(i, j - 1);
	}
	if (i > 0 && j > 0)  {
		F(i - 1, j - 1);
	}
	if (i > 0)  {
		F(i - 1, j);
	}
	if (j < ny)  {
		F(i, j + 1);
	}
	if (i < nx && j < ny)  {
		F(i + 1, j + 1);
	}
	return m;
}

int get_len_msr_off_diag(int nx, int ny) {
	int m = 0;
	int i;
	int j;
	for (i = 0; i <= nx; ++i) {
		for (j = 0; j <= ny; ++j) {
			m += get_off_diag(nx, ny, i, j);
		}
	}
	return m;
}

int allocate_msr_matrix(int nx, int ny, double** p_A, int** p_I) {
	int diag_len = (nx + 1) * (ny + 1);
	int off_diag = get_len_msr_off_diag(nx, ny);
	int len = diag_len + off_diag + 1;
	double* A = nullptr;
	int* I = nullptr;
	A = new double[len];
	if (A == nullptr) {
		return 1;
	}
	I = new int[len];
	if (I == nullptr) {
		return 2;
	}
	*p_A = A;
	*p_I = I;
	return 0;
}

void fill_I(int nx, int ny, int* I) {
	int n = (nx + 1) * (ny + 1);
	int l;
	int i;
	int j;
	int m;
	int r = n + 1;
	for (l = 0; l < n; ++l) {
		l2ij(nx, ny, i, j, l); 
		I[l] = r;
		m = get_off_diag(nx, ny, i, j, I + r);
		r += m;
	}
	I[l] = r;
}

void fill_A_ij(int nx, int ny, double hx, double hy, int i, int j, double* A_diag, double* A_off_diag) {
	double s = hx * hy;
	int l;
	
	if (i > 0 && i < nx && j > 0 && j < ny) {
		*A_diag = 6 * s / 12;
		for (l = 0; l < 6; ++l) {
			A_off_diag[l] = s / 12;
		}
	} 
	else if (j == 0 && i > 0 && i < nx) {
		*A_diag = 3 * s / 12;
		A_off_diag[0] = 1 * s / 24;
		A_off_diag[1] = 1 * s / 24;
		A_off_diag[2] = 2 * s / 24;
		A_off_diag[3] = 2 * s / 24;
	}
	else if (j == ny && i > 0 && i < nx) {
		*A_diag = 3 * s / 12;
		A_off_diag[0] = 1 * s / 24;
		A_off_diag[1] = 2 * s / 24;
		A_off_diag[2] = 2 * s / 24;
		A_off_diag[3] = 1 * s / 24;
	}
	else if (i == 0 && j > 0 && j < ny) {
		*A_diag = 3 * s / 12;
		A_off_diag[0] = 2 * s / 24;
		A_off_diag[1] = 1 * s / 24;
		A_off_diag[2] = 1 * s / 24;
		A_off_diag[3] = 2 * s / 24;
	}
	else if (i == nx && j > 0 && j < ny) {
		*A_diag = 3 * s / 12;
		A_off_diag[0] = 1 * s / 24;
		A_off_diag[1] = 2 * s / 24;
		A_off_diag[2] = 2 * s / 24;
		A_off_diag[3] = 1 * s / 24;
	}
	else if (i == 0 && j == 0) {
		*A_diag = 2 * s / 12;
		A_off_diag[0] = 1 * s / 24;
		A_off_diag[1] = 1 * s / 24;
		A_off_diag[2] = 2 * s / 24;
	}
	else if (i == nx && j == ny) {
		*A_diag = 2 * s / 12;
		A_off_diag[0] = 1 * s / 24;
		A_off_diag[1] = 2 * s / 24;
		A_off_diag[2] = 1 * s / 24;
	}
	else if (i == 0 && j == ny) {
		*A_diag = 1 * s / 12;
		A_off_diag[0] = 1 * s / 24;
		A_off_diag[1] = 1 * s / 24;
	}
	else if (i == nx && j == 0) {
		*A_diag = 1 * s / 12;
		A_off_diag[0] = 1 * s / 24;
		A_off_diag[1] = 1 * s / 24;
	}
}

void fill_A(int nx, int ny, double hx, double hy, int* I, double* A, int p, int k) {
	int i;
	int j;
	int l;
	int i1;
	int i2;
	int n = (nx + 1) * (ny + 1);
	thread_rows(n , p, k, i1, i2);
	
	for (l = i1; l < i2; ++l) {
		l2ij(nx, ny, i, j, l);
		double* A_diag = A + l;
        double* A_off_diag = A + I[l];
        fill_A_ij(nx, ny, hx, hy, i, j, A_diag, A_off_diag);
	}
	
	ReduceSum(p);
}

#define FF(I, J) (f(x0 + (I) * hx, y0 + (J) * hy) + (((std::abs(i -  double(nx / 2)) < 0.25) && (std::abs(j - double(ny / 2)) < 0.25)) ? p_mid * 0.1 * f_abs : 0))

double F_IJ(int nx, int ny, double hx, double hy, double x0, double y0, int i, int j, double (*f)(double, double), int p_mid, double f_abs) {
	double w = hx * hy / 192;
	
	if (i > 0 && i < nx && j > 0 && j < ny) {
		return w * (36 * FF(i, j) 
		+ 20 * (FF(i + 0.5, j) + FF(i, j - 0.5) + FF(i - 0.5, j - 0.5) + FF(i - 0.5, j) + FF(i, j + 0.5) + FF(i + 0.5, j + 0.5))
		+ 4 * (FF(i + 0.5, j - 0.5) + FF(i - 0.5, j - 1) + FF(i - 1,j - 0.5) + FF(i - 0.5, j + 0.5) + FF(i + 0.5, j + 1) + FF(i + 1, j + 0.5))
		+ 2 * (FF(i + 1, j) + FF(i, j - 1) + FF(i - 1, j - 1) + FF(i - 1, j) + FF(i, j + 1) + FF(i + 1, j + 1)));
	}
	if (i > 0 && i < nx && j == 0) {
		return w * (18 * FF(i, j) 
		+ 10 * (FF(i + 0.5, j) + FF(i - 0.5, j))
		+ 20 * (FF(i, j + 0.5) + FF(i + 0.5, j + 0.5))
		+ 4 * (FF(i - 0.5, j + 0.5) + FF(i + 0.5, j + 1) + FF(i + 1, j + 0.5))
		+ 1 * (FF(i - 1, j) + FF(i + 1, j))
		+ 2 * (FF(i, j + 1) + FF(i + 1, j + 1)));
	}
	if (i > 0 && i < nx && j == ny) {
        return w * (18 * FF(i, j) 
		+ 10 * (FF(i - 0.5, j) + FF(i + 0.5, j))
		+ 20 * (FF(i, j - 0.5) + FF(i - 0.5, j - 0.5))
		+ 4 * (FF(i + 0.5,j - 0.5) + FF(i - 0.5, j - 1) + FF(i - 1, j - 0.5)) 
		+ 1 * (FF(i - 1, j) + FF(i + 1, j)) 
		+ 2 * (FF(i, j - 1) + FF(i - 1, j - 1)));
    }
    if (i == 0 && j > 0 && j < ny) {
        return w * (18 * FF(i, j) 
		+ 10 * (FF(i, j - 0.5) + FF(i, j + 0.5)) 
		+ 20 * (FF(i + 0.5, j) + FF(i + 0.5, j + 0.5))
		+ 4 * (FF(i + 0.5, j - 0.5) + FF(i + 0.5, j + 1) + FF(i + 1, j + 0.5)) 
		+ 1 * (FF(i, j - 1) + FF(i, j + 1))
		+ 2 * (FF(i + 1, j) + FF(i + 1, j + 1)));     
    }
    if (i == nx && j > 0 && j < ny) {
        return w * (18 * FF(i, j) 
		+ 10 * (FF(i, j - 0.5) + FF(i, j + 0.5))
		+ 20 * (FF(i - 0.5, j) + FF(i - 0.5, j - 0.5))
		+ 4 * (FF(i - 0.5, j - 1) + FF(i - 1, j - 0.5) + FF(i - 0.5, j + 0.5))
		+ 1 * (FF(i, j - 1) + FF(i, j + 1))
		+ 2 * (FF(i - 1, j) + FF(i - 1, j - 1)));
    }
    if (i == 0 && j == 0) {
        return w * (12 * FF(i, j) 
		+ 10 * (FF(i + 0.5, j) + FF(i, j + 0.5)) 
		+ 20 * (FF(i + 0.5, j + 0.5))
		+ 4 * (FF(i + 1, j + 0.5) + FF(i + 0.5, j + 1))
		+ 1 * (FF(i + 1, j) + FF(i, j + 1))
		+ 2 * (FF(i + 1, j + 1)));     
    }    
    if (i == nx && j == ny) {
        return w * (12 * FF(i, j)
		+ 10 * (FF(i - 0.5, j) + FF(i, j - 0.5))
		+ 20 * (FF(i - 0.5, j - 0.5))
		+ 4 * (FF(i - 0.5, j - 1) + FF(i - 1, j - 0.5))
		+ 1 * (FF(i, j - 1) + FF(i - 1, j))
		+ 2 * (FF(i - 1, j - 1)));   
    }
    if (i == 0 && j == ny) {
        return w * (6 * FF(i, j)
		+ 10 * (FF(i + 0.5, j) + FF(i, j - 0.5))
		+ 4 * (FF(i + 0.5, j - 0.5)) 
		+ 1 * (FF(i+1,j) + FF(i,j-1))); 
    }
    if (i == nx && j == 0) {
        return w * (6 * FF(i, j)
		+ 10 * (FF(i - 0.5, j) + FF(i, j + 0.5))
		+ 4 * (FF(i - 0.5, j + 0.5))
		+ 1 * (FF(i - 1, j) + FF(i, j + 1)));      
    }
    return 1e308;
}

void fill_B(int nx, int ny, double hx, double hy, double x0, double y0, double* b, double (*f)(double, double), int p, int k, int p_mid, double f_abs) {
	int n = (nx + 1) * (ny + 1);
	int i1;
	int i2;
	int l;
	int i;
	int j;
	thread_rows(n, p, k, i1, i2);
	
	for (l = i1; l < i2; ++l) {
        l2ij(nx, ny, i, j, l);
        b[l] = F_IJ(nx, ny, hx, hy, x0, y0, i, j, f, p_mid, f_abs);
    }
    
	ReduceSum(p);
}

double r1(int nx, int ny, double hx, double hy, double x0, double y0, double* x, double (*f)(double, double), int p, int k) {
	int n = (nx + 1) * (ny + 1);
	int i1;
	int i2;
	int l;
	int i;
	int j;
	double r1_residual = -1;
	
    thread_rows(n, p, k, i1, i2);
	
	for (l = i1; l < i2; ++l) {
        l2ij(nx, ny, i, j, l);
		
		if (i == nx || j == ny) {
            continue;
        }
        
		r1_residual = std::fmax(r1_residual, std::abs(f(x0 + (i + 2. / 3) * hx, y0 + (j + 1. / 3) * hy) - (x[l] + x[l + 1] + x[l + 1 + nx + 1]) / 3));
		r1_residual = std::fmax(r1_residual, std::abs(f(x0 + (i + 1. / 3) * hx, y0 + (j + 2. / 3) * hy) - (x[l] + x[l + nx + 1] + x[l + 1 + nx + 1]) / 3));
	}
	
	ReduceMax(p, &r1_residual, 1);
	return r1_residual;
}

double r2(int nx, int ny, double hx, double hy, double x0, double y0, double* x, double (*f)(double, double), int p, int k) {
	int n = (nx + 1) * (ny + 1);
	int i1;
	int i2;
	int l;
	int i;
	int j;
	double r2_residual = 0;
	
    thread_rows(n, p, k, i1, i2);
	
	for (l = i1; l < i2; ++l) {
        l2ij(nx, ny, i, j, l);
		
		if (i == nx || j == ny) {
            continue;
        }
        
		r2_residual += std::abs(f(x0 + (i + 2. / 3) * hx, y0 + (j + 1. / 3) * hy) - (x[l] + x[l + 1] + x[l + 1 + nx + 1]) / 3);
		r2_residual += std::abs(f(x0 + (i + 1. / 3) * hx, y0 + (j + 2. / 3) * hy) - (x[l] + x[l + nx + 1] + x[l + 1 + nx + 1]) / 3);
	}
	
	r2_residual = ReduceSumDet(p, k, r2_residual);
	return r2_residual * hx * hy / 2;;
}

double r3(int nx, int ny, double hx, double hy, double x0, double y0, double* x, double (*f)(double, double), int p, int k) {
	int n = (nx + 1) * (ny + 1);
	int i1;
	int i2;
	int l;
	int i;
	int j;
	double r3_residual = -1;
	
    thread_rows(n, p, k, i1, i2);
	
	for (l = i1; l < i2; ++l) {
        l2ij(nx, ny, i, j, l);
        
		r3_residual = std::fmax(r3_residual, std::abs(f(x0 + i * hx, y0 + j * hy) - x[l]));
	}
	
	ReduceMax(p, &r3_residual, 1);
	return r3_residual;
}

double r4(int nx, int ny, double hx, double hy, double x0, double y0, double* x, double (*f)(double, double), int p, int k) {
	int n = (nx + 1) * (ny + 1);
	int i1;
	int i2;
	int l;
	int i;
	int j;
	double r4_residual = 0;
	
    thread_rows(n, p, k, i1, i2);
	
	for (l = i1; l < i2; ++l) {
        l2ij(nx, ny, i, j, l);
        
		r4_residual += std::abs(f(x0 + i * hx, y0 + j * hy) - x[l]);
	}
	
	r4_residual = ReduceSumDet(p, k, r4_residual);
	return r4_residual * hx * hy;
}

double GetCpuTime() {
	struct rusage buf;
	getrusage(RUSAGE_THREAD, &buf);
	return buf.ru_utime.tv_sec + buf.ru_utime.tv_usec / 1.e6;
}

void* thread_func(void* arg) {
    thread_data* data = (thread_data*) arg;
    pthread_cond_t *cond = data->cond;
    pthread_mutex_t *mutex = data->mutex;
	
    cpu_set_t cpu;
    CPU_ZERO(&cpu);
    int n_cpus = get_nprocs();
    int cpu_id = n_cpus - 1 - (data->k % n_cpus);
    CPU_SET(cpu_id, &cpu);
    pthread_t tid = pthread_self();
	
    pthread_setaffinity_np(tid, sizeof(cpu), &cpu);
	
	while (!(data->finish)) {
		if (!(data->ready)) {
			double a = data->a;
			double b = data->b;
			double c = data->c;
			double d = data->d;
			int nx = data->nx;
			int ny = data->ny;
			double (*f)(double, double) = data->f;
			
			int p_mid = data->p_mid;
			double f_abs = data->f_abs;
			
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
			
			int n = (nx + 1) * (ny + 1);
			double hx = (b - a) / nx;
			double hy = (d - c) / ny;
			int maxsteps = 300;
			
			fill_A(nx, ny, hx, hy, I, A, p, k);
			fill_B(nx, ny, hx, hy, a, c, B, f, p, k, p_mid, f_abs); 
			
			data->t1 = GetCpuTime();
			data->it = minimal_residual_msr_matrix_full(n, A, I, B, x, r, u, v, eps, maxit, maxsteps, p, k);
			data->t1 = GetCpuTime() - data->t1;
			
			data->t2 = GetCpuTime();
			data->r1 = r1(nx, ny, hx, hy, a, c, x, f, p, k);
			data->r2 = r2(nx, ny, hx, hy, a, c, x, f, p, k);
			data->r3 = r3(nx, ny, hx, hy, a, c, x, f, p, k);
			data->r4 = r4(nx, ny, hx, hy, a, c, x, f, p, k);
			data->t2 = GetCpuTime() - data->t2;
			
			data->ready = true;
            ReduceSum(p);

            pthread_mutex_lock(mutex);
            while (data->ready) {
                pthread_cond_wait(cond, mutex);
            }
            pthread_mutex_unlock(mutex);
		}
	}

	return nullptr;
}
