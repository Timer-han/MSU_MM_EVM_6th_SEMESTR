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
    int pi
);
int minimal_residual_msr_matrix_full(
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
    int maxsteps,
    int p,
    int pi
);
void thread_rows(int n, int p, int pi, int &i1, int &i2);

double scalar_product(int n, double *x, double *y, int p, int pi);
void mult_sub_vector(int n, double *x, double *y, double t, int p, int pi);

void apply_preconditioner_msr_matrix(
    int n,
    double *A,
    int * /*I*/,
    double *v,
    double *r,
    int p,
    int pi
);
void mult_msr_matrix_vector(
    int n,
    double *A,
    int *I,
    double *x,
    double *y,
    int p,
    int pi
);

void ij2l(int nx, int /*ny*/, int i, int j, int &l);
void l2ij(int nx, int /*ny*/, int &i, int &j, int l);

int get_len_msr(int nx, int ny);
int get_all_diag(int nx, int ny, int i, int j, int *I_ij = nullptr);
int get_len_msr_all_diag(int nx, int ny);

int allocate_msr_matrix(int nx, int ny, double **p_A, int **p_I);

void fill_I(int nx, int ny, int *I);
void fill_A_ij(
    int nx,
    int ny,
    double hx,
    double hy,
    int i,
    int j,
    double *A_diag,
    double *A_off_diag
);
void fill_A(
    int nx,
    int ny,
    double hx,
    double hy,
    int *I,
    double *A,
    int p,
    int pi
);


double F_IJ(
	int nx,
	int ny,
	double hx,
	double hy,
	double x0,
	double y0,
	int i,
	int j,
	double (*f)(double, double),
	int p_mid,
	double f_abs
);
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
    int pi,
	int p_mid,
	double f_abs
);


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
    int pi
);
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
    int pi
);
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
    int pi
);
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
    int pi
);

double GetCpuTime ();
void* thread_func (void* args);

class thread_data {
public:
	int p_mid;
	double f_abs;

    double a;
    double b;
    double c;
    double d;
    int nx;
    int ny;
	double (*f)(double, double);
	double eps;
	int max_it;
	int p;
	int pi;
	
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
	
	pthread_cond_t* cond;
	pthread_mutex_t* mutex;
	bool ready;
	bool finish;
};
