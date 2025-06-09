class thread_data {
public:
    pthread_t tid;
    thread_data *pthread_args;

    double a;
    double b;
    double c;
    double d;
    int nx;
    int ny;
    int func_id;
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

    bool *is_closing;
    bool *is_running;
    bool *is_calculated;


    pthread_mutex_t *mutex;
    pthread_cond_t *cond;
  
    // int is_running = 1;
};


double f_0(double /* x */, double /* y */);
double f_1(double x, double /* y */);
double f_2(double /* x */, double y);
double f_3(double x, double y);
double f_4(double x, double y);
double f_5(double x, double y);
double f_6(double x, double y);
double f_7(double x, double y);

void select_func(int k, double (*&f)(double, double));

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
    double (*f)(double, double)
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
    int pi
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