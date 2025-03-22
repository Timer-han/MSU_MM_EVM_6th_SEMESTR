Метод минимальных невязок
``` cpp
int minim_resid_msr_matrix
(
	int n,
	double *A,
	int *I,
	double *b,
	double *x, /*нач знач, а потом результат*/
	double *r,
	double *u,
	double *v,
	double eps,
	int max_it,
	int p,
	int a
) {
	double prec, b_norm2, tau;
	int it;
	b_norm = scalar_product (n, b, b, p, k);
	prec = b_corm2 * eps * eps;
	// r = Ax
	matrix_mult_vector_msr(n, A, I, x, r, p, k);
	// 1 точка синхронизации
	
}
```