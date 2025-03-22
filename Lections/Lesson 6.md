Метод минимальных невязок

$$
\tau = \frac{()}{()}
$$
``` cpp
int minim_resid_msr_matrix(
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
	int k
) {
	double prec, b_norm2, tau;
	int it;
	b_norm = scalar_product (n, b, b, p, k);
	prec = b_corm2 * eps * eps;
	// r = Ax
	matrix_mult_vector_msr (n, A, I, x, r, p, k);
	// 1 точка синхронизации
	// r = Ax - b, r-=b, r-=1*b
	mult_sub_vector(n, r, b, 1., p, k);
	// r -= 1.*b, 1 точка синхронизации
	for (it = 0; it < max_it; it++) {
		// Mr = r - решить систему
		apply_precenditiour_msr_matrix (
			n, A, I, v, r, p, k
		); // 1 точка синхронизациии
		// u = Av, u = AM^(-1)r
		matrix_mult_vector_msr(n, A, I, v, u, p, k);
		// 1 точка синхронизации
		// c_1 = (u, r)
		c1 = scalar_product(n, u, r, p, k);
		// c2 = (u, u)
		c2 = scalar_product(n, u, u, p, k);
		if (c1 < prec || c2 < prec) {
			break;
		}
		tau = c1/c2;
		x -= tau * v;
		mult_sub_vector(n, x, r, tau, p, k);
		// 1 точка синхронизации
		// r -= tau * u
		mult_sub_vector(n, r, u, tau, p, k);
		// 1 точка синхронизации
	}
	if (it >= max_it) {
		return -1;
	}
	return it;
}
```