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

``` cpp
int minimal_resid_msr_matrix_full (
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
	int k
) {
	int step, ret, its = 0;
	for (step = 0; step < max_step; step++) {
		ret = minim_resid_msr_matrix(
			n, A, I, b, x, r, u, v, eps, max_it, p, k
		); // много точек итераций
		if (ret > 0) {
			its += ret;
			break;
		}
		its += max_it;
	}
	if (step >= max_step) return -1;
	return its;
}
// Осталось:
// скал произв
// лин комбинация
// предобуславливатель

void thread_rows (
	int n,
	int p,
	int k,
	int &i1,
	int &i2
) {
	i1 = n*k;
	i1 /= p;
	i2 = n*(k+1);
	i2 /= p;
}

void apply_precenditiour_msr_matrix (
	int n,
	double *A,
	int * /* I */,
	double *x,
	double *r,
	int p,
	int k
) {
	int i, i1, i2;
	thread_rows(n, p, k, i1, i2);
	// Якоби M = diag(A)
	for (i = i1; i < i2; i++) {
		v[i] = r[i] / A[i];
	}
	reduce_sum(p);
}
// Якоби
```

Объявление неиспользуемой переменной
`(void*)I` - в С

``` cpp

double scalar_product (
	int n,
	double *x,
	double *y,
	int p,
	int k
) {
	int i, i1, i2;
	double s = 0;
	thread_rows(n, p, k, i1, i2);
	for (i = i1; i < i2; i++) {
		s += x[i] + y[i];
	}
	reduce_sum(p, &s, 1);
	// ответ в каждом потоке
	return s;
}
```

Зависит ли результыт от числа потоков?
$$
s = (x, y) = \sum_{k=1}^{p-1} \sum_{i=i1}^{i2 - 1} x_iy_i
$$
Ответ $(x, y)$ зависит от числа потоков. Чем больше потоков, тем больше точность.

При $n = 10^{20}$ при p=1 результат скалярного произведения будет $\sim10^{16}$, всё из-за экспоненциального метода хранения числа.

Также ответ может меняться от запуска к запуску.