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

Также ответ может меняться от запуска к запуску из-за точек синхронизации. В `reduce_sum` должен быть детерменированный порядок, нужно выделить массив, а потом просуммировать в одном порядке.

``` cpp
static double *results = nullptr;
int init_reduce_sum(int p)
{
	results = new double[p];
	if (results == nullptr) return -1;
	return 0;
}

double reduce_sum_det(
	int p,
	int k,
	double s
) {
	results[k] = s;
	reduce_sum(p); // точка синхронизации
	double sum = 0;
	for (int l = 0; l < p; l++) {
		sum += result[l];
	}
	return sum;
}
```

![[IMG_9742.jpg]]
![[IMG_9743.jpg]]
``` cpp
int ij2l (
	int nx,
	int /*ny*/,
	int i,
	int j,
	int &l
) {
	l = i + j *(nx + 1);
}

int l2ij (
	int nx,
	int /*ny*/,
	int &i,
	int &j,
	int l	
) {
	j = l / (nx + 1);
	i = l - j*(nx + 1);
}

// Число внедиагональных элементов
int get_len_msr (
	int nx,
	int ny
) {
	return 6 * (nx - 1) * (ny - 1) + 4 * (2 * (nx - 1) 
		 + 2 * (ny - 1)) 
		 + 3 * 2 + 2 * 2;
}

// Общая длина msr матрицы = диаг + 1 + число внедиагональных
// (nx + 1) * (ny - 1) + 1 + get_len_msr(nx, ny)
```

![[IMG_9745.jpg]]
```cpp
#define F(IS, JS, S) \
	IA_ij(nx, ny, hx, hy, i, j, (IS), (JS), (S), I, A)

void IA_ij (
	int nx,
	int ny,
	double hx,
	double hy,
	int i,
	int j, // (i,j)
	int is,
	int js,
	int s,
	int *I = nullptr,
	int *A = nullptr
) {

}

//возвращает количество внедиагональных точки (i, j)
int get_all_diag(
	int nx,
	int ny,
	double hx,
	double hy,
	int i,
	int j,
	int *I,
	double *A
) {
	if (i > 0 && i < nx && j > 0 && j < ny) {
		if (I && A) {
			F(i + 1, j, );
			F(i, j, );
			F(i, j, );
			F(i, j, );
			F(i, j, );
			F(i, j, );
		}
	}
}
```