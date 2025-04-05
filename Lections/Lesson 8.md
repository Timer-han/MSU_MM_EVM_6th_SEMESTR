Можно взять сетку мельче и посчитать интеграл по ней.
Рассмотрели треугольник, нарисовали средние линии, вычислили скалярное произведение
$$
(\varphi A, f)=\frac{1}{192}h_x h_y (6 f_A + f_B + f_C + 10 f_D + 4f_E + 10f_F)
$$
![[IMG_9821.jpg]]
![[IMG_9822.jpg]]
Из этого хорошо видно, что объём пирамиды равен 1, а значит всё хорошо
``` cpp
#define F(I, J) (f(a + (I) * h_x, c + (J) * h_y))

double F_ij (
	int nx,
	int ny,
	double hx,	
	double hy,
	double a,
	double c,
	int i,	
	int j,
	double (*f) (double x, double y)	
) {
	double w = hx * hy / 192; // ввели вес
	if (i > 0 && i < nx && i > 0 && i < ny) {
		return w * (36 * F(i, j) 
				+ 20 * 
				(
					F(i + 0.5, j) + F(i, j - 0.5)
					+ F(i - 0.5, j - 0.5) + F(i - 0.5, j)
					+ F(i, j + 0.5) + F(i + 0.5, j + 0.5)
				)
				+ 4 *
				(
					F(i + 0.5, j - 0.5) + F(i - 0.5, j)
					+ F(i - 1, j - 0.5) + F(i - 0.5, j + 0.5)
					+ F(i + 0.5, j + 1) + F(i + 1, j + 0.5)
				)
				+ 2 *
				(
					F(i + 1, j) + F(i, j - 1)
					+ F(i - 1, j - 1) + F(i - 1, j)
					+ F(i, j + 1) + F(i + 1, j + 1)
				)
			);
	}

	if (i > 0 && i < n && j == 0) { //  нижняя сторона
		return w *
		(
			18 * F(i, j) 
		)
	}
}
```