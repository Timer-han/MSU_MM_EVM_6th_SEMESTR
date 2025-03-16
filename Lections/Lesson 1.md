Заполнение матрицы с помощью MPI
p=3
n \* n
![[IMG_9437.jpg]]
Есть глобальная и локальные нумерации матриц для каждого процесса
Для получения блочного номера делим на $m$
Строки локальной матрицы расположены в глобальной по модулю p. 

Распишем преобразования номеров:
``` c++
// local to global
int l2g (
	int n,
	int m,
	int p,
	int k,
	int i_loc
) { 
	// Номер блочной строки локальный
	int i_loc_m = i_loc / m;
	// Номер блочной строки глобальной
	int i_glob_m = i_loc_m * p + k;
	return i_glob_m * m + i_loc % m;
}

// global to local
int g2l (
	int n,
	int m,
	int p,
	int k,
	int i_glob
) {
	// номер блочной строки глобальной
	int i_glob_m = i_glob / m;
	// номер блочной строки локальной
	int i_loc_m = i_glob_m / p;
	return i_loc_m * m + i_glob % m;
}

// максимальное число строк локальной матрицы
int get_max_rows(
	int n,
	int m,
	int p
) {
	int b = (n + m - 1) / m;
	return (b + p - 1) / p;
}

// число блочных строк на процесс
int get_rows(
	int n,
	int m,
	int p,
	int k
) {
	int b = (n + m - 1) / m;
	return b % p <= k ? b / p : b / p + 1;
}

// в каком процессе лежит строка?
int get_k(
	int n,
	int m,
	int p,
	int i_glob
) {
	int i_glob_m = i_glob / m;
	return i_glob_m % p;
}
```


Теперь реализуем инициализацию матрицы:
``` c++
//a_ij = f(s, n, i, j)
void initmatrix(
	double *a,
	int n,
	int m,
	int p,
	int k,
	double (*f)(int s, int n, int i, int j)
) {
	int i_loc, j_loc, i_glob, j_glob, rows;
	// сколько строк в процессе
	rows = get_rows(n, m, p, k);
	for (i_loc = 0, i_loc < rows; i_loc++) {
		i_glob = l2g(n, m, p, k, i_loc);
		for (j_loc = 0, j_loc < n, j_loc++) {
			j_glob = j_loc;
			a[i_loc *  n + j_loc] = (*f)(s, n, i_glob, j_glob);
		}
	}
}
```

Чтение матрицы. Один процесс читает в свой буффер, отправляет:
``` c++
int read_matrix(
	double *a,
	int n,
	int m,
	int p,
	int k,
	const char *name,
	double *buf, // буффер - блочная строка n * m
	MPI_Comm com
) {
	int main_k = 0; // кто читает файл
	FILE *fp = nullptr;
	int err = 0;
	if (k == main_k) {
		fp = fopen("r", name);
		if (fp == nullptr) err = 1;
	}
	MPI_Bcast(&err, 1, MPI_INT, main_k, com);
	if (err) return err; // во всех процессах
	
	memset(buf, 0, n * m * sizeof(double));
	
	// число блочных строк
	int b, max_b = (n + m - 1) / m;
	for (b = 0, b < max_b, b++) {
		// владелец строки
		int owner = b % p;
		int rows = b * m + m <= n ? m : n - b * m;
		// rows = min (m, n - b * m)
		
		// лок номер строки
		int b_loc = b / p;
		if (k == main_k) {
			err += read_array(fp, buf, n * rows);
			// ошибки обрабатываем потом
			if (owner == main_k) {
				// владалец - главный, копируем строку на место
				memcpy(a + b_loc * n * m, buf, n * rows);
			} else {
				// надо отправить процессу owner
				MPI_Send(
					buf,
					n * rows,
					MPI_DOUBLE,
					owner, 
					0 /*tag*/,
					com
				);
			}
		} else {
			if (owner == k) {
				MPI_Status &st;
				MPI_Recv(
					a + b_loc * n * m,
					n * rows,
					MPI_DOUBLE,
					main_k,
					0, //tag
					com,
					&st
				)
			}
		}
	}
	
	if (k == main_k) {
		fclose(fp);
		fp = nullptr;
	}
	MPI_Bcast(&err, 1, MPI_INT, main_k, com);
	if (err) return err;
	return 0;
}

int read_array(FILE *fp, double *a, int len)
{
	for (int i = 0, i < len, i++) {
		if (fscanf(fp, "%lf", a + i) != 1) return -2;
	}
	return 0;
}


// печать матрицы
void print_matrix(
	double * a,
	int n,
	int m,
	int p,
	int k,
	double *buf, // n * m блочная строка
	int max_prtint,
	MPI_Comm com
) {
	int main_k = 0; // только 0 в большинстве систем
	int b, max_b = (n + m - 1) / m;
	int printed_rows = 0;
	for (b = 0; b < max_b; b++) {
		int owner = b % p; // где эта строка
		int rows = min(m, n - b * m);
		int b_loc = b / p;
		if (k == main_k) {
			if (owner == main_k) {
				// печать массива, который есть локально
				printed_rows += print_array(
					a + b_loc * n * m,
					n,
					rows,
					printed_rows,
					max_print
				);
			}
			else {
				// главный должен получить от владельца
				MPI_Status st;
				MPI_Recv(
					buf,
					n*rows,
					MPI_DOUBLE,
					owner,
					0, //tag
					com,
					&st
				);
				printed_rows += print_array(
					buf,
					n
					rows,
					printed_rows,
					max_print
				)
			}
		} else {
			// остальные процессы
			if (k == owner) {
				MPI_Send(
					a + b_loc * n * m,
					n * rows,
					MPI_DOUBLE,
					main_k,
					0, //tag
					com
				)
			}
		}
	}
}

// печать прямоугольной матрицы с адресом a, длиной строки n, числом строк rows
// число напечатаных строк <= max_print
// возвращает число напечатаных строк
int print_array(
	double *a,
	int n,
	int rows,
	int printed_rows,
	int max_print
) {
	// число печатаемых столбцов
	int p_n = (n > max_print ? max_print : n);
	if (printed_rows >= max_print) return 0;
	
	// количество печатаемых строк
	int p_r = (printed_rows + rows <= max_print ? rows : max_print - rows);
	
	for (int i = 0; i < p_r; i++) { // по строкам
		for (int j = 0; j < p_n; j++) { // по столбцам
			printf(" %10.3e", a[i * n + j]);
		}
		printf("\n");
	}
	return p_r;
}
```


# Отчёт по практикуму. MPI версия
1) Описать разделение данных. 
	Свои данные.
	Локальная и глобальная нумерация.
2) Формулы в локальной нумерации.
3) Явно написать точки коммуникации.
	a. Количество точек ~4-6n
	b. Объём перессылаемых данных ~n\*n
	- Коллективный обмен: 1 блочную строку отправить всем
		Рассылает всем за логарифмическое количество данных, то есть:
		a. 1
		b. 1 строка (n * m)
	- Попарный обмен
		i-й отправляет i+1-му
		a. 1
		b. 1 строка
	- Попарный обмен: циклом отправил главный элемент всем
		a. p
		b. p * строку
	- Аналогично MPI_Reduce + MPI_Bcast - 2й обмен, а MPI_Allreduce - 1й
4) Сложность
	- Число операций S(n, m, p)
	- Число обменов C(n, m, p)
	- Объём обменов V(n, m, p)
	$S(n, m, p) = 1/p(c_1 * n^3 + c_2 * n^2*m + c_3 * n*m^2 + c_4 * m^3) + O(n^2 + nm + m^2)$
	$S(n, m, 1)$
	$S(n, n, p)$
	$S(n, 1, p)$