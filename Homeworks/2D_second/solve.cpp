#include <stdio.h>
#include <math.h>

#include "matrix_Gram.h"
#include "for_pthread.h"
#include "solve.h"

static double *results = nullptr;

void thread_rows (int n, int p, int k, int &i1, int &i2)
{
  i1 = n * k;
  i1 /= p; // i1 = n * k / p;
  i2 = n * (k + 1);
  i2 /= p; // i2 = n * (k + 1) / p;
}

double max_f_ab (double a, double b, double c, double d, int func_id, double (*f) (double, double))
{
  double p1 = 0, p2 = 0, p3 = 0, p4 = 0, max = 0;
  switch (func_id)
    {
      case 0:
        return 1.;
        break;
      case 1:
      case 2:
      case 3:
        return f (b, d);
        break;
      case 4:
      case 5:
        p1 = f (a, c); max = p1;
        p2 = f (b, c); max = (p2 > max ? p2 : max);
        p3 = f (a, d); max = (p3 > max ? p3 : max);
        p4 = f (b, d); max = (p4 > max ? p4 : max);
        return max;
        break;
      case 6:
        if (c * d <= 0)
          {
            p1 = f (a, 0); max = p1;
            p2 = f (b, 0); max = (p2 > max ? p2 : max);
            return max;
          }
        else
          {
            p1 = f (a, c); max = p1;
            p2 = f (b, c); max = (p2 > max ? p2 : max);
            p3 = f (a, d); max = (p3 > max ? p3 : max);
            p4 = f (b, d); max = (p4 > max ? p4 : max);
            return max;
          }
        break;
      case 7:
        if (a * b <= 0)
          {
            if (c * d <= 0)
              {
                return f (0, 0);
              }
            else
              {
                return (d > 0 ? f(0, c) : f(0, d));
              }
          }
        else
          {
            if (c * d <= 0)
              {
                return (a > 0 ? f(a, 0) : f(b, 0));
              }
            else
              {
                p1 = f (a, c); max = p1;
                p2 = f (b, c); max = (p2 > max ? p2 : max);
                p3 = f (a, d); max = (p3 > max ? p3 : max);
                p4 = f (b, d); max = (p4 > max ? p4 : max);
                return max;
              }
          }
        break;
    }
  return 1e308;
}

void matrix_mult_vector_msr (int n, double *A, int *I, double *x, double *y, int p, int k)
{
  int i, j, l, J;
  int i1, i2;
  double s;
  thread_rows (n, p, k, i1, i2);
  /*if (k == 0)
    {
      printf ("INIT vector x:\n");
      print_mas (n, x, 9); // <- ####
    }*/
  for (i = i1; i < i2; i++)
    {
      s = A[i] * x[i];
      l = I[i+1] - I[i];
      J = I[i];
      for (j = 0; j < l; j++)
        s += A[J + j] * x[I[J + j]];
      y[i] = s;
    }
  reduce_sum (p);
  /*if (k == 0)
    {
      printf ("RES vector y:\n");
      print_mas (n, y, 9); // <- ####
    }*/
}
  
int minimal_errors_msr_matrix (int n, double *A, int *I, double *b,
                               double *x, /* in-out */
                               double *r, double *u, double *v,
                               double eps, int max_it, double omega, int p, int k)
{
  double prec, b_norm2, tau, c1, c2;
  int it;
  b_norm2 = scalar_product (n, b, b, p, k);
  prec = b_norm2 * eps * eps;
  
  // r = Ax - b
  /*if (k == 0)
    {
      printf ("b_norm2 %10.3e\n", b_norm2);
      printf ("eps %10.3e\n", eps);
      printf ("b_norm2 * eps %10.3e\n", b_norm2 * eps);
      printf ("BEGIN: Vector x:\n");
      print_mas (n, x, 9); // <- ####
      printf ("BEGIN: Vector b:\n");
      print_mas (n, b, 9); // <- ####
    }*/
  matrix_mult_vector_msr (n, A, I, x, r, p, k); // r = Ax, 1 т. синх.
  /*if (k == 0)
    {
      printf ("BEGIN: Vector r:\n");
      print_mas (n, r, 9); // <- ####
    }*/
  mult_sub_vector(n, r, b, 1., p, k); // r -= 1 * b, 1 т. синх.
  /*if (k == 0)
    {
      printf ("BEGIN: Vector r:\n");
      print_mas (n, r, 9); // <- ####
    }*/
  
  for (it = 0; it < max_it; it++)
    {
      /*printf ("MIN_ERR: IT = %d ##################################\n", it);
      if (k == 0)
        {
          printf ("Vector r:\n");
          print_mas (n, r, 9); // <- ####
        }*/
      // v: Mv = r
      apply_preconditional_msr_matrix (n, A, I, v, r, omega, p, k); // 1 т. синх.
      /*if (k == 0)
        {
          printf ("Vector v:\n");
          print_mas (n, v, 9); // <- ####
        }*/
      
      // u = Av [= AM^(-1)r]
      matrix_mult_vector_msr (n, A, I, v, u, p, k); // 1 т. синх.
      /*if (k == 0)
        {
          printf ("Vector u:\n");
          print_mas (n, u, 9); // <- ####
        }*/
      
      // c_1 = (v, r)
      c1 = scalar_product(n, v, r, p, k);
      
      // c2 = (u, v)
      c2 = scalar_product(n, u, v, p, k);
      /*if (k == 0)
        {
          printf ("c1 %10.3e c2 %10.3e prec %10.3e\n", c1, c2, prec);
        }*/
      /*if (k == 0)
        {
          //printf ("Vector x:\n");
          print_mas (n, x, 9); // <- ####
        }*/
      if (c1 < prec || c2 < prec) {
        //if (k == 0)printf ("END:\n");
        break;
      }
      tau = c1/c2;
      /*if (k == 0)
        {
          printf ("tau %10.3e\n", tau);
        }*/
      // x -= tau * v;
      mult_sub_vector(n, x, v, tau, p, k); // 1 т. синх.
      
      // r -= tau * u
      mult_sub_vector(n, r, u, tau, p, k); // 1 т. синх.
    }
  if (it >= max_it)
    {
      return -1;
    }
  return it;
}

int minimal_errors_msr_matrix_full (int n, double *A, int *I, double *b,
                                   double *x, /* in-out */
                                   double *r, double *u, double *v, double eps, 
                                   int max_it, int max_step, double omega, 
                                   int p, int k)
{
  int step, ret, its = 0;
  for (step = 0; step < max_step; step++)
    {
      // много точек итераций
      ret = minimal_errors_msr_matrix (n, A, I, b, x, r, u, v, eps, max_it, omega, p, k);
      if (ret >= 0)
        {
          its += ret;
          break;
        }
      its += max_it;
    }
  if (step >= max_step)
    return -1;
  return its;
}


/*void apply_preconditional_msr_matrix (int n, double *A, int *, double *v,
                                        double *r, double, int p, int k)
{
  int i, i1, i2;
  thread_rows (n, p, k, i1, i2);
  // Якоби M = diag(A)
  for (i = i1; i < i2; i++) {
    v[i] = r[i] / A[i];
  }
  reduce_sum(p);
}*/

void apply_preconditional_msr_matrix (int n, double *A, int *I, double *v,
                                      double *r, double omega, int p, int k)
{
  (void)I;
  int i, i1, i2;
  int J, s_len, j;
  double res;
  thread_rows (n, p, k, i1, i2);
  // Метод верхней релаксации M^-1 = k(2-k)*(D+kR)^-1*D*(D+kL)^-1
  for (i = i1; i < i2; i++)
    {
      //printf ("###i = %d in [%d, %d)\n", i, i1, i2);
      res = r[i];
      J = I[i];
      s_len = I[i+1] - I[i];
      //printf ("   b[i] = %.1lf/24, J = %d, j_len = %d\n", res*24, J, s_len);
      //res -= A[i] * v[i];
      //printf ("A_ii = %lf, x[i] = %lf\n", A[i], v[i]);
      for (int s = 0; s < s_len; s++)
        {
          j = I[J + s];
          if (i1 <= j && j < i)
            {
              //printf ("   A[%d, %d] = %.1lf/24, x[%d] = %.1lf\n", i, j, A[J + s]*24, j, v[j]);
              res -= omega * A[J + s] * v[j];
            }
        }
      //printf ("   res = %.1lf/24 / A[i] = %.1lf/24 = %.3lf\n", res*24, A[i]*24, res / A[i]);
      v[i] = res / A[i];
    }
  //if(k ==0){printf ("1: "); print_mas (n, v, 9);}
  for (i = i1; i < i2; i++)
    { 
      //v[i]=r[i];
      v[i] *= omega * (2 - omega) * A[i];
    }
  //if(k ==0){printf ("2: "); print_mas (n, v, 9);}
  for (i = i2 - 1; i >= i1; i--)
    {
      res = v[i];
      J = I[i];
      s_len = I[i+1] - I[i];
      for (int s = 0; s < s_len; s++)
        {
          j = I[J + s];
          if (i <= j && j < i2)
            {
              res -= omega * A[J + s] * v[j];
            }
        }
      v[i] = res / A[i];
    }
  //if(k ==0){printf ("3: "); print_mas (n, v, 9);}
  reduce_sum(p);
}

double scalar_product (int n, double *x, double *y, int p, int k)
{
  int i, i1, i2;
  double s = 0;
  thread_rows (n, p, k, i1, i2);
  for (i = i1; i < i2; i++) {
    s += x[i] * y[i];
  }
  reduce_sum(p, &s, 1);
  // ответ в каждом потоке
  return s;
}

void mult_sub_vector (int n, double *x, double *y, double alpha, int p, int k)
{
  int i, i1, i2;
  thread_rows (n, p, k, i1, i2);
  for (i = i1; i < i2; i++) {
    x[i] -= y[i] * alpha;
  }
  reduce_sum(p);
}


int init_reduce_sum (int p)
{
  results = new double[p];
  if (results == nullptr)
    return -1;
  return 0;
}

int delete_reduce_sum ()
{
  //if (results == nullptr)
  //  return -1;
  delete[] results;
  return 0;
}

double reduce_sum_det (int p, int k, double s)
{
  results[k] = s;
  reduce_sum(p); // точка синхронизации
  double sum = 0;
  for (int l = 0; l < p; l++) {
    sum += results[l];
  }
  return sum;
}



// ###############################################################
// ###############################################################
// ###############################################################


int IA_ij (int n_x, int n_y, double hx, double hy, int i, int j, int is, int js, int s, int *I, double *A)
{
  (void) I;
  double Sq = hx * hy;
  int l, ls;
  ij2l (n_x, n_y, i, j, l);
  ij2l (n_x, n_y, is, js, ls);
  
  if (I != nullptr)
    I[s] = ls;
  if (i > 0 && i < n_x && j > 0 && j < n_y) // #
    {
      if (l == ls)
        A[s] = 6 * 1 / 12. * Sq;
      else
        A[s] = 2 * 1 / 24. * Sq;
    }
  if (i > 0 && i < n_x && j == 0) // _
    {
      if (l == ls)
        A[s] = 3 * 1 / 12. * Sq;
      else
      if (s == 0 || s == 1)
        A[s] = 1 * 1 / 24. * Sq;
      else
      if (s == 2 || s == 3)
        A[s] = 2 * 1 / 24. * Sq;
      else
        return -1;
    }
  if (i > 0 && i < n_x && j == n_y) // ^
    {
      if (l == ls)
        A[s] = 3 * 1 / 12. * Sq;
      else
      if (s == 0 || s == 3)
        A[s] = 1 * 1 / 24. * Sq;
      else
      if (s == 1 || s == 2)
        A[s] = 2 * 1 / 24. * Sq;
      else
        return -1;
    }
  if (i == 0 && j > 0 && j < n_y) // <
    {
      if (l == ls)
        A[s] = 3 * 1 / 12. * Sq;
      else
      if (s == 0 || s == 3)
        A[s] = 2 * 1 / 24. * Sq;
      else
      if (s == 1 || s == 2)
        A[s] = 1 * 1 / 24. * Sq;
      else
        return -1;
    }
  if (i == n_x && j > 0 && j < n_y) // >
    {
      if (l == ls)
        A[s] = 3 * 1 / 12. * Sq;
      else
      if (s == 0 || s == 3)
        A[s] = 1 * 1 / 24. * Sq;
      else
      if (s == 1 || s == 2)
        A[s] = 2 * 1 / 24. * Sq;
      else
        return -1;
    }
  if (i == 0 && j == 0) // <_
    {
      if (l == ls)
        A[s] = 2 * 1 / 12. * Sq;
      else
      if (s == 0 || s == 1)
        A[s] = 1 * 1 / 24. * Sq;
      else
      if (s == 2)
        A[s] = 2 * 1 / 24. * Sq;
      else
        return -1;
    }
  if (i == n_x && j == n_y) // >^
    {
      if (l == ls)
        A[s] = 2 * 1 / 12. * Sq;
      else
      if (s == 0 || s == 2)
        A[s] = 1 * 1 / 24. * Sq;
      else
      if (s == 1)
        A[s] = 2 * 1 / 24. * Sq;
      else
        return -1;
    }
  if ((i == n_x && j == 0) || (i == 0 && j == n_y)) // <^ or >_
    {
      if (l == ls)
        A[s] = 1 * 1 / 12. * Sq;
      else
      if (s == 0 || s == 1)
        A[s] = 1 * 1 / 24. * Sq;
      else
        return -1;
    }
  return 0;
}



// ###############################################################
// ###############################################################
// ###############################################################


int get_len_msr (int n_x, int n_y)
{
  return 6 * (n_x - 1) * (n_y - 1) + 4 * (2 * (n_x - 1) 
                                        + 2 * (n_y - 1)) + 3 * 2 + 2 * 2;
}

void ij2l (int n_x, int /*n_y*/, int i, int j, int &l)
{
  l = i + j * (n_x + 1);
}

void l2ij (int n_x, int /*n_y*/, int &i, int &j, int l)
{
  j = l / (n_x + 1);
  i = l - j * (n_x + 1);
}

#define F(IS, JS, S) \
        IA_ij (nx, ny, hx, hy, i, j, (IS), (JS), (S), I, A)

int get_off_diag (int nx, int ny, double hx, double hy, int i, int j,
                  int *I, double *A)
{
  int s = 0;
  if (i < nx)           { if (A != nullptr) F (i + 1, j,     s); s++; }
  if           (j > 0)  { if (A != nullptr) F (i,     j - 1, s); s++; }
  if (i > 0  && j > 0)  { if (A != nullptr) F (i - 1, j - 1, s); s++; }
  if (i > 0)            { if (A != nullptr) F (i - 1, j,     s); s++; }
  if           (j < ny) { if (A != nullptr) F (i,     j + 1, s); s++; }
  if (i < nx && j < ny) { if (A != nullptr) F (i + 1, j + 1, s); s++; }
  return s; // количество диагональных элементов
}

//возвращает количество внедиагональных точек (i, j)
int get_len_msr_off_diag (int nx, int ny)
{
  double hx = 0, hy = 0; int i, j, res = 0;
  for (i = 0; i <= nx; i++)
    for (j = 0; j <= ny; j++)
      res += get_off_diag (nx, ny, hx, hy, i, j);
  return res;
}

void get_diag (int nx, int ny, double hx, double hy, int i, int j, int* /*I*/, double *A)
{
  IA_ij (nx, ny, hx, hy, i, j, i, j, 0, nullptr, A);
}

void fill_I (int nx, int ny, double hx, double hy, int *I)
{
  int i, j, l, N = (nx + 1) * (ny + 1);
  int r = N + 1;
  for (l = 0; l < N; l++)
    {
      l2ij (nx, ny, i, j, l);
      int s = get_off_diag (nx, ny, hx, hy, i, j);
      I[l] = r;
      r += s;
    }
  I[l] = r;
}

int fill_IA (int nx, int ny, double hx, double hy, int *I, double *A, int p, int k)
{
  int i, j, l, l1, l2, N = (nx + 1) * (ny + 1), r, s, t;
  int err = 0, len = 0;
  thread_rows (N, p, k, l1, l2);
  for (l = l1; l < l2; l++)
    {
      r = I[l];
      s = I[l+1] - I[l];
      l2ij (nx, ny, i, j, l);
      get_diag (nx, ny, hx, hy, i, j, I,  A + l);
      t = get_off_diag (nx, ny, hx, hy, i, j, I + r, A + r);
      if (t != s)
        {
          err = 1; break;
        }
      len += s;
    }
  reduce_sum (p, &err, 1);
  if (err != 0)
    return -1;
  reduce_sum (p, &len, 1);
  if (I[N] != (N + 1) + len)
    return -2;
  return 0;
}

void print_MSR (int nx, int ny, int *I, double *A, int pr)
{
  int m = get_len_msr(nx, ny) + 1 + (nx + 1) * (ny + 1);
  m = (m < pr ? m : pr);
  for (int i = 0; i < m; i++)
    printf (" %4d", I[i]);
  printf ("\n");
  for (int i = 0; i < m; i++)
    printf (" %4.1lf", A[i] * 24);
  printf ("\n");
}

void print_MSR_mat (int nx, int ny, int *I, double *A)
{
  int r, s;
  for (int l = 0; l < (nx + 1) * (ny + 1); l++)
    {
      r = I[l];
      s = I[l+1] - r;
      printf (" {%2d - %4.1lf}", l, A[l]*24);
      for (int i_s = 0; i_s < s; i_s++)
        {
          printf (" {%2d - %4.1lf}", I[r + i_s], A[r + i_s]*24);
        }
      printf ("\n");
    }
}

void print_mas (int n, double *x, int pr)
{
  int m = (n < pr ? n : pr);
  for (int i = 0; i < m; i++)
    printf (" %10.3e", x[i]);
  printf ("\n");
}



// ###############################################################
// ###############################################################
// ###############################################################


double p_f (int nx, int ny, double a, double b, double c, double d, double *x_mas, double x, double y)
{
  double hx = (b - a) / nx;
  double hy = (d - c) / ny;
  
  //int i = fmod(x - a, hx);
  //int j = fmod(y - c, hy);
  int i = floor((x - a)/hx);
  int j = floor((y - c)/hy);
  if (i < 0) i = 0;
  if (i >= nx) i = nx - 1;
  if (j < 0) j = 0;
  if (j >= ny) j = ny - 1;
  //printf ("\n%10.3e %10.3e %d %d\n", x, y, i, j);
  //printf ("%10.3e %10.3e %10.3e %10.3e\n", x, y, x - i * hx - a, y - j * hy - c);
  if (x - i * hx - a > y - j * hy - c) // _> triangle
    {
      int l1, l2, l3;
      ij2l (nx, ny, i, j, l1);
      ij2l (nx, ny, i+1, j+1, l2);
      ij2l (nx, ny, i+1, j, l3);
      return x_mas[l1] * ((i+1) * hx + a - x)/hx + x_mas[l2] * (y - j * hy - c)/hy + x_mas[l3] * ((x - i * hx - a)/hx - (y - j * hy - c)/hy);
    }
  else // <^ triangle
    {
      int l1, l2, l3;
      ij2l (nx, ny, i, j, l1);
      ij2l (nx, ny, i+1, j+1, l2);
      ij2l (nx, ny, i, j+1, l3);
      //printf ("### %10.3e %10.3e\n", x_mas[l3], x_mas[l2]);
      //printf ("### %10.3e\n", x_mas[l1]);
      //printf ("@@@ %10.3e %10.3e\n", ((y - j * hy - c)/hy - (x - i * hx - a)/hx), (x - i * hx - a)/hx);
      //printf ("@@@ %10.3e\n", (c + (j+1) * hy - y)/hy);
      return x_mas[l1] * (c + (j+1) * hy - y)/hy + x_mas[l2] * (x - i * hx - a)/hx + x_mas[l3] * ((y - j * hy - c)/hy - (x - i * hx - a)/hx);
    }
}

#define P_F(X, Y) \
        p_f (nx, ny, a, b, c, d, x_mas, (X), (Y))

double residual_C_norm_error (int nx, int ny, double a, double b, double c, double d, double *x_mas, double (*f) (double x, double y), int p, int k)
{
  int l1, l2;
  double x, y;
  double hx = (b - a) / nx;
  double hy = (d - c) / ny;
  thread_rows (2 * nx * ny, p, k, l1, l2);
  int i, i2, j;
  double diff = -1;
  
  double max = -1;
  for (int l = l1; l < l2; l++)
    {
      j = l / (2 * nx);
      i2 = l - j * (2 * nx);
      i = i2 / 2;
      if (i2 - i * 2 == 0) // _> triangle
        {
          x = a + (3 * i + 1) * hx / 3;
          y = c + (3 * j + 2) * hy / 3;
        }
      else
        {
          x = a + (3 * i + 2) * hx / 3;
          y = c + (3 * j + 1) * hy / 3;
        }
      diff = fabs (f(x, y) - P_F(x, y));
      if (max < diff)
        max = diff;
    }
  reduce_max (p, &max, 1);
  return max;
}

double residual_L1_norm_error (int nx, int ny, double a, double b, double c, double d, double *x_mas, double (*f) (double x, double y), int p, int k)
{
  int l1, l2;
  double x, y;
  double hx = (b - a) / nx;
  double hy = (d - c) / ny;
  thread_rows (2 * nx * ny, p, k, l1, l2);
  int i, i2, j;
  double diff = -1;
  
  double sum = 0;
  for (int l = l1; l < l2; l++)
    {
      j = l / (2 * nx);
      i2 = l - j * (2 * nx);
      i = i2 / 2;
      if (i2 - i * 2 == 0) // _> triangle
        {
          x = a + (3 * i + 1) * hx / 3;
          y = c + (3 * j + 2) * hy / 3;
        }
      else
        {
          x = a + (3 * i + 2) * hx / 3;
          y = c + (3 * j + 1) * hy / 3;
        }
      diff = fabs (f(x, y) - P_F(x, y));
      sum += diff * hx * hy / 2.;
    }
  reduce_sum (p, &sum, 1);
  return sum;
}

double residual_C_norm_diff (int nx, int ny, double a, double b, double c, double d, double *x_mas, double (*f) (double x, double y), int p, int k)
{
  double hx = (b - a) / nx;
  double hy = (d - c) / ny;
  double x, y;
  int j1, j2;
  double diff = -1;
  thread_rows (ny + 1, p, k, j1, j2);

  double max = -1;
  for (int j = j1; j < j2; j++)
    {
      for (int i = 0; i <= nx; i++)
        {
          x = a + i * hx;
          y = c + j * hy;
          diff = fabs (f(x, y) - P_F(x, y));
        }
      if (max < diff)
        max = diff;
    }
  reduce_max (p, &max, 1);
  return max;
}

double residual_L1_norm_diff (int nx, int ny, double a, double b, double c, double d, double *x_mas, double (*f) (double x, double y), int p, int k)
{
  double hx = (b - a) / nx;
  double hy = (d - c) / ny;
  double x, y;
  int j1, j2;
  thread_rows (ny + 1, p, k, j1, j2);
  double diff = -1;
  
  double sum = 0;
  for (int j = j1; j < j2; j++)
    {
      for (int i = 0; i <= nx; i++)
        {
          x = a + i * hx;
          y = c + j * hy;
          diff = fabs (f(x, y) - P_F(x, y));
          sum += diff * hx * hy;
        }
    }
  reduce_sum (p, &sum, 1);
  return sum;
}
  
     
