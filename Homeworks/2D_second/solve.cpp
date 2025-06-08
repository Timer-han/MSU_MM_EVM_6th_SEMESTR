#include <stdio.h>
#include <math.h>

#include "msr_matrix.h"
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

