#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include <sys/sysinfo.h>
#include <sys/time.h>

#include "thread_work.h"
#include "solve.h"

double get_time()
{
  struct timeval time;
  gettimeofday (&time, 0);
  return time.tv_sec + time.tv_usec*1.e-6;
}

int process_args (thread_data *a)
{
  if (a->error_flag != 0)
    {
      if ((0 < a->p) && a[0].error_type == io_status::non_applicable)
        {
          printf("This method is not applicable\n");
          return 1;
        }
      return 1;
    }
  return 0;
}

void thread_data::set_func ()
{
  switch (func_id)
    {
      case 0: f = f_0; break;
      case 1: f = f_1; break;
      case 2: f = f_2; break;
      case 3: f = f_3; break;
      case 4: f = f_4; break;
      case 5: f = f_5; break;
      case 6: f = f_6; break;
      case 7: f = f_7; break;
    }
}

void reduce_sum (int p, int *a, int n)
{
  static pthread_mutex_t m = PTHREAD_MUTEX_INITIALIZER;
  static pthread_cond_t c_in = PTHREAD_COND_INITIALIZER;
  static pthread_cond_t c_out = PTHREAD_COND_INITIALIZER;
  static int t_in = 0;
  static int t_out = 0;
  static int *r = nullptr;
  int i;
  if (p <= 1)
    return;
  pthread_mutex_lock(&m);
  if (r == nullptr)
    r = a;
  else
    for (i = 0; i < n; i++)
      r[i] += a[i];
  t_in++;
  if (t_in >= p)
    {
      t_out = 0;
      pthread_cond_broadcast(&c_in);
    }
  else
    while (t_in < p)
      pthread_cond_wait (&c_in, &m);
  
  if (r != a)
    for (i = 0; i < n; i++)
      a[i] = r[i];
  t_out++;
  if (t_out >= p)
    {
      t_in = 0;
      r = nullptr;
      pthread_cond_broadcast(&c_out);
    }
  else
    while (t_out < p)
      pthread_cond_wait (&c_out, &m);
  pthread_mutex_unlock(&m);
}

void reduce_sum (int p, double *a, int n)
{
  static pthread_mutex_t m = PTHREAD_MUTEX_INITIALIZER;
  static pthread_cond_t c_in = PTHREAD_COND_INITIALIZER;
  static pthread_cond_t c_out = PTHREAD_COND_INITIALIZER;
  static int t_in = 0;
  static int t_out = 0;
  static double *r = nullptr;
  int i;
  if (p <= 1)
    return;
  pthread_mutex_lock(&m);
  if (r == nullptr)
    r = a;
  else
    for (i = 0; i < n; i++)
      r[i] += a[i];
  t_in++;
  if (t_in >= p)
    {
      t_out = 0;
      pthread_cond_broadcast(&c_in);
    }
  else
    while (t_in < p)
      pthread_cond_wait (&c_in, &m);
  
  if (r != a)
    for (i = 0; i < n; i++)
      a[i] = r[i];
  t_out++;
  if (t_out >= p)
    {
      t_in = 0;
      r = nullptr;
      pthread_cond_broadcast(&c_out);
    }
  else
    while (t_out < p)
      pthread_cond_wait (&c_out, &m);
  pthread_mutex_unlock(&m);
}

void reduce_max (int p, double *a, int n)
{
  static pthread_mutex_t m = PTHREAD_MUTEX_INITIALIZER;
  static pthread_cond_t c_in = PTHREAD_COND_INITIALIZER;
  static pthread_cond_t c_out = PTHREAD_COND_INITIALIZER;
  static int t_in = 0;
  static int t_out = 0;
  static double *r = nullptr;
  int i;
  if (p <= 1)
    return;
  pthread_mutex_lock(&m);
  if (r == nullptr)
    r = a;
  else
    for (i = 0; i < n; i++)
      if (r[i] < a[i])
        r[i] = a[i];
  t_in++;
  if (t_in >= p)
    {
      t_out = 0;
      pthread_cond_broadcast(&c_in);
    }
  else
    while (t_in < p)
      pthread_cond_wait (&c_in, &m);
  
  if (r != a)
    for (i = 0; i < n; i++)
      a[i] = r[i];
  t_out++;
  if (t_out >= p)
    {
      t_in = 0;
      r = nullptr;
      pthread_cond_broadcast(&c_out);
    }
  else
    while (t_out < p)
      pthread_cond_wait (&c_out, &m);
  pthread_mutex_unlock(&m);
}

void *thread_func (void* ptr)
{
  thread_data *pthread_arg = (thread_data *) ptr;
  
  /* Reading data from thread_data */
  int k = pthread_arg->pi;
  int p = pthread_arg->p;
  double eps = pthread_arg->eps;
  int max_it = pthread_arg->max_it;
  
  int *I = pthread_arg->I;
  double *A, *B, *x, *r, *u, *v;
  int nx, ny;
  int func_id;
  double (*f) (double, double);
  double max_abs_f;
  int disturbance;
  double hx, hy;
  
  double a = pthread_arg->a;
  double b = pthread_arg->b;
  double c = pthread_arg->c;
  double d = pthread_arg->d;
  
  /* Reading data from thread_data */
  
  int p_calc = p - 1;
  int k_calc = k - 1;
  int k_main_calc = 0;
  
  /* INIT MUTEX, CPU */
  pthread_mutex_t &mutex_gui_kernel = *pthread_arg->mutex_gui_kernel;
  pthread_cond_t &cond_gui_kernel = *pthread_arg->cond_gui_kernel;
  
  int &need_to_calculation = *pthread_arg->need_to_calculation;
  int &calculation_is_ready = *pthread_arg->calculation_is_ready;
  int &is_closing = *pthread_arg->is_closing;
  
  thread_data *main_arg = pthread_arg->pthread_args;
  pthread_t tid = pthread_arg->tid;
  
  cpu_set_t cpu;
  CPU_ZERO (&cpu);
  
  int n_cpus = get_nprocs();
  int cpu_id = n_cpus - 1 - (k_calc % n_cpus);
  CPU_SET(cpu_id, &cpu);
  pthread_setaffinity_np(tid, sizeof(cpu), &cpu);
  /* INIT MUTEX, CPU */
  
  int it = 0;
  double t;
  int max_step = 100;
  
  /* Mutex waiting */
  pthread_mutex_lock (&mutex_gui_kernel);
  //printf("[k = %d]: thr_func_setup: need_to_calc = %d\n", k, need_to_calculation);
  //printf("[k = %d]: thr_func_setup: is_closing   = %d\n", k, is_closing);
  //printf("[k = %d]: thr_func_setup: calc_is_ready = %d\n", k, calculation_is_ready);
  while (!is_closing && !need_to_calculation)
    {
      //printf("[k = %d]: thr_func_setup: is_waiting\n", k);
      pthread_cond_wait (&cond_gui_kernel, &mutex_gui_kernel);
    }
  pthread_mutex_unlock(&mutex_gui_kernel);
  /* Mutex waiting */
    
  if (k_calc == k_main_calc)
    {
      if (is_closing)
        main_arg->is_running = 0;
    }
  reduce_sum (p_calc);
  pthread_arg->is_running = main_arg->is_running;
      
  
  
  
  while (pthread_arg->is_running)
    {
      if (k_calc == k_main_calc) printf("[k = %d]: start calculation\n", k);
        
      /* Reading mutable data from thread_data */
      I = pthread_arg->I = main_arg->I;
      A = pthread_arg->A = main_arg->A;
      B = pthread_arg->B = main_arg->B;
      x = pthread_arg->x = main_arg->x;
      r = pthread_arg->r = main_arg->r;
      u = pthread_arg->u = main_arg->u;
      v = pthread_arg->v = main_arg->v;
      
      func_id = pthread_arg->func_id = main_arg->func_id;
      pthread_arg->set_func ();
      f = pthread_arg->f;
      max_abs_f = max_f_ab (a, b, c, d, func_id, f);
      nx = pthread_arg->nx = main_arg->nx;
      ny = pthread_arg->ny = main_arg->ny;
      disturbance = pthread_arg->disturbance = main_arg->disturbance;
      /* Reading mutable data from thread_data */
      
      //if (k_calc == k_main_calc) printf("[k = %d]: begin: func_id = %d\n", k, pthread_arg->func_id);
      
      // Вычисление структуры MSR матрицы и B. Заполнение
      hx = (b - a) / nx;
      hy = (d - c) / ny;
      
      int i1, i2;
      thread_rows((nx + 1) * (ny + 1), p_calc, k_calc, i1, i2);
      
      memset (I + i1, 0, (i2 - i1) * sizeof(int));
      memset (A + i1, 0, (i2 - i1) * sizeof(double));
      memset (B + i1, 0, (i2 - i1) * sizeof(double));
      memset (x + i1, 0, (i2 - i1) * sizeof(double));
      memset (r + i1, 0, (i2 - i1) * sizeof(double));
      memset (u + i1, 0, (i2 - i1) * sizeof(double));
      memset (v + i1, 0, (i2 - i1) * sizeof(double));
      
      reduce_sum (p_calc);
      
      if (k_calc == k_main_calc)
        {
          fill_I (nx, ny, hx, hy, I);
        }
        
      reduce_sum (p_calc);
      
      memset (I + I[i1], 0, (I[i2] - I[i1]) * sizeof(int));
      memset (A + I[i1], 0, (I[i2] - I[i1]) * sizeof(double));
      
      pthread_arg->error_flag = fill_IA (nx, ny, hx, hy, I, A, p_calc, k_calc);
      
      reduce_sum (p_calc, &pthread_arg->error_flag, 1);
      
      init_b (nx, ny, a, b, c, d, f, B, p_calc, k_calc, disturbance, max_abs_f);
      
      /*if (k == k_main_calc)
        {
          print_MSR_mat (nx, ny, I, A); // <- ####
          print_B (nx, ny, B, 42); // <- ####
        }*/
        
      /*if (k == 0)
        for (int i = 0; i < (nx + 1) * (ny + 1); i++)
          x[i] = 0;*/
      
      if (k_calc == k_main_calc)
        {
          t = get_time ();
        }
      it = minimal_errors_msr_matrix_full ((nx+1)*(ny+1), A, I, B, x // in-out 
                                          , r, u, v, eps, max_it, max_step, p_calc, k_calc);
      if (k_calc == k_main_calc)
        {
          pthread_arg->t1 = get_time () - t;
        }
      
      if (it < 0)
        {
          pthread_arg->error_type = io_status::non_applicable;
          pthread_arg->error_flag = 1;
          return nullptr;
        }
      
      
      /*if (debug && k_calc == k_main_calc)
        {
          //printf ("Result vector x:\n");
          //print_B (nx, ny, x, (nx+1)*(ny+1)); // <- ####
          int n = 6;
          double dx = (b-a)/n, dy = (d-c)/n;
          for (int j = n; j >= 0; j--)
            {
              for (int i = 0; i <= n; i++)
                printf(" [%9.2e]", p_f (nx, ny, a, b, c, d, x, a + i*dx, c + j*dy));
              printf("\n");
            }
        }*/
      
      if (k_calc == k_main_calc)
        {
          t = get_time ();
        }
      pthread_arg->r1 = residual_C_norm_error (nx, ny, a, b, c, d, x, f, p_calc, k_calc);
      pthread_arg->r2 = residual_L1_norm_error (nx, ny, a, b, c, d, x, f, p_calc, k_calc);
      pthread_arg->r3 = residual_C_norm_diff (nx, ny, a, b, c, d, x, f, p_calc, k_calc);
      pthread_arg->r4 = residual_L1_norm_diff (nx, ny, a, b, c, d, x, f, p_calc, k_calc);
      if (k_calc == k_main_calc)
        {
          pthread_arg->t2 = get_time () - t;
        }
      
      pthread_arg->it = it;
      
      /* Mutex sending */
      //pthread_mutex_lock (&mutex_gui_kernel);
      reduce_sum (p_calc);
      if (k_calc == k_main_calc) printf("[k = %d]: end calculation\n", k);
      if (k_calc == k_main_calc)
        {
          calculation_is_ready = 1;
          //if (!is_closing)
          need_to_calculation = 0;
        }
      //if (k_calc == k_main_calc) printf("[k = %d]: calc_is_ready = %d, need_to_calc = %d\n", k, calculation_is_ready, need_to_calculation);
      
      //pthread_cond_broadcast(&cond_gui_kernel);
      //pthread_mutex_unlock(&mutex_gui_kernel);
      /* Mutex sending */
      
      /* Mutex waiting */
      pthread_mutex_lock (&mutex_gui_kernel);
      //printf("[k = %d]: begin: need_to_calculation = %d\n", k, need_to_calculation);
      while (!is_closing && !need_to_calculation)
        {
          pthread_cond_wait (&cond_gui_kernel, &mutex_gui_kernel);
        }
      pthread_mutex_unlock (&mutex_gui_kernel);
      reduce_sum (p_calc);
      
      if (k_calc == k_main_calc)
        {
          if (is_closing)
            main_arg->is_running = 0;
        }
      reduce_sum (p_calc);
      pthread_arg->is_running = main_arg->is_running;
      printf("[k = %d]: pthread_arg->is_running = %d\n", k, pthread_arg->is_running);
      
      /*pthread_mutex_lock (&mutex_gui_kernel);
      if (k_calc == k_main_calc)
        {
          //else
          need_to_calculation = 0;
          calculation_is_ready = 1;
        }
      if (k_calc == k_main_calc)
      pthread_mutex_unlock(&mutex_gui_kernel);*/
      /* Mutex waiting */
    }
  
  pthread_arg->error_type = io_status::success;
  return nullptr;
}

void *single_func (void* ptr)
{
  thread_data *pthread_arg = (thread_data *) ptr;
  
  /* Reading data from thread_data */
  double eps = pthread_arg->eps;
  int max_it = pthread_arg->max_it;
  
  int *I = pthread_arg->I;
  double *A, *B, *x, *r, *u, *v;
  double a, b, c, d;
  int nx, ny;
  double (*f) (double, double);
  
  /* Reading data from thread_data */
  
  /* INIT MUTEX, CPU */
  
  int &need_to_calculation = *pthread_arg->need_to_calculation;
  int &calculation_is_ready = *pthread_arg->calculation_is_ready;
  
  thread_data *main_arg = pthread_arg->pthread_args;
  
  
  int p_calc = 1;
  int k_calc = 0;
  int k_main_calc = 0;
  int it = 0;
  double t;
  int max_step = 100, func_id;
  double hx, hy;
  int disturbance;
  double max_abs_f;
  
  a = pthread_arg->a;
  b = pthread_arg->b;
  c = pthread_arg->c;
  d = pthread_arg->d;
  
    
  /* Reading mutable data from thread_data */
  I = pthread_arg->I = main_arg->I;
  A = pthread_arg->A = main_arg->A;
  B = pthread_arg->B = main_arg->B;
  x = pthread_arg->x = main_arg->x;
  r = pthread_arg->r = main_arg->r;
  u = pthread_arg->u = main_arg->u;
  v = pthread_arg->v = main_arg->v;
  
  nx = pthread_arg->nx = main_arg->nx;
  ny = pthread_arg->ny = main_arg->ny;
  
  func_id = pthread_arg->func_id = main_arg->func_id;
  pthread_arg->set_func ();
  f = pthread_arg->f;
  max_abs_f = max_f_ab (a, b, c, d, func_id, f);
  disturbance = pthread_arg->disturbance = main_arg->disturbance;
  /* Reading mutable data from thread_data */
  
  //printf("[k = 0]: begin: func_id = %d\n", k, pthread_arg->func_id);
  
  // Вычисление структуры MSR матрицы и B. Заполнение
  hx = (b - a) / nx;
  hy = (d - c) / ny;
  
  int i1, i2;
  thread_rows((nx + 1) * (ny + 1), p_calc, k_calc, i1, i2);
  
  memset (I + i1, 0, (i2 - i1) * sizeof(int));
  memset (A + i1, 0, (i2 - i1) * sizeof(double));
  memset (B + i1, 0, (i2 - i1) * sizeof(double));
  memset (x + i1, 0, (i2 - i1) * sizeof(double));
  memset (r + i1, 0, (i2 - i1) * sizeof(double));
  memset (u + i1, 0, (i2 - i1) * sizeof(double));
  memset (v + i1, 0, (i2 - i1) * sizeof(double));
  
  reduce_sum (p_calc);
  
  if (k_calc == k_main_calc)
    {
      fill_I (nx, ny, hx, hy, I);
    }
    
  reduce_sum (p_calc);
  
  memset (I + I[i1], 0, (I[i2] - I[i1]) * sizeof(int));
  memset (A + I[i1], 0, (I[i2] - I[i1]) * sizeof(double));
  
  pthread_arg->error_flag = fill_IA (nx, ny, hx, hy, I, A, p_calc, k_calc);
  
  reduce_sum (p_calc, &pthread_arg->error_flag, 1);
  
  init_b (nx, ny, a, b, c, d, f, B, p_calc, k_calc, disturbance, max_abs_f);
  
  /*if (k == k_main_calc)
    {
      //print_MSR_mat (nx, ny, I, A); // <- ####
      //print_B (nx, ny, B, 42); // <- ####
    }*/
  // Окончание заполнения массива и вывода  
    
  // ###############
  /*if (k == 0)
    for (int i = 0; i < (nx + 1) * (ny + 1); i++)
      x[i] = 0;*/
  
  if (k_calc == k_main_calc)
    {
      t = get_time ();
    }
  it = minimal_errors_msr_matrix_full ((nx+1)*(ny+1), A, I, B, x // in-out 
                                       , r, u, v, eps, max_it, max_step, p_calc, k_calc);
  if (k_calc == k_main_calc)
    {
      pthread_arg->t1 = get_time () - t;
    }
  
  if (it < 0)
    {
      pthread_arg->error_type = io_status::non_applicable;
      pthread_arg->error_flag = 1;
      return nullptr;
    }
  
  
  /*if (debug && k_calc == k_main_calc)
    {
      //printf ("Result vector x:\n");
      //print_B (nx, ny, x, (nx+1)*(ny+1)); // <- ####
      int n = 6;
      double dx = (b-a)/n, dy = (d-c)/n;
      for (int j = n; j >= 0; j--)
        {
          for (int i = 0; i <= n; i++)
            printf(" [%9.2e]", p_f (nx, ny, a, b, c, d, x, a + i*dx, c + j*dy));
          printf("\n");
        }
    }*/
  
  t = get_time ();
  
  pthread_arg->r1 = residual_C_norm_error (nx, ny, a, b, c, d, x, f, p_calc, k_calc);
  pthread_arg->r2 = residual_L1_norm_error (nx, ny, a, b, c, d, x, f, p_calc, k_calc);
  pthread_arg->r3 = residual_C_norm_diff (nx, ny, a, b, c, d, x, f, p_calc, k_calc);
  pthread_arg->r4 = residual_L1_norm_diff (nx, ny, a, b, c, d, x, f, p_calc, k_calc);
  
  pthread_arg->t2 = get_time () - t;
  
  pthread_arg->it = it;
  
  //printf("[k = 0]: end calculation\n", k);                            , r, u, v, eps, max_it, max_step, omega, p_calc, k_calc);
  if (k_calc == k_main_calc)
    {
      calculation_is_ready = 1;
      need_to_calculation = 0;
    }
  printf("[k = %d]: calculation_is_ready = %d, need_to_calculation = %d\n", k_calc, calculation_is_ready, need_to_calculation);
  
  printf("[k = %d]: begin: need_to_calculation = %d\n", k_calc, need_to_calculation);
  
  printf("[k = %d]: start calculation\n", k_calc);
  
  pthread_arg->error_type = io_status::success;
  return nullptr;
}
