#ifndef FOR_PTHREAD_H
#define FOR_PTHREAD_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include <sys/sysinfo.h>
#include <sys/time.h>


enum class io_status
{
  undefine,
  error_open,
  error_read,
  error_buf,
  error_mem,
  degenerate_data,
  error_incorrect_parametres,
  non_applicable,
  success,
};

class thread_data
{
public:
  pthread_t tid;
        
  double a = 0;
  double b = 0;
  double c = 1;
  double d = 1;
  int nx = 0;
  int ny = 0;
  int func_id = 0;
  int pi = 0;
  int p = 1;
  double eps = 1e-14;
  int max_it = 100000;
  
  int *I = 0;
  double *A = 0;
  double *B = 0;
  double *x = 0;
  double *r = 0;
  double *u = 0;
  double *v = 0;
  
  double (*f) (double, double);
  
  int it = 0;
  double t1 = 0;
  double t2 = 0;
  double r1 = 0, r2 = 0, r3 = 0, r4 = 0;
  int *mas_k= 0;
  int count = 0;
  io_status error_type = io_status::undefine;
  int error_flag = 0;
  
  int disturbance = 0;
  
  thread_data *pthread_args;
  void set_func ();
  int *need_to_calculation;
  int *calculation_is_ready;
  int *is_closing;
  
  pthread_mutex_t *mutex_gui_kernel;
  pthread_cond_t *cond_gui_kernel;
  
  int is_running = 1;
};

double get_time();

int process_args (thread_data *a);

void reduce_sum (int p, int *a = nullptr, int n = 0);
void reduce_sum (int p, double *a, int n);
void reduce_max (int p, double *a, int n);

void *thread_func (void* ptr);
void *single_func (void* ptr);

#endif // FOR_PTHREAD_H

