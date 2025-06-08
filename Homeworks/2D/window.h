#ifndef WINDOW_H
#define WINDOW_H

#include <QtWidgets/QtWidgets>
#include <pthread.h>
#include <sys/sysinfo.h>
#include "for_pthread.h"

class Window : public QWidget
{
  Q_OBJECT

private:
  char *prog_name = nullptr;

  double a, b, c, d;
  double eps = 1e-15;
  double omega;
  
  size_t nx = 1, ny = 1;
  size_t mx = 1, my = 1;
  int p = 2;
  int max_it = 1, err = 0;
  int it = 0;
  double t1 = 0, t2 = 0;
  double r1 = 0, r2 = 0, r3 = 0, r4 = 0;
  
  int func_id = 0;
  int zoom = 0;
  int disturbance = 0;
  double (*f) (double, double);
  int view_id = 0;
  int print_graph_func = 1;
  int print_graph_aprox = 0;
  int print_graph_diff = 0;
  
  int func_id_print = 0;
  size_t nx_print = 1;
  size_t ny_print = 1;
  int disturbance_print = 0;
  double (*f_print) (double, double);
  
  pthread_mutex_t mutex_gui_kernel = PTHREAD_MUTEX_INITIALIZER;
  pthread_cond_t cond_gui_kernel = PTHREAD_COND_INITIALIZER;
  
  Args *pthread_args = nullptr;
  int *I_calc = nullptr;
  double *A_calc = nullptr, *B_calc = nullptr, *x_calc = nullptr;
  double *r_calc = nullptr, *u_calc = nullptr, *v_calc = nullptr;
  
  int *I_print = nullptr;
  double *A_print = nullptr, *B_print = nullptr, *x_print = nullptr;
  double *r_print = nullptr, *u_print = nullptr, *v_print = nullptr;
  
  int need_to_calculate = 0;
  int calculation_is_ready = 1;
  int is_closing = 0;
  
  
  const char *f_name;
  //size_t n;
  //size_t n_max;
  size_t m_max = 10'000;
  //size_t n_mas = 40'000'000;
  
  double *x_mas = nullptr;
  
  
  int is_painting = 0;
  QTimer *timer;
  
  int status = 0;

public:
  Window (QWidget *parent, int argc, char *argv[]);
  ~Window ();

  QSize minimumSizeHint () const;
  QSize sizeHint () const;
  int get_status ();
  int parse_command_line (int argc, char *argv[]);
  void set_func ();
  void set_func_print ();
  void set_view ();
  void timerIsTimeout ();
  QPointF l2g (double x_loc, double y_loc);
  QColor l2g_colour (double z_loc, double z_min, double z_max);
  
  int enough_first_memory ();
  int enough_memory_for_calc ();
  void set_memory_from_calc_to_print ();
  void delete_memory_for_print ();
  int delete_last_memory ();
  
  int create_threads ();
  int close_threads ();
  
  void recalculate ();
  void requst_to_calculate ();
  
public slots:
  void change_func ();
  void change_view ();
  void change_zoom_up ();
  void change_zoom_down ();
  void change_nx_ny_down ();
  void change_nx_ny_up ();
  void change_mx_my_down ();
  void change_mx_my_up ();
  void change_disturbance_down ();
  void change_disturbance_up ();

protected:
  void paintEvent (QPaintEvent *event);
};

#endif // WINDOW_H
