#include <QPainter>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "window.h"
#include "for_pthread.h"
#include "msr_matrix.h"
#include "f.h"
#include "solve.h"
#include "init_b.h"
#include "residual.h"

#define OMEGA 1

#define MONITOR_SIZE_W 1920
#define DEFAULT_A -10
#define DEFAULT_B 10
#define DEFAULT_C -10
#define DEFAULT_D 10
#define DEFAULT_N 10

#define L2G(X,Y) (l2g ((X), (Y)))
#define L2GCOLOUR(Z) (l2g_colour ((Z), min_z, max_z))

Window::Window (QWidget *parent, int argc, char *argv[])
  : QWidget (parent)
{
  //printf("[k = 0]: WINDOW_CONSTRUCTOR\n");
  timer = new QTimer(this);
  connect (timer, &QTimer::timeout, this, &Window::timerIsTimeout);
  a = DEFAULT_A;
  b = DEFAULT_B;
  nx = DEFAULT_N;
  ny = DEFAULT_N;
  mx = DEFAULT_N;
  my = DEFAULT_N;

  func_id = 0;
  
  if (parse_command_line (argc, argv))
    {
      //printf("Usage : %s a b c d nx ny mx my k eps max_it p [omega]\n", argv[0]);
      QMessageBox::warning (0, "Wrong input arguments!", 
                            "Usage ./<a.out> a b c d nx ny mx my func_id eps max_it p [omega]");
      status = -1;
      return;
    }
  
  if (enough_first_memory ())
    {
      QMessageBox::warning (0, "Error", "Not enough memory!");
      status = -2;
      return;
    }
  
  if (create_threads ())
    {
      status = -3;
      return;
    }
  
  set_func ();
  set_view ();
  
  requst_to_calculate ();
  //printf("[k = 0]: WINDOW_CONSTRUCTOR end\n");
}


Window::~Window ()
{
  //printf("[k = 0]: WINDOW_DESTRUCTOR\n");

  if (!calculation_is_ready)
    {
      QMessageBox::information (0, "Message", "Application is closing.");
    }
  //while (need_to_calculate);
  //printf("@@@111\n");
  
  pthread_mutex_lock (&mutex_gui_kernel);
  is_closing = 1;
  pthread_cond_broadcast (&cond_gui_kernel);
  pthread_mutex_unlock (&mutex_gui_kernel);
  
  close_threads ();
  printf("@@@\n");
  delete_memory_for_print ();
  set_memory_from_calc_to_print ();
  delete_memory_for_print ();
  //printf("@@@\n");
  delete_last_memory();
  //printf("[k = 0]: WINDOW_DESTRUCTOR end\n");
}

QSize Window::minimumSizeHint () const
{
  return QSize (100, 100);
}

QSize Window::sizeHint () const
{
  return QSize (1000, 1000);
}

int Window::get_status ()
{
  return status;
}

int Window::parse_command_line (int argc, char *argv[])
{
  /* Чтение входных данных */
  if (!((argc == 13 || argc == 14)
        && sscanf(argv[1], "%lf", &a) == 1
        && sscanf(argv[2], "%lf", &b) == 1 // && a <= b
        && sscanf(argv[3], "%lf", &c) == 1
        && sscanf(argv[4], "%lf", &d) == 1 // && c <= d
        && b - a > 1.e-6
        && d - c > 1.e-6
        && sscanf(argv[5], "%zu", &nx) == 1 && (nx >= 1)
        && sscanf(argv[6], "%zu", &ny) == 1 && (ny >= 1)
        && sscanf(argv[7], "%zu", &mx) == 1 && (mx >= 1)
        && sscanf(argv[8], "%zu", &my) == 1 && (my >= 1)
        && sscanf(argv[9], "%d", &func_id) == 1 && (0 <= func_id) && (func_id <= 7)
        && sscanf(argv[10], "%lf", &eps) == 1 && (eps >= 0)
        && sscanf(argv[11], "%d", &max_it) == 1 && (max_it >= 0)
        && sscanf(argv[12], "%d", &p) == 1 && (p >= 1)
       ) )
    return -1;
  
  omega = OMEGA;
  
  func_id_print = func_id;
  nx_print = nx;
  ny_print = ny;
  disturbance_print = disturbance;
  set_func_print ();
  
  if ((argc == 14)
      && !(sscanf(argv[13], "%lf", &omega) == 1 && (0 <= omega) && (omega <= 2) ))
    {
      return -1;
    }
  prog_name = argv[0];
  
  /*if (n > n_mas)
    n = n_mas;
  if (n < 2)
    n = 2;*/
    
  // Окончание чтения входных данных
  return 0;
}

// ##############################################
// ################### Memory ###################

int Window::enough_first_memory ()
{
  // Выделение общей памяти
  if (   (pthread_args = new Args[p]) == nullptr
      || init_reduce_sum (p) != 0
     )
    {
      if (pthread_args != nullptr) delete[] pthread_args;
      delete_reduce_sum ();
      return -2;
    }
  return 0;
}

int Window::enough_memory_for_calc ()
{
  //printf("[k = 0]: enough_memory_for_calc\n");
  // Выделение общей памяти
  if (   (I_calc = new int[(nx + 1) * (ny + 1) + 1 + get_len_msr (nx, ny)]) == nullptr
      || (A_calc = new double[(nx + 1) * (ny + 1) + 1 + get_len_msr (nx, ny)]) == nullptr
      || (B_calc = new double[(nx + 1) * (ny + 1)]) == nullptr
      || (x_calc = new double[(nx + 1) * (ny + 1)]) == nullptr
      || (r_calc = new double[(nx + 1) * (ny + 1)]) == nullptr
      || (u_calc = new double[(nx + 1) * (ny + 1)]) == nullptr
      || (v_calc = new double[(nx + 1) * (ny + 1)]) == nullptr
     )
    {
      if (I_calc != nullptr) delete[] I_calc;
      if (A_calc != nullptr) delete[] A_calc;
      if (B_calc != nullptr) delete[] B_calc;
      if (x_calc != nullptr) delete[] x_calc;
      if (r_calc != nullptr) delete[] r_calc;
      if (u_calc != nullptr) delete[] u_calc;
      if (v_calc != nullptr) delete[] v_calc;
      return -2;
    }
  pthread_args[0].I = I_calc;
  pthread_args[0].A = A_calc;
  pthread_args[0].B = B_calc;
  pthread_args[0].x = x_calc;
  pthread_args[0].r = r_calc;
  pthread_args[0].u = u_calc;
  pthread_args[0].v = v_calc;
  //printf("[k = 0]: enough_memory_for_calc end\n");
  return 0;
}

void Window::set_memory_from_calc_to_print ()
{
  //printf("[k = 0]: set_memory_from_calc_to_print\n");
  I_print = I_calc; I_calc = nullptr;
  A_print = A_calc; A_calc = nullptr;
  B_print = B_calc; B_calc = nullptr;
  x_print = x_calc; x_calc = nullptr;
  r_print = r_calc; r_calc = nullptr;
  u_print = u_calc; u_calc = nullptr;
  v_print = v_calc; v_calc = nullptr;
  //printf("[k = 0]: set_memory_from_calc_to_print end\n");
}

void Window::delete_memory_for_print ()
{
  //printf("[k = 0]: delete_memory_for_print\n");
  /*if (I_print != nullptr) delete[] I_print;
  if (A_print != nullptr) delete[] A_print;
  if (B_print != nullptr) delete[] B_print;
  if (x_print != nullptr) delete[] x_print;
  if (r_print != nullptr) delete[] r_print;
  if (u_print != nullptr) delete[] u_print;
  if (v_print != nullptr) delete[] v_print;*/
  if (I_print != nullptr) { delete[] I_print; I_print = nullptr; }
  if (A_print != nullptr) { delete[] A_print; A_print = nullptr; }
  if (B_print != nullptr) { delete[] B_print; B_print = nullptr; }
  if (x_print != nullptr) { delete[] x_print; x_print = nullptr; }
  if (r_print != nullptr) { delete[] r_print; r_print = nullptr; }
  if (u_print != nullptr) { delete[] u_print; u_print = nullptr; }
  if (v_print != nullptr) { delete[] v_print; v_print = nullptr; }
  //printf("[k = 0]: delete_memory_for_print end\n");
}

int Window::delete_last_memory ()
{
  if (pthread_args != nullptr) 
    delete[] pthread_args;
  delete_reduce_sum ();
  return 0;
}

// ################### Memory ###################
// ##############################################

// ###############################################
// ################### Threads ###################

int Window::create_threads ()
{
  char mes[40] = {};
  for (int k = 0; k < p; k++)
    {
      pthread_args[k].a = a;
      pthread_args[k].b = b;
      pthread_args[k].c = c;
      pthread_args[k].d = d;
      pthread_args[k].nx = nx;
      pthread_args[k].ny = ny;
      pthread_args[k].func_id = func_id;
      pthread_args[k].k = k;
      pthread_args[k].p = p;
      pthread_args[k].eps = eps;
      pthread_args[k].max_it = max_it;
      pthread_args[k].omega = omega;
      
      pthread_args[k].set_func();
      pthread_args[k].mutex_gui_kernel = &mutex_gui_kernel;
      pthread_args[k].cond_gui_kernel = &cond_gui_kernel;
      
      pthread_args[k].need_to_calculation = &need_to_calculate;
      pthread_args[k].calculation_is_ready = &calculation_is_ready;
      pthread_args[k].is_closing = &is_closing;
      
      pthread_args[k].pthread_args = pthread_args;
    }
  for (int k = 1; k < p; k++)
    {
      if (pthread_create(&pthread_args[k].tid, nullptr, thread_func, pthread_args + k))
        {
          sprintf(mes, "Cannot create thread %d\n", k);
          QMessageBox::warning (0, "Error", mes);
          abort ();
        }
    }
  pthread_args[0].tid = pthread_self();

  return 0;
}

int Window::close_threads ()
{
  for (int k = 1; k < p; k++)
    pthread_join (pthread_args[k].tid, nullptr);
  process_args(pthread_args);
  if (p > 1 && pthread_args[1].error_type != io_status::success)
    {
      QMessageBox::warning (0, "Error", 
                            "Runtime error");
      return -5;
    }
  return 0;
}

// ################### Threads ###################
// ###############################################

void Window::requst_to_calculate ()
{
  if (p > 1)
    {
      if (enough_memory_for_calc ())
        {
          QMessageBox::warning (0, "Error", "Not enough memory!");
        }
      else
        {
          pthread_mutex_lock(&mutex_gui_kernel);
          calculation_is_ready = 0;
          need_to_calculate = 1;
          pthread_cond_broadcast(&cond_gui_kernel);
          pthread_mutex_unlock(&mutex_gui_kernel);
          timer->start (50);
          
          /*pthread_mutex_lock (&mutex_gui_kernel);
          need_to_calculate = 1;
          pthread_cond_broadcast (&cond_gui_kernel);
          pthread_mutex_unlock (&mutex_gui_kernel);
          
          pthread_mutex_lock (&mutex_gui_kernel);
          while (!have_calculated)
          {
              pthread_cond_wait (&cond_gui_kernel, &mutex_gui_kernel);
          }
          pthread_mutex_unlock (&mutex_gui_kernel);*/
        }
    }
  else
    {
      if (enough_memory_for_calc ())
        {
          QMessageBox::warning (0, "Error", "Not enough memory!");
        }
      else
        {
          single_func (pthread_args);
          x_mas = pthread_args[0].x;
          delete_memory_for_print ();
          set_memory_from_calc_to_print ();
          update ();
        }
    }
}

void Window::set_func ()
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

void Window::set_func_print ()
{
  switch (func_id_print)
    {
      case 0: f_print = f_0; f_name = "k = 0; f (x) = 1"; break;
      case 1: f_print = f_1; f_name = "k = 1; f (x) = x"; break;
      case 2: f_print = f_2; f_name = "k = 2; f (x) = y"; break;
      case 3: f_print = f_3; f_name = "k = 3; f (x) = x + y"; break;
      case 4: f_print = f_4; f_name = "k = 4; f (x) = sqrt(x^2 + y^2)"; break;
      case 5: f_print = f_5; f_name = "k = 5; f (x) = x^2 + y^2"; break;
      case 6: f_print = f_6; f_name = "k = 6; f (x) = exp (x^2 - y^2)"; break;
      case 7: f_print = f_7; f_name = "k = 7; f (x) = 1 / (25 * (x^2 + y^2) + 1)"; break;
    }
}

/// change current function for drawing
void Window::change_func ()
{
  if (!calculation_is_ready)
    {
      QMessageBox::information (0, "Message", "Calculation is running.");
      return;
    }
  // calculation_is_ready = 0;
  func_id = (func_id + 1) % 8;
  set_func ();
  pthread_args[0].func_id = func_id;
  pthread_args[0].set_func ();
  
  requst_to_calculate ();
}

void Window::change_nx_ny_up ()
{
  if (!calculation_is_ready)
    {
      QMessageBox::information (0, "Message", "Calculation is running.");
      return;
    }
  //printf("[k = 0]: change_nx_ny_up: %ld %ld\n", n, n_mas);
  calculation_is_ready = 0;
  /*if (nx >= n_max || ny >= n_max)
    return;
  
  if (2 * nx > n_max)
    nx = n_max;
  else
    nx *= 2;
  if (2 * ny > n_max)
    ny = n_max;
  else
    ny *= 2;*/
  /*if (mx * 2  >= m_max || my * 2 >= m_max)
    return;*/
  nx *= 2;
  ny *= 2;
  pthread_args[0].nx = nx;
  pthread_args[0].ny = ny;
  requst_to_calculate ();
}

void Window::change_nx_ny_down ()
{
  if (!calculation_is_ready)
    {
      QMessageBox::information (0, "Message", "Calculation is running.");
      return;
    }
  if (nx >= 2 || ny >= 2)
    {
      calculation_is_ready = 0;
      if (nx >= 2)
        nx /= 2;
      if (ny >= 2)
        ny /= 2;
      pthread_args[0].nx = nx;
      pthread_args[0].ny = ny;
      requst_to_calculate ();
    }
}

void Window::change_mx_my_up ()
{
  /*if (is_painting)
    {
      QMessageBox::warning (0, "Message", "Please, wait");
      return;
    }*/
  if (mx * 2  >= m_max || my * 2 >= m_max)
    return;
  mx *= 2;
  my *= 2;
  update ();
}

void Window::change_mx_my_down ()
{
  if (mx >= 2 && my >= 2)
    {
      mx /= 2;
      my /= 2;
      update ();
    }
}

void Window::set_view ()
{
  switch (view_id)
    {
      case 0:
        print_graph_func = 1;
        print_graph_aprox = 0;
        print_graph_diff = 0;
        break;
      case 1:
        print_graph_func = 0;
        print_graph_aprox = 1;
        print_graph_diff = 0;
        break;
      case 2:
        print_graph_func = 0;
        print_graph_aprox = 0;
        print_graph_diff = 1;
        break;
    }
}

void Window::change_view ()
{
  view_id = (view_id + 1) % 3;
  set_view ();
  update ();
}

void Window::change_zoom_up ()
{
  if (zoom < (int) (8 * sizeof(int)-1))
    {
      zoom += 1;
      update ();
    }
}
void Window::change_zoom_down ()
{
  if (zoom > 0)
    {
      zoom -= 1;
      update ();
    }
}

void Window::change_disturbance_up ()
{
  if (!calculation_is_ready)
    {
      QMessageBox::information (0, "Message", "Calculation is running.");
      return;
    }
  disturbance += 1;
  pthread_args[0].disturbance = disturbance;
  requst_to_calculate (); 
}
void Window::change_disturbance_down ()
{
  if (!calculation_is_ready)
    {
      QMessageBox::information (0, "Message", "Calculation is running.");
      return;
    }
  disturbance -= 1;
  pthread_args[0].disturbance = disturbance;
  requst_to_calculate ();
}

/*void Window::recalculate ()
{
  static pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
  static pthread_cond_t cond = PTHREAD_COND_INITIALIZER;
  pthread_mutex_lock (&mutex);
  need_to_calculate = 1;
  pthread_cond_broadcast(&cond);
  pthread_mutex_unlock(&mutex);
  
  pthread_mutex_lock (&mutex);
  while (!calculation_is_ready)
  {
      pthread_cond_wait(&cond,&mutex);
  }
  pthread_mutex_unlock (&mutex);
}*/
  

QPointF Window::l2g (double x_loc, double y_loc)
{
  double a1 = 0.5 * (a + b - (b - a) / (1ll << zoom));
  double b1 = 0.5 * (a + b + (b - a) / (1ll << zoom));
  double c1 = 0.5 * (c + d - (d - c) / (1ll << zoom));
  double d1 = 0.5 * (c + d + (d - c) / (1ll << zoom));
  double x_gl = (x_loc - a1) / (b1 - a1) * width ();
  double y_gl = (d1 - y_loc) / (d1 - c1) * height ();
  return QPointF (x_gl, y_gl);
}
QColor Window::l2g_colour (double z_loc, double z_min, double z_max)
{
  int z_gl = (z_loc - z_min) / (z_max - z_min) * 255;
  if (z_gl < 0)
    z_gl = 0;
  if (z_gl > 255)
    z_gl = 255;
  return QColor (100/*127 - z_gl/2*/,z_gl/*z_gl*/,255 - z_gl,255);
}

void Window::timerIsTimeout ()
{
  if (calculation_is_ready)
    {
      timer->stop();
      need_to_calculate = 0;
      x_mas = pthread_args[0].x;
      if (p >= 2)
        {
          pthread_args[0].it = pthread_args[1].it;
          pthread_args[0].t1 = pthread_args[1].t1;
          pthread_args[0].t2 = pthread_args[1].t2;
          pthread_args[0].r1 = pthread_args[1].r1;
          pthread_args[0].r2 = pthread_args[1].r2;
          pthread_args[0].r3 = pthread_args[1].r3;
          pthread_args[0].r4 = pthread_args[1].r4;
        }
      delete_memory_for_print ();
      set_memory_from_calc_to_print ();
      
      func_id_print = func_id;
      set_func_print ();
      nx_print = nx;
      ny_print = ny;
      disturbance_print = disturbance;
      
      update();
    }
}

/// render graph
void Window::paintEvent (QPaintEvent * /* event */)
{
  //is_painting = 1;
  if (x_mas == nullptr && !print_graph_func)
    return;
  
  QPainter painter (this);
  QPen pen_black (Qt::black, 0, Qt::SolidLine);
  painter.setPen (pen_black);
  
  double x1, x2, y1, y2, z;
  double max_z = 0, min_z = 0;
  double delta_x = (b - a) / (1ll << zoom) / mx;
  double delta_y = (d - c) / (1ll << zoom) / my;
  double delta_z;
  double a1 = 0.5 * (a + b - (b - a) / (1ll << zoom));
  double b1 = 0.5 * (a + b + (b - a) / (1ll << zoom));
  double c1 = 0.5 * (c + d - (d - c) / (1ll << zoom));
  double d1 = 0.5 * (c + d + (d - c) / (1ll << zoom));
  
  double max_abs_f = max_f_ab (a, b, c, d, func_id, f_print);
  int is_valid;
  
  char buf[100];
    
  // calculate min and max for current function
  if (print_graph_func)
    {
      max_z = min_z = 0;
      is_valid = 0;
      for (x1 = a1; x1 - b1 < 1.e-6; x1 += delta_x)
        {
          for (y1 = c1; y1 - d1 < 1.e-6; y1 += delta_y)
            {
              //z = f (x1, y1);
              z = f_dist (nx_print, ny_print, a, b, c, d, (x1-a)/(b-a)*nx_print, (y1-c)/(d-c)*ny_print, f_print, disturbance_print, max_abs_f);
              if (!is_valid || z < min_z)
                {
                  min_z = z;
                  is_valid = 1;
                  //printf("fu = %10.3e, min_z1 = %10.3e, x1 = %10.3e\n", y1, min_y, x1);
                }
              if (!is_valid || z > max_z)
                {
                  max_z = z;
                  is_valid = 1;
                  //printf("fu = %10.3e, max_z1 = %10.3e, x1 = %10.3e\n", y1, max_y1, x1);
                }
            }
        }
    }
  else
  if (print_graph_aprox)
    {
      max_z = min_z = 0;
      is_valid = 0;
      for (x1 = a1; x1 - b1 < 1.e-6; x1 += delta_x)
        {
          for (y1 = c1; y1 - d1 < 1.e-6; y1 += delta_y)
            {
              z = p_f (nx_print, ny_print, a, b, c, d, x_mas, x1, y1);
              if (!is_valid || z < min_z)
                {
                  min_z = z;
                  is_valid = 1;
                  //printf("a1 = %10.3e, min_y1 = %10.3e, x1 = %10.3e\n", 
                  //       y1, min_y, x1);
                }
              if (!is_valid || z > max_z)
                {
                  max_z = z;
                  is_valid = 1;
                }
            }
        }
    }
  else
  if (print_graph_diff)
    {
      max_z = min_z = 0;
      for (x1 = a1; x1 - b1 < 1.e-6; x1 += delta_x)
        {
          for (y1 = c1; y1 - d1 < 1.e-6; y1 += delta_y)
            {
              z = f_dist (nx_print, ny_print, a, b, c, d, (x1-a)/(b-a)*nx_print, (y1-c)/(d-c)*ny_print, f_print, disturbance_print, max_abs_f)
                 - p_f (nx_print, ny_print, a, b, c, d, x_mas, x1, y1);
              if (!is_valid || z < min_z)
                {
                  min_z = z;
                  //printf("a1 = %10.3e, min_y1 = %10.3e, x1 = %10.3e\n", 
                  //       y1, min_y, x1);
                }
              if (!is_valid || z > max_z)
                {
                  max_z = z;
                  is_valid = 1;
                }
            }
        }
    }
  
  if (fabs (max_z - min_z) < 1e-10)
    delta_z = 1e-12;
  else
    delta_z = 0.01 * (max_z - min_z);
    
  min_z -= delta_z;
  max_z += delta_z;
  
  
  QPointF points[4];

  painter.setPen(QPen(Qt::NoPen));
  if (print_graph_func)
    {
      x1 = a1;
      for (x2 = x1 + delta_x; x2 - b1 < 1.e-6; x2 += delta_x) 
        {
          y1 = c1;
          for (y2 = y1 + delta_y; y2 - d1 < 1.e-6; y2 += delta_y) 
            {
              z = f_print (x1 + 1./3*delta_x, y1 + 2./3*delta_y);
              painter.setBrush (L2GCOLOUR(z));
              points[0] = L2G(x1, y1);
              points[1] = L2G(x1, y2);
              points[2] = L2G(x2, y2);
              painter.drawPolygon(points, 3);
              
              z = f_print (x1 + 2./3*delta_x, y1 + 1./3*delta_y);
              
              painter.setBrush (L2GCOLOUR(z));
              points[0] = L2G(x1, y1);
              points[1] = L2G(x2, y1);
              points[2] = L2G(x2, y2);
              painter.drawPolygon(points, 3);

              y1 = y2;
            }
          x1 = x2;
        }
    }
  else
  if (print_graph_aprox)
    {
      x1 = a1;
      for (x2 = x1 + delta_x; x2 - b1 < 1.e-6; x2 += delta_x) 
        {
          y1 = c1;
          for (y2 = y1 + delta_y; y2 - d1 < 1.e-6; y2 += delta_y) 
            {
              z = p_f (nx_print,ny_print, a,b,c,d, x_mas, x1 + 1./3*delta_x, y1 + 2./3*delta_y);
              painter.setBrush (L2GCOLOUR(z));
              points[0] = L2G(x1, y1);
              points[1] = L2G(x1, y2);
              points[2] = L2G(x2, y2);
              painter.drawPolygon(points, 3);
              
              z = p_f (nx_print,ny_print, a,b,c,d, x_mas, x1 + 2./3*delta_x, y1 + 1./3*delta_y);
              painter.setBrush (L2GCOLOUR(z));
              points[0] = L2G(x1, y1);
              points[1] = L2G(x2, y1);
              points[2] = L2G(x2, y2);
              painter.drawPolygon(points, 3);

              y1 = y2;
            }
          x1 = x2;
        }
    }
  else
  if (print_graph_diff)
    {
      x1 = a1;
      for (x2 = x1 + delta_x; x2 - b1 < 1.e-6; x2 += delta_x) 
        {
          y1 = c1;
          for (y2 = y1 + delta_y; y2 - d1 < 1.e-6; y2 += delta_y) 
            {
              z = f_print (x1 + 1./3*delta_x, y1 + 2./3*delta_y) 
                - p_f (nx_print,ny_print, a,b,c,d, x_mas, x1 + 1./3*delta_x, y1 + 2./3*delta_y);
              painter.setBrush (L2GCOLOUR(z));
              points[0] = L2G(x1, y1);
              points[1] = L2G(x1, y2);
              points[2] = L2G(x2, y2);
              painter.drawPolygon(points, 3);
              //fabs
              z = f_print (x1 + 2./3*delta_x, y1 + 1./3*delta_y) 
                - p_f (nx_print,ny_print, a,b,c,d, x_mas, x1 + 2./3*delta_x, y1 + 1./3*delta_y);
              painter.setBrush (L2GCOLOUR(z));
              // local coords are converted to draw coords
              points[0] = L2G(x1, y1);
              points[1] = L2G(x2, y1);
              points[2] = L2G(x2, y2);
              painter.drawPolygon(points, 3);

              y1 = y2;
            }
          x1 = x2;
        }
    }
  
  printf("+------------------------------------------------+\n");
  //printf("min_y1 = %10.3e, max_y1 = %10.3e %d\n", min_y, max_y, width ());
  printf("Max_of_func = %10.3e\n", max_abs_f);
  printf("Max_fabs = %10.3e\n", (max_z > -min_z ? max_z : -min_z));
  
  painter.setPen (pen_black);
  painter.setBrush (QColor(255, 255, 255, 64));
  points[0] = QPointF (width() - 3, 3);
  points[1] = QPointF (width() - 3, 25);
  points[2] = QPointF (width() - 110, 25);
  points[3] = QPointF (width() - 110, 3);
  painter.drawPolygon(points, 4);
  
  //painter.drawText (width() - 120, 20, "Legend:");
  if (print_graph_func)
    painter.drawText (width() - 105, 20, "Function");
  else
  if (print_graph_aprox)
    painter.drawText (width() - 105, 20, "Aproximation");
  else
  if (print_graph_diff)
    painter.drawText (width() - 105, 20, "Difference");
  
  int len;
  if (func_id_print == 7)
    len = 310;
  else
  if (func_id_print == 4)
    len = 238;
  else
  if (func_id_print == 6)
    len = 232;
  else
    len = 212;
  painter.setBrush (QColor(255, 255, 255, 64));
  points[0] = QPointF (3, 3);
  points[1] = QPointF (3, 185);
  points[2] = QPointF (len, 185);
  points[3] = QPointF (len, 3);
  painter.drawPolygon(points, 4);
  
  // render function name
  painter.setPen ("black");
  painter.drawText (5, 20, f_name);
  
  sprintf(buf, "Max_of_func = %10.3e", max_abs_f);
  painter.drawText (5, 40, buf);
  sprintf(buf, "Max_fabs = %10.3e", (max_z > -min_z ? max_z : -min_z));
  painter.drawText (5, 60, buf);
  sprintf(buf, "Fmax = %10.3e", max_z);
  painter.drawText (5, 80, buf);
  sprintf(buf, "Fmin = %10.3e", min_z);
  painter.drawText (5, 100, buf);
  sprintf(buf, "Zoom: s = %d", zoom);
  painter.drawText (5, 120, buf);
  sprintf(buf, "nx/ny: %ld/%ld", nx_print, ny_print);
  painter.drawText (5, 140, buf);
  sprintf(buf, "mx/my: %zu/%zu", mx, my);
  painter.drawText (5, 160, buf);
  sprintf(buf, "disturbance: p = %d", disturbance_print);
  painter.drawText (5, 180, buf);
  
  // draw axis
  painter.setPen (pen_black);
  if (c1 <= 0 && 0 <= d1)
    painter.drawLine (L2G(a1, 0), L2G(b1, 0));
  if (a1 <= 0 && 0 <= b1)
    painter.drawLine (L2G(0, c1), L2G(0, d1));
  
  //sprintf(buf, "(%10.3e, %10.3e)", x2, y2);
  //painter.drawText (L2G(x2, y2), buf);
  //printf("pr_d2: x2 = %10.3e, y2 = %10.3e\n", x2, y2);
  
  
  /*if (print_graph_1 && n <= n_mas_1)
    {
      painter.setPen (pen_blue);
      for (size_t i = 0; i < n; i++) 
        {
          if (a1 <= x_t01[i] && x_t01[i] <= b1)
            y1 = Pf_01 (x_t01[i], a, b, n, x_t01, alpha);
          // local coords are converted to draw coords
          painter.drawEllipse (L2G(x_t01[i], y1), 3, 3);
        }
    }
  if (print_graph_2)
    {
      painter.setPen (pen_red);
      for (size_t i = 0; i < n; i++) 
        {
          if (a1 <= x_t02[i] && x_t02[i] <= b1)
            y1 = Pf_02 (x_t02[i], a, b, n, x_t02, mas_4n);
          // local coords are converted to draw coords
          painter.drawEllipse (L2G(x_t02[i], y1), 3, 3);
        }
    }*/

  it = pthread_args[0].it;
  t1 = pthread_args[0].t1;
  t2 = pthread_args[0].t2;
  r1 = pthread_args[0].r1;
  r2 = pthread_args[0].r2;
  r3 = pthread_args[0].r3;
  r4 = pthread_args[0].r4;
  printf (
"%s : Task = %d R1 = %e R2 = %e R3 = %e R4 = %e T1 = %.2f T2 = %.2f\
 It = %d E = %e K = %d Nx = %zu Ny = %zu P = %d\n",
prog_name, 7, r1, r2, r3, r4, t1, t2, it, eps, func_id_print, nx_print, ny_print, p);
  printf (
"%s : Task = %d R1 = %e R2 = %e R3 = %e R4 = %e T1 = %.2f T2 = %.2f\
 It = %d E = %e K = %d Nx = %zu Ny = %zu P = %d\n",
prog_name, 7, r1, r2, r3, r4, t1, t2, it, eps, func_id, nx, ny, p);
  printf("+------------------------------------------------+\n");
  /*printf (
"%s : Task = %d R1 = %e R2 = %e R3 = %e R4 = %e T1 = %.2f T2 = %.2f\
 It = %d E = %e K = %d Nx = %d Ny = %d P = %d\n",
prog_name, 7, r1, r2, r3, r4, t1, t2, it, eps, func_id, nx, ny, p);*/
  
  //is_painting = 0;
}
