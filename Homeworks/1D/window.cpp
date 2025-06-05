#include <QPainter>
#include <stdio.h>
#include <cmath>

#include "window.h"

#define EPS 1.e-6
#define DEFAULT_A -10
#define DEFAULT_B 10
#define DEFAULT_N 10
#define L2G(X,Y) (l2g ((X), (Y), min_y, max_y))


template <typename T>
int sgn (T val)
{
  return (T(0) < val) - (val < T(0));
}

static
double f_0 (double /* x */)
{
  return 1;
}

static
double f_1 (double x)
{
  return x;
}

static
double f_2 (double x)
{
  return x * x;
}

static
double f_3 (double x)
{
  return x * x * x;
}

static
double f_4 (double x)
{
  return x * x * x * x;
}

static
double f_5 (double x)
{
  return std::exp (x);
}

static
double f_6 (double x)
{
  return 1 / (25 * x * x + 1);
}

static
double get_f (
  double (*f) (double),
  double x,
  double a,
  double b,
  int n,
  double max_y,
  double min_y,
  int p)
{
  double df = std::max(fabs (max_y), fabs (min_y)) * 0.1 * p;
  if (fabs (x - a - (b - a) / n * (n / 2)) > EPS)
    {
      df = 0.;
    }
  return f (x) + df;
}

static
double derivative (
  double (*f) (double),
  double x,
  double dx)
{
  if (f == f_1)
    return 1;
  if (f == f_2)
    return 2 * x;
  if (f == f_3)
    return 3 * x * x;
  if (f == f_4)
    return 4 * x * x * x;
  if (f == f_5)
    return std::exp (x);
  if (f == f_6)
    return -50 * x / (25 * x * x + 1) / (25 * x * x + 1);
  return (f (x + dx) - f (x - dx)) / (2 * dx);
}

static
void fill_newton_coef_calc (
  double (*f) (double),
  double a,
  double b,
  int n,
  double *x,
  double *result,
  double max_y,
  double min_y,
  int p)
{
  if (n > 50) return;
  double dx = (b - a) / n;
  for (int i = 0; i < 2 * (n + 1); i++)
    {
      x[i] = a + (i / 2) * dx;
      result[i] = get_f (f, x[i], a, b, n, max_y, min_y, p);
      // printf("f(%6.1e)=%6.1e\n", x[i], result[i]);
    }
  // printf("\n");

  for (int i = 2 * n; i > 0; i--)
    {
      if (i % 2 == 0)
        result[i] = (result[i] - result[i - 1]) / (x[i] - x[i - 1]);
      else
        result[i] = derivative (f, x[i], dx * 1e-5);
      // printf("c[%d]=%6.1e\n", i, result[i]);
    }
  // printf("\n");

  for (int i = 1; i < 2 * n; i++)
    {
      for (int j = 2 * n; j >= i; j--)
        {
          result[j + 1] = (result[j + 1] - result[j]) / (x[j + 1] - x[j - i]);
          // printf("c[%d]=%8.3e | ", j + 1, result[j + 1]);
        }
      // printf("\n");
    }
}



static
double get_lagrange_value (
  double point,
  int n,
  double *x,
  double *result)
{
  double value = result[2 * n - 1];
  for (int i = 2 * n; i >= 0; i--)
    {
      value = result[i] + (point - x[i]) * value;
    }
  return value;
}

static
void fill_qubic_coeff_arr(
  double (*f) (double),
  int n,
  double a,
  double b,
  double *x,
  double *c,
  double *d,
  double max_y,
  double min_y,
  int p)
{
	int i;
	double d1, d2, dx = (b - a) / (n - 1);
	double f1, f2, f_1, f_2;

  for (i = 0; i < n + 3; i++)
    {
      x[i] = a + (i) * dx;
    }

	for (i = 1; i < n - 1; i++)
    {
      f_1 = get_f (f, x[i],     a, b, n, max_y, min_y, p);
      f_2 = get_f (f, x[i + 1], a, b, n, max_y, min_y, p);
      d1 = -(f_1 - f_2) / dx;
      f_2 = get_f (f, x[i - 1], a, b, n, max_y, min_y, p);
      d2 = -(f_2 - f_1) / dx;
      if (d1 * d2 > 0) d[i] = sgn (d2) * std::min(fabs(d1), fabs(d2));
      else d[i] = 0.0;
    }

	f1 = f(x[0] - dx);
	f2 = f(x[n - 1] + dx);

	d1 = (f(x[0]) - f1) / dx;
	d2 = (f(x[1]) - f(x[0])) / dx;
	if (d1 * d2 > 0) d[0] = sgn (d2) * std::min(fabs(d1), fabs(d2));
	else d[0] = 0.0;

	d1 = (f(x[n - 1]) - f(x[n - 2])) / dx;
	d2 = (f2 - f(x[n - 1])) / dx;
	if (d1 * d2 > 0) d[n - 1] = sgn (d2) * std::min(fabs(d1), fabs(d2));
	else d[n - 1] = 0.0;

	for (i = 0; i < n - 1; i++)
    {
      c[4 * i] = get_f (f, x[i], a, b, n, max_y, min_y, p);
      c[4 * i + 1] = d[i];
      f_1 = get_f (f, x[i + 1], a, b, n, max_y, min_y, p);
      f_2 = get_f (f, x[i], a, b, n, max_y, min_y, p);
      d1 = (f_1 - f_2) / dx;
      c[4 * i + 2] = (3 * d1 - 2 * d[i] - d[i + 1]) / dx;
      c[4 * i + 3] = (d[i] + d[i + 1] - 2 * d1) / dx / dx;
      
      // c[4 * i] = f(x[i]);
      // c[4 * i + 1] = d[i];
      // d1 = (f(x[i + 1]) - f(x[i])) / dx;
      // c[4 * i + 2] = (3 * d1 - 2 * d[i] - d[i + 1]) / dx;
      // c[4 * i + 3] = (d[i] + d[i + 1] - 2 * d1) / dx / dx;
    }
}

static
double get_appr_cubic_value(double point, int n, double *c, double *x)
{
	int i = 0, dn = n;

  
  if (point < x[0] + EPS) dn = i + 2;
  if (point > x[n - 1] - EPS) i = dn - 2;
  while (dn - i > 4)
  {
    if (fabs (point - x[(dn + i) / 2]) < EPS) break;
    else if (point < x[(dn + i) / 2]) dn = (dn + i) / 2;
    else i = (dn + i) / 2;
  }
  
  for (; i < dn; i++)
    if (point <= x[i + 1] && point > x[i])
      break;

	return c[4 * i]
       + c[4 * i + 1] * (point - x[i])
       + c[4 * i + 2] * (point - x[i]) * (point - x[i])
       + c[4 * i + 3] * (point - x[i]) * (point - x[i]) * (point - x[i]);
}

Window::Window (QWidget *parent)
  : QWidget (parent)
{
  a = DEFAULT_A;
  b = DEFAULT_B;
  n = DEFAULT_N;

  func_id = 0;
  graph_type = 0;
  p = 0;

}

Window::~Window ()
{
  if (lagrange_x) delete[] lagrange_x;
  if (lagrange_c) delete[] lagrange_c;
  if (cubic_x) delete[] cubic_x;
  if (cubic_c) delete[] cubic_c;
  if (cubic_d) delete[] cubic_d;
}

QSize Window::minimumSizeHint () const
{
  return QSize (100, 100);
}

QSize Window::sizeHint () const
{
  return QSize (1000, 1000);
}

int Window::parse_command_line (int argc, char *argv[])
{
  if (argc == 1)
    return 0;

  if (argc == 2)
    return -1;

  if (   sscanf (argv[1], "%lf", &a) != 1
      || sscanf (argv[2], "%lf", &b) != 1
      || b - a < 1.e-6
      || (argc > 3 && sscanf (argv[3], "%d", &n) != 1)
      || n <= 0
      || (argc > 4 && sscanf (argv[4], "%d", &func_id) != 1)
      || func_id < 0)
    return -2;

  if (n > width ())
    n = width ();
  if (n < 4)
    n = 5;
  
  n--;

  lagrange_x = new double[2 * (n + 3)];
  lagrange_c = new double[2 * (n + 3)];
  cubic_x = new double[2 * (n + 3)];
  cubic_d = new double[2 * (n + 3)];
  cubic_c = new double[4 * (n + 3)];

  memset(lagrange_x, 0, sizeof(double) * 2 * (n + 3));
  memset(lagrange_c, 0, sizeof(double) * 2 * (n + 3));
  memset(cubic_x, 0, sizeof(double) * 2 * (n + 3));
  memset(cubic_d, 0, sizeof(double) * 2 * (n + 3));
  memset(cubic_c, 0, sizeof(double) * 4 * (n + 3));

  graph_type--;
  func_id--;
  change_func ();
  change_graph_type ();

  return 0;
}

void Window::increase_p ()
{
  p++;

  graph_type--;
  func_id--;
  change_func ();
  change_graph_type ();
}

void Window::decrease_p ()
{
  p--;

  graph_type--;
  func_id--;
  change_func ();
  change_graph_type ();
}

void Window::expand_bounds ()
{
  a *= 2;
  b *= 2;
  s++;

  graph_type--;
  func_id--;
  change_func ();
  change_graph_type ();
}

void Window::compress_bounds ()
{
  a /= 2;
  b /= 2;
  s--;

  graph_type--;
  func_id--;
  change_func ();
  change_graph_type ();
}

void Window::expand_n ()
{
  n = n * 2 + 1;

  if (lagrange_x) delete[] lagrange_x;
  if (lagrange_c) delete[] lagrange_c;
  if (cubic_x) delete[] cubic_x;
  if (cubic_c) delete[] cubic_c;
  if (cubic_d) delete[] cubic_d;

  lagrange_x = new double[2 * (n + 3)];
  lagrange_c = new double[2 * (n + 3)];
  cubic_x = new double[2 * (n + 3)];
  cubic_d = new double[2 * (n + 3)];
  cubic_c = new double[4 * (n + 3)];

  memset(lagrange_x, 0, sizeof(double) * 2 * (n + 3));
  memset(lagrange_c, 0, sizeof(double) * 2 * (n + 3));
  memset(cubic_x, 0, sizeof(double) * 2 * (n + 3));
  memset(cubic_d, 0, sizeof(double) * 2 * (n + 3));
  memset(cubic_c, 0, sizeof(double) * 4 * (n + 3));

  graph_type--;
  func_id--;
  change_func ();
  change_graph_type ();
}

void Window::compress_n ()
{
  if (n <= 4)
    return;
  n = (n + 1) / 2 - 1;
  
  if (lagrange_x) delete[] lagrange_x;
  if (lagrange_c) delete[] lagrange_c;
  if (cubic_x) delete[] cubic_x;
  if (cubic_c) delete[] cubic_c;
  if (cubic_d) delete[] cubic_d;
  
  lagrange_x = new double[2 * (n + 3)];
  lagrange_c = new double[2 * (n + 3)];
  cubic_x = new double[2 * (n + 3)];
  cubic_d = new double[2 * (n + 3)];
  cubic_c = new double[4 * (n + 3)];

  memset(lagrange_x, 0, sizeof(double) * 2 * (n + 3));
  memset(lagrange_c, 0, sizeof(double) * 2 * (n + 3));
  memset(cubic_x, 0, sizeof(double) * 2 * (n + 3));
  memset(cubic_d, 0, sizeof(double) * 2 * (n + 3));
  memset(cubic_c, 0, sizeof(double) * 4 * (n + 3));

  graph_type--;
  func_id--;

  change_func ();
  change_graph_type ();
}

void Window::change_graph_type ()
{
  graph_type = (graph_type + 1) % 4;

  switch (graph_type)
    {
      case 0:
      f_mode = "Lagrange";
      break;
      case 1:
      f_mode = "Cubic";
      break;
      case 2:
      f_mode = "Lagrange & Cubic";
        break;
      case 3:
        f_mode = "Residual";
        break;
    }
  
  update ();
}

/// change current function for drawing
void Window::change_func ()
{
  func_id = (func_id + 1) % 7;

  switch (func_id)
    {
      case 0:
        f_name = "f (x) = 1";
        f = f_0;
        break;
      case 1:
        f_name = "f (x) = x";
        f = f_1;
        break;
      case 2:
        f_name = "f (x) = x ^ 2";
        f = f_2;
        break;
      case 3:
        f_name = "f (x) = x ^ 3";
        f = f_3;
        break;
      case 4:
        f_name = "f (x) = x ^ 4";
        f = f_4;
        break;
      case 5:
        f_name = "f (x) = exp (x)";
        f = f_5;
        break;
      case 6:
        f_name = "f (x) = 1 / (25 * x ^ 2 + 1)";
        f = f_6;
        break;
    }

  fill_newton_coef_calc (f, a, b, n, lagrange_x, lagrange_c, max_y, min_y, p);
  fill_qubic_coeff_arr(f, n + 1, a, b, cubic_x, cubic_c, cubic_d, max_y, min_y, p);
  
  update ();
}

QPointF Window::l2g (double x_loc, double y_loc, double y_min, double y_max)
{
  double x_gl = (x_loc - a) / (b - a) * width ();
  double y_gl = (y_max - y_loc) / (y_max - y_min) * height ();
  return QPointF (x_gl, y_gl);
}


// draw lagrange approximated line for graph
void Window::paintLagrange (
  QPainter &painter, QPen pen)
{
  double x1, x2, y1, y2;

  if (n > 50) return;

  QPen old_pen = painter.pen ();
  painter.setPen (pen);

  
  x1 = a;
  y1 = get_lagrange_value (x1, n, lagrange_x, lagrange_c);
  for (x2 = x1 + (b - a) / width (); x2 - b < 1.e-6; x2 += (b - a) / width ()) 
    {
      y2 = get_lagrange_value (x2, n, lagrange_x, lagrange_c);

      painter.drawLine (L2G(x1, y1), L2G(x2, y2));

      x1 = x2, y1 = y2;
    }
  x2 = b;
  y2 = get_lagrange_value (x2, n, lagrange_x, lagrange_c);
  painter.drawLine (L2G(x1, y1), L2G(x2, y2));

  painter.setPen (old_pen);
}

void Window::paintQubic (
  QPainter &painter, QPen pen)
{
  double x1, x2, y1, y2;

  painter.setPen (pen);
  
  x1 = a;
  y1 = f(a);
  for (x2 = x1 + (b - a) / width (); x2 - b < 1.e-6; x2 += (b - a) / width ()) 
    {
      y2 = get_appr_cubic_value (x2, n + 1, cubic_c, cubic_x);
      // printf("coord: %lf %lf\n", x2, y2);
      // local coords are converted to draw coords
      painter.drawLine (L2G(x1, y1), L2G(x2, y2));

      x1 = x2, y1 = y2;
    }
}

void Window::paintLagrResidual (
  QPainter &painter, QPen pen, bool pr)
{
  if (n > 50) return;
  double x1, x2, y1, y2;
  double bd = (b - a) * 20 / (width ());
  auto old_pen = painter.pen ();

  res1 = 0;
  painter.setPen (pen);
  
  x1 = a;
  y1 = 0.;
  for (x2 = x1 + (b - a) / width (); x2 - b < 1.e-6; x2 += (b - a) / width ()) 
    {
      y2 = fabs(
        f(x2) - get_lagrange_value (x2, n, lagrange_x, lagrange_c)
      );
      if (x2 > a + bd && x2 < b - bd) res1 = std::max(res1, y2);
      if (pr) painter.drawLine (L2G(x1, y1), L2G(x2, y2));

      x1 = x2, y1 = y2;
    }

  painter.setPen (old_pen);
}

void Window::paintCubiResidual (
  QPainter &painter, QPen pen, bool pr)
{
  double x1, x2, y1, y2;
  double bd = (b - a) * 20 / (width ());
  auto old_pen = painter.pen ();
  
  res2 = 0;
  painter.setPen (pen);

  
  x1 = a;
  y1 = 0;
  for (x2 = x1 + (b - a) / width (); x2 - b < 1.e-6; x2 += (b - a) / width ()) 
    {
      y2 = fabs(
        f(x2) - get_appr_cubic_value (x2, n + 1, cubic_c, cubic_x)
      );
      if (x2 > a + bd && x2 < b - bd) res2 = std::max(res2, y2);
      if (pr) painter.drawLine (L2G(x1, y1), L2G(x2, y2));

      x1 = x2, y1 = y2;
    }
  
  painter.setPen (old_pen);
}

/// render graph
void Window::paintEvent (QPaintEvent * /* event */)
{  
  QPainter painter (this);
  double x1, x2, y1, y2;
  double delta_y;
  QPen pen_black(Qt::black, 1, Qt::SolidLine); 
  QPen pen_red(Qt::red, 1, Qt::SolidLine); 
  QPen pen_green(QColor(0, 150, 0), 1, Qt::SolidLine);
  QPen pen_blue(QColor(0, 0, 200), 1, Qt::SolidLine);

  painter.setPen (pen_black);

  // calculate min and max for current function
  max_y = min_y = 0;
  for (x1 = a; x1 - b < 1.e-6; x1 += (b - a) / width ())
    {
      y1 = f (x1);
      if (y1 < min_y)
        min_y = y1;
      if (y1 > max_y)
        max_y = y1;
    }

  delta_y = 0.01 * (max_y - min_y);
  min_y -= delta_y;
  max_y += delta_y;

  // draw source line for graph
  if (graph_type != 3)
    {
    x1 = a;
    y1 = f (x1);
    for (x2 = x1 + (b - a) / width (); x2 - b < 1.e-6; x2 += (b - a) / width ()) 
      {
        y2 = f (x2);
        // local coords are converted to draw coords
        painter.drawLine (L2G(x1, y1), L2G(x2, y2));

        x1 = x2, y1 = y2;
      }
    x2 = b;
    y2 = f (x2);
    painter.drawLine (L2G(x1, y1), L2G(x2, y2));
    }
  
  bool lprint = false;
  bool cprint = false;
  switch (graph_type)
    {
      case 0:
        paintLagrange (painter, pen_green);
        break;
      case 1:
        paintQubic (painter, pen_blue);
        break;
      case 2:
        paintLagrange (painter, pen_green);
        paintQubic (painter, pen_blue);
        break;
      case 3:
        lprint = true;
        cprint = true;
        break;
    }
  paintLagrResidual (painter, pen_green, lprint);
  paintCubiResidual (painter, pen_blue, cprint);
  
  

  // draw axis
  painter.setPen (pen_red);
  painter.drawLine (L2G(a, 0), L2G(b, 0));
  painter.drawLine (L2G(0, min_y), L2G(0, max_y));

  // render function name
  painter.setPen ("blue");
  painter.drawText (3, 20, f_name);
  painter.drawText (3, 40, f_mode);
  painter.drawText (3, 60, QString("n = %1").arg(n + 1));
  painter.drawText (3, 80, QString("a = %1, b = %2").arg(a).arg(b));
  painter.drawText (3, 100, QString("p = %1").arg(p));
  painter.setPen (pen_green);
  painter.drawText (3, 120, QString("res1 = %1").arg(res1));
  painter.setPen (pen_blue);
  painter.drawText (3, 140, QString("res1 = %1").arg(res2));

}
