#ifndef WINDOW_H
#define WINDOW_H

#include <QtWidgets/QtWidgets>

class Window : public QWidget
{
  Q_OBJECT

private:
  int func_id;
  int graph_type;
  const char *f_name;
  const char *f_mode;
  double a;
  double b;
  double res1;
  double res2;
  int s = 1;
  int p = 0;
  int n;
  double (*f) (double);
  double *lagrange_x;
  double *lagrange_c;
  double *cubic_x;
  double *cubic_c;
  double *cubic_d;
  double min_y = 0;
  double max_y = 0;
  double max_s = 0;
  double min_s = 0;

public:
  QPainter painter;

  Window (QWidget *parent);
  ~Window ();

  QSize minimumSizeHint () const;
  QSize sizeHint () const;

  int parse_command_line (int argc, char *argv[]);
  QPointF l2g (double x_loc, double y_loc, double y_min, double y_max);
public slots:

  void expand_n ();
  void compress_n ();

  void increase_p ();
  void decrease_p ();

  void expand_bounds ();
  void compress_bounds ();

  void change_func ();
  void change_graph_type ();

  void paintLagrange (QPainter &painter, QPen pen);
  void paintQubic (QPainter &painter, QPen pen);
  void paintLagrResidual (QPainter &painter, QPen pen, bool pr = false);
  void paintCubiResidual (QPainter &painter, QPen pen, bool pr = false);

protected:
  void paintEvent (QPaintEvent *event);
};

#endif
