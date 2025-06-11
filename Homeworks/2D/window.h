#ifndef WINDOW_H
#define WINDOW_H

#include <QtWidgets/QtWidgets>

#include "functions.h"

class qtMainWindow : public

class Window : public QWidget {
    Q_OBJECT

private:
	QWidget *widget;
    const char* f_name;
	char* file_name;
    int type_of_graph;
	
	double f_max;
	double f_min;
	double f_abs;
	
    double a;
    double b;
	double c;
    double d;
    int nx;
	int ny;
	int mx;
	int my;
	int func_id;
	double eps;
	int maxit;
	int p_thread;
	int s;
    int p;
	
	double* A;
	int* I;
	double* x;
	double* B;

	int nx_print = -1;
	int ny_print = -1;
	double a_print;
	double b_print;
	double c_print;
	double d_print;
	// double* A_print = nullptr;
	// int* I_print = nullptr;
	double* x_print = nullptr;
	// double* B_print = nullptr;
	// double* r_print = nullptr;
	// double* u_print = nullptr;
	// double* v_print = nullptr;
	
	double* r;
	double* u;
	double* v;
	
	double r1;
	double r2;
	double r3;
	double r4;
	
	double t1;
	double t2;
	
	int it;
	
    double (*f)(double, double);
	thread_data* thrd_data;
	pthread_t* tid;
	
	pthread_cond_t cond;
	pthread_mutex_t mutex;

public:
    Window(QWidget *parent);
	~Window();

    QSize minimumSizeHint() const;
    QSize sizeHint() const;

    int parse_command_line(int argc, char *argv[]);
    QPointF l2g(double x_loc, double y_loc, double y_min, double y_max);
    
public slots:
	void update_image_bounds();
	void select_func(int func_id);
	int read_args(char* argv[]);
	bool is_threads_ready();
	void waiting_threads();
	int memory_realloc();
	void update_thread_data();
	void save_prev_results();


    void change_f(); // 0
    void change_show_mode(); // 1
    void increase_visible_area(); // 2
    void decrease_visible_area(); // 3
    void increase_triangulation(); // 4
    void decrease_triangulation(); // 5
    void increase_protrusion(); // 6
    void decrease_protrusion(); // 7
	void increase_m(); // 8
	void decrease_m(); // 9
	void close();
    
protected:
	
	QPointF l2g(double x_loc, double y_loc);
	void draw_triangle(QPointF p_1, QPointF p_2, QPointF p_3, QPainter* painter, QColor color);
	QColor get_graph_color(double value, double max_value, double min_value);
	void draw_f(QPainter* painter);
	double get_aprx_value(double x, double y);
	void draw_aprx(QPainter* painter);
	void draw_res(QPainter* painter);
    void paintEvent(QPaintEvent *event);
	void closeEvent(QCloseEvent *event);
};

double f_0(double /* x */, double /* y */);
double f_1(double x, double /* y */);
double f_2(double /* x */, double y);
double f_3(double x, double y);
double f_4(double x, double y);
double f_5(double x, double y);
double f_6(double x, double y);
double f_7(double x, double y);

#endif
