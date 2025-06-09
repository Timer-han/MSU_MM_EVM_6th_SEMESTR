#ifndef WINDOW_H
#define WINDOW_H

#include <QtWidgets/QtWidgets>

#include "functions.h"

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
	thread_data* data_arr;
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
	void update_f_min_max_abs();
	void select_func(int func_id);
	int ReadData(char* argv[]);
	bool threads_are_ready();
	void waiting_threads();
	int update_memory();
	void update_thread_data();
	
    void ChangeFunc(); // 0
    void ChangeTypeOfGraph(); // 1
    void ExtendArea(); // 2
    void CompressArea(); // 3
    void IncreaseN(); // 4
    void DecreaseN(); // 5
    void IncreaseFuncMiddle(); // 6
    void DecreaseFuncMiddle(); // 7
	void IncreaseM(); // 8
	void DecreaseM(); // 9
	void Finish();
    
protected:
	
	QPointF l2g(double x_loc, double y_loc);
	void DrawTriangle(QPointF p_1, QPointF p_2, QPointF p_3, QPainter* painter, QColor color);
	QColor GetColor(double value, double max_value, double min_value);
	void DrawFunc(QPainter* painter);
	double ApproximationPointValue(double x, double y);
	void DrawApproximation(QPainter* painter);
	void DrawResidual(QPainter* painter);
    void paintEvent(QPaintEvent *event);
};

#endif
