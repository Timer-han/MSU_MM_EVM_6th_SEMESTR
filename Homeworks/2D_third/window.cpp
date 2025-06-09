#include <QPainter>
#include <stdio.h>

#include "window.h"
#include "reduce.h"

#define L2G(X, Y) (l2g((X), (Y), y_min, y_max))

double f_0(double /* x */, double /* y */)
{
	return 1;
}

double f_1(double x, double /* y */)
{
	return x;
}

double f_2(double /* x */, double y)
{
	return y;
}

double f_3(double x, double y)
{
	return x + y;
}

double f_4(double x, double y)
{
	return sqrt(x * x + y * y);
}

double f_5(double x, double y)
{
	return x * x + y * y;
}

double f_6(double x, double y)
{
	return exp(x * x - y * y);
}

double f_7(double x, double y)
{
	return 1. / (25 *(x * x + y * y) + 1);
}

void Window::update_f_min_max_abs() {	
	double hx = (b - a) / mx;
	double hy = (d - c) / my;
	
	f_max = f(a + hx * (1. / 3), c + hy * (2. / 3));
	f_min = f(a + hx * (1. / 3), c + hy * (2. / 3));
	
	for (int i = 0; i < mx; ++i) {
		for (int j = 0; j < my; ++j) {
			f_max = std::max(f_max, f(a + hx * (i + 1. / 3), c + hy * (j + 2. / 3)));
			f_min = std::min(f_min, f(a + hx * (i + 1. / 3), c + hy * (j + 2. / 3)));
			
			f_max = std::max(f_max, f(a + hx * (i + 2. / 3), c + hy * (j + 1. / 3)));
			f_min = std::min(f_min, f(a + hx * (i + 2. / 3), c + hy * (j + 1. / 3)));
		}
	}
	
	f_abs = std::max(std::abs(f_max), std::abs(f_min)); 
	
	if (f_max - f_min <= 0) {
		f_min -= 1.;
	}
	// double delta_f = 0.01 * (f_max - f_min);
	// f_max += delta_f;
	// f_min -= delta_f;
}

void Window::select_func(int func_id) {
	switch (func_id) {
        case 0:
            f_name = "f (x) = 1";
            f = f_0;
            break;
        case 1:
            f_name = "f (x) = x";
            f = f_1;
            break;
        case 2:
            f_name = "f (x) = y";
            f = f_2;
            break;
        case 3:
            f_name = "f (x) = x + y";
            f = f_3;
            break;
        case 4:
            f_name = "f (x) = sqrt(x * x + y * y)";
            f = f_4;
            break;
        case 5:
            f_name = "f (x) = x * x + y * y";
            f = f_5;
            break;
        case 6:
            f_name = "f (x) = exp(x * x - y * y)";
            f = f_6;
            break;
		case 7:
            f_name = "f (x) = 1. / (25 *(x * x + y * y) + 1)";
            f = f_7;
            break;
    }
}

Window::Window(QWidget *parent) : QWidget(parent) {
	widget = parent;
	
	s = 0;
	p = 0;
	func_id = 0;
	select_func(func_id);
	
	A = nullptr;
	I = nullptr;
	x = nullptr;
	B = nullptr;
	r = nullptr;
	u = nullptr;
	v = nullptr;
	
	cond = PTHREAD_COND_INITIALIZER;
	mutex = PTHREAD_MUTEX_INITIALIZER;
}

Window::~Window() {
	if (A) delete[] A;
	if (I) delete[] I;
	if (x) delete[] x;
	if (B) delete[] B;
	if (r) delete[] r;
	if (u) delete[] u;
	if (v) delete[] v;
	if (data_arr) delete[] data_arr;
	if (tid) delete[] tid;
}

void Window::Finish() {
	if (threads_are_ready()) {
		widget->close();
	} else {
        QMessageBox::warning(0, "Warning!", "The calculations are not finished yet!");
	}
}

QSize Window::minimumSizeHint() const {
	return QSize(100, 100);
}

QSize Window::sizeHint() const {
	return QSize(1000, 1000);
}

int Window::ReadData(char* argv[]) {
	file_name = argv[0];
	if (sscanf(argv[1], "%lf", &a) != 1 || sscanf(argv[2], "%lf", &b) != 1 || b - a < 1.e-6 
		|| sscanf(argv[3], "%lf", &c) != 1 || sscanf(argv[4], "%lf", &d) != 1 || d - c < 1.e-6 
		|| sscanf(argv[5], "%d", &nx) != 1 || sscanf(argv[6], "%d", &ny) != 1 || nx <= 0 || ny <= 0
		|| sscanf(argv[7], "%d", &mx) != 1 || sscanf(argv[8], "%d", &my) != 1 || mx <= 0 || my <= 0
		|| sscanf(argv[9], "%d", &func_id) != 1 || func_id < 0 || 7 < func_id
		|| sscanf(argv[10], "%lf", &eps) != 1 
		|| sscanf(argv[11], "%d", &maxit) != 1
		|| sscanf(argv[12], "%d", &p_thread) != 1) {
        return -1;
	}
	return 0;
}

bool Window::threads_are_ready() {
	for (int i = 0; i < p_thread; ++i) {
		if (!data_arr[i].ready) {
			return false;
		}
	}
	
	return true;
}

void Window::waiting_threads() {
    if (threads_are_ready()) {
        printf ("%s : Task = %d R1 = %e R2 = %e R3 = %e R4 = %e T1 = %.2f T2 = %.2f It = %d E = %e K = %d Nx = %d Ny = %d P = %d\n", 
				file_name, 1, data_arr->r1, data_arr->r2, data_arr->r3, data_arr->r4, data_arr->t1, data_arr->t2, data_arr->it, eps, func_id, nx, ny, p);
        update();
    } else {
        QTimer::singleShot(200, this, &Window::waiting_threads);
    }
}

int Window::update_memory() {
    if (A) delete[] A;
	if (I) delete[] I;
    if (x) delete[] x;
    if (B) delete[] B;
    if (r) delete[] r;
    if (u) delete[] u;
    if (v) delete[] v;

    if (allocate_msr_matrix(nx, ny, &A,  &I) != 0) {
		return -1;
	}
	
    int n = (nx + 1) * (ny + 1);
	
    fill_I(nx, ny, I);
	x = new double[n];
	B = new double[n];
	r = new double[n];
	u = new double[n];
	v = new double[n];
	
	for (int i = 0; i < n; ++i) {
		x[i] = 0;
		r[i] = 0;
		u[i] = 0;
		v[i] = 0;
	}
	
	return 0;
}

void Window::update_thread_data() {
	for (int i = 0; i < p_thread; ++i) {
        data_arr[i].a = a;
        data_arr[i].b = b;
        data_arr[i].c = c;
        data_arr[i].d = d;
		
		data_arr[i].nx = nx;
		data_arr[i].ny = ny;
		
		data_arr[i].f = f;
		data_arr[i].p_mid = p;
		data_arr[i].f_abs = f_abs;
		
		data_arr[i].A = A;
		data_arr[i].I = I;
		data_arr[i].x = x;
		data_arr[i].B = B;
		
		data_arr[i].r = r;
		data_arr[i].u = u;
		data_arr[i].v = v;
    }
}

int Window::parse_command_line(int argc, char* argv[]) {
    if (argc != 13) {
        return -1;
    }

    if (ReadData(argv)) {
        return -2;
	}
	
	InitReduceSum(p_thread);
	
	update_memory();
	
	select_func(func_id);
	type_of_graph = 0;
	update_f_min_max_abs();
	
	data_arr = new thread_data[p_thread];
    tid = new pthread_t[p_thread];
    for (int i = 0; i < p_thread; ++i) {
        data_arr[i].a = a;
        data_arr[i].b = b;
        data_arr[i].c = c;
        data_arr[i].d = d;
		
		data_arr[i].nx = nx;
		data_arr[i].ny = ny;
		
		data_arr[i].f = f;
		data_arr[i].p_mid = p;
		data_arr[i].f_abs = f_abs;
		
		data_arr[i].eps = eps;
		data_arr[i].maxit = maxit;
		data_arr[i].p = p_thread;
		data_arr[i].k = i;
		
		data_arr[i].A = A;
		data_arr[i].I = I;
		data_arr[i].x = x;
		data_arr[i].B = B;
		
		data_arr[i].r = r;
		data_arr[i].u = u;
		data_arr[i].v = v;
		
        data_arr[i].mutex = &mutex;
        data_arr[i].cond = &cond;
		
        data_arr[i].ready = false;
        data_arr[i].finish = false;
    }
    
	for (int i = 0; i < p_thread; ++i) {
		pthread_create(&tid[i], 0, thread_func, data_arr + i);
	}
	
	QTimer::singleShot(200, this, &Window::waiting_threads);
	
	return 0;
}

void Window::ChangeFunc() {
	if (threads_are_ready()) {
		func_id = (func_id + 1) % 8;
		select_func(func_id);
		update_memory();
		update_f_min_max_abs();
		update_thread_data();

        for (int i = 0; i < p_thread; ++i) {
            data_arr[i].ready = false;
        }

		pthread_cond_broadcast(&cond);

        QTimer::singleShot(200, this, &Window::waiting_threads);

    } else {
        QMessageBox::warning(0, "Warning!", "The calculations are not finished yet!");
    }
}

void Window::ChangeTypeOfGraph() {
    type_of_graph  = (type_of_graph + 1) % 3;
    update();
}

void Window::ExtendArea() {
	if (threads_are_ready()) {
		double mid = (a + b) / 2;
		double len = (b - a);
		a = mid - len;
		b = mid + len;
		
		mid = (c + d) / 2;
		len = (d - c);
		c = mid - len;
		d = mid + len;
		
		++s;
		
		update_memory();
		update_f_min_max_abs();
		update_thread_data();
		
        for (int i = 0; i < p_thread; ++i) {
            data_arr[i].ready = false;
        }

		pthread_cond_broadcast(&cond);

        QTimer::singleShot(200, this, &Window::waiting_threads);

    } else {
        QMessageBox::warning(0, "Warning!", "The calculations are not finished yet!");
    }
}

void Window::CompressArea() {
	if (threads_are_ready()) {
		double mid = (a + b) / 2;
		double len = (b - a);
		a = mid - len / 4;
		b = mid + len / 4;
		
		mid = (c + d) / 2;
		len = (d - c);
		c = mid - len / 4;
		d = mid + len / 4;
		
		--s;
		
		update_memory();
		update_f_min_max_abs();
		update_thread_data();
		
        for (int i = 0; i < p_thread; ++i) {
            data_arr[i].ready = false;
        }

		pthread_cond_broadcast(&cond);

        QTimer::singleShot(200, this, &Window::waiting_threads);

    } else {
        QMessageBox::warning(0, "Warning!", "The calculations are not finished yet!");
    }
}

void Window::IncreaseN() {
	if (threads_are_ready()) {
		nx = 2 * nx;
		ny = 2 * ny;
		
		update_memory();
		update_f_min_max_abs();
		update_thread_data();
		
        for (int i = 0; i < p_thread; ++i) {
            data_arr[i].ready = false;
        }

		pthread_cond_broadcast(&cond);

        QTimer::singleShot(200, this, &Window::waiting_threads);
		
    } else {
        QMessageBox::warning(0, "Warning!", "The calculations are not finished yet!");
    } 
}

void Window::DecreaseN() {
	if (threads_are_ready()) {
		nx = nx / 2;
		ny = ny / 2;
		
		nx = std::max(nx, 1);
		ny = std::max(ny, 1);
		
		update_memory();
		update_f_min_max_abs();
		update_thread_data();
		
		for (int i = 0; i < p_thread; ++i) {
			data_arr[i].nx = nx;
			data_arr[i].ny = ny;
			data_arr[i].ready = false;
		}

		pthread_cond_broadcast(&cond);

        QTimer::singleShot(200, this, &Window::waiting_threads);

    } else {
        QMessageBox::warning(0, "Warning!", "The calculations are not finished yet!");
    }
}

void Window::IncreaseFuncMiddle() {
    if (threads_are_ready()) {
		++p;
		
		update_memory();
		update_f_min_max_abs();
		update_thread_data();
		
		for (int i = 0; i < p_thread; ++i) {
			data_arr[i].ready = false;
		}

		pthread_cond_broadcast(&cond);

        QTimer::singleShot(200, this, &Window::waiting_threads);

    } else {
        QMessageBox::warning(0, "Warning!", "The calculations are not finished yet!");
    }
}

void Window::DecreaseFuncMiddle() {
	if (threads_are_ready()) {
		--p;
		
		update_memory();
		update_f_min_max_abs();
		update_thread_data();
		
		for (int i = 0; i < p_thread; ++i) {
			data_arr[i].ready = false;
		}

		pthread_cond_broadcast(&cond);

        QTimer::singleShot(200, this, &Window::waiting_threads);

    } else {
        QMessageBox::warning(0, "Warning!", "The calculations are not finished yet!");
    }
}

void Window::IncreaseM() {
	mx = 2 * mx;
	my = 2 * my;
	
	mx = std::min(1000, mx);
	my = std::min(1000, my);
	
	update_f_min_max_abs();
	
	update();
}

void Window::DecreaseM() {
	mx = mx / 2;
	my = my / 2;
	
	mx = std::max(1, mx);
	my = std::max(1, my);
	
	update_f_min_max_abs();
	
	update();
}

QPointF Window::l2g(double x_loc, double y_loc) {
	double x_gl = (x_loc - a) / (b - a) * width();
	double y_gl = (d - y_loc) / (d - c) * height();
	return QPointF(x_gl, y_gl);
}

void Window::DrawTriangle(QPointF p_1, QPointF p_2, QPointF p_3, QPainter* painter, QColor color) {
	QPainterPath path;
	path.moveTo(p_1);
	path.lineTo(p_2);
	path.lineTo(p_3);
	path.closeSubpath();
	
	painter->fillPath(path, QBrush{color});
}

QColor Window::GetColor(double value, double max_value, double min_value) {
	if (value < min_value) {
		value = min_value;
	}
	if (value > max_value) {
		value = max_value;
	}
	
	double coeff;
	if (max_value <= min_value) {
		coeff = 1;
	} else {
		coeff = (value - min_value) / (max_value - min_value);
	}
	
	if (type_of_graph == 0) {
		return QColor(0, 0, 255 * coeff);
	}
	if (type_of_graph == 1) {
		return QColor(0, 255 * coeff, 0);
	}
	return QColor(255 * coeff, 0, 0);
}

void Window::DrawFunc(QPainter* painter) {
	double hx = (b - a) / mx;
    double hy = (d - c) / my;
    
	QColor color;
	
    for (int i = 0; i < mx; ++i) {
        for (int j = 0; j < my; ++j) {
			color = GetColor(f(a + hx * (i + 1. / 3), c + hy * (j + 2. / 3)), f_max, f_min);
			DrawTriangle(l2g(a + hx * i, c + hy * j), l2g(a + hx * i, c + hy * (j + 1)), l2g(a + hx * (i + 1), c + hy * (j + 1)), painter, color);
			
			color = GetColor(f(a + hx * (i + 2. / 3), c + hy * (j + 1. / 3)), f_max, f_min);
			DrawTriangle(l2g(a + hx * i, c + hy * j), l2g(a + hx * (i + 1), c + hy * j), l2g(a + hx * (i + 1), c + hy * (j + 1)), painter, color);
        }
    }
    
    int line = 0;
	char max_mod_buf[20];
	double max_mod = std::max(std::abs(f_max), std::abs(f_min));
	
    painter->setPen("white");
    painter->drawText(10, line += 20, (std::string("k = ") + std::to_string(func_id)).c_str());
    painter->drawText(10, line += 20, f_name);
	
	sprintf(max_mod_buf, "%10.3e", max_mod);
	painter->drawText(10, line += 20, (std::string("|F|_max = ") + std::string(max_mod_buf)).c_str());
	printf("|F|_max =%10.3e\n", max_mod);
	
    painter->drawText(10, line += 20, (std::string("nx = ") + std::to_string(nx)).c_str());
	painter->drawText(10, line += 20, (std::string("ny = ") + std::to_string(ny)).c_str());
    painter->drawText(10, line += 20, (std::string("mx = ") + std::to_string(mx)).c_str());
	painter->drawText(10, line += 20, (std::string("my = ") + std::to_string(my)).c_str());
    painter->drawText(10, line += 20, (std::string("s = ") + std::to_string(s)).c_str());
    painter->drawText(10, line += 20, (std::string("p = ") + std::to_string(p)).c_str());
}

double Window::ApproximationPointValue(double x_in, double y_in) {
	int i = (int)((x_in - a) / (b - a) * nx);
	int j = (int)((y_in - c) / (d - c) * ny);
	
	i = std::min(std::max(i, 0), nx - 1);
    j = std::min(std::max(j, 0), ny - 1);
	
	int l;
    ij2l(nx, ny, i, j, l);
	
	if (x_in - a - ((double)i) / nx * (b - a) <= 0 && y_in - c - ((double)j) / ny * (d - c) <= 0) {
		return x[l];
	}
	
	if (x_in - a - ((double)i) / nx * (b - a) <= 0) {
		double point_0_0 = x[l];
		double point_0_1 = x[l + nx + 1];
		return point_0_0 + (point_0_1 - point_0_0) * ((y_in - c - ((double)j) / ny * (d - c)) / ((d - c) / ny));
	}
	
	if (y_in - c - ((double)j) / ny * (d - c) <= 0) {
		double point_0_0 = x[l];
		double point_1_0 = x[l + 1];
		return point_0_0 + (point_1_0 - point_0_0) * ((x_in - a - ((double)i) / nx * (b - a)) / ((b - a) / nx));
	}
	
	double coeff = (y_in - c - ((double)j) / ny * (d - c)) / (x_in - a - ((double)i) / nx * (b - a));
	
	if (coeff > ((d - c) / ny) / ((b - a) / nx)) {
		double point_0_0 = x[l];
		double point_0_1 = x[l + nx + 1];
		double point_1_1 = x[l + 1 + nx + 1];
		double point_1_0 = point_0_0 + point_1_1 - point_0_1;
		
		double temp_x = point_0_0 + (point_1_0 - point_0_0) * ((x_in - a - ((double)i) / nx * (b - a)) / ((b - a) / nx));
		double temp_y = point_0_0 + (point_0_1 - point_0_0) * ((y_in - c - ((double)j) / ny * (d - c)) / ((d - c) / ny));
		
		return temp_x + temp_y - point_0_0;
	}
	
	double point_0_0 = x[l];
	double point_1_0 = x[l + 1];
	double point_1_1 = x[l + 1 + nx + 1];
	double point_0_1 = point_0_0 + point_1_1 - point_1_0;
	
	double temp_x = point_0_0 + (point_1_0 - point_0_0) * ((x_in - a - ((double)i) / nx * (b - a)) / ((b - a) / nx));
	double temp_y = point_0_0 + (point_0_1 - point_0_0) * ((y_in - c - ((double)j) / ny * (d - c)) / ((d - c) / ny));
	
	return temp_x + temp_y - point_0_0;
}

void Window::DrawApproximation(QPainter* painter) {
	double hx = (b - a) / mx;
    double hy = (d - c) / my;
	
	double app_max = ApproximationPointValue(a + hx * (1. / 3), c + hy * (2. / 3));
	double app_min = ApproximationPointValue(a + hx * (1. / 3), c + hy * (2. / 3));
	
	QColor color;
	double temp;
	
	for (int i = 0; i < mx; ++i) {
		for (int j = 0; j < my; ++j) {
			temp = ApproximationPointValue(a + hx * (i + 1. / 3), c + hy * (j + 2. / 3));
			app_max = std::max(app_max, temp);
			app_min = std::min(app_min, temp);
			
			temp = ApproximationPointValue(a + hx * (i + 2. / 3), c + hy * (j + 1. / 3));
			app_max = std::max(app_max, temp);
			app_min = std::min(app_min, temp);
		}
	}
	
	for (int i = 0; i < mx; ++i) {
        for (int j = 0; j < my; ++j) {
			color = GetColor(ApproximationPointValue(a + hx * (i + 1. / 3), c + hy * (j + 2. / 3)), app_max, app_min);
			DrawTriangle(l2g(a + hx * i, c + hy * j), l2g(a + hx * i, c + hy * (j + 1)), l2g(a + hx * (i + 1), c + hy * (j + 1)), painter, color);
			
			color = GetColor(ApproximationPointValue(a + hx * (i + 2. / 3), c + hy * (j + 1. / 3)), app_max, app_min);
			DrawTriangle(l2g(a + hx * i, c + hy * j), l2g(a + hx * (i + 1), c + hy * j), l2g(a + hx * (i + 1), c + hy * (j + 1)), painter, color);
        }
    }
    
	int line = 0;
	char max_mod_buf[20];
	double max_mod = std::max(std::abs(app_max), std::abs(app_min));
	
    painter->setPen("white");
    painter->drawText(10, line += 20, (std::string("k = ") + std::to_string(func_id)).c_str());
    painter->drawText(10, line += 20, f_name);
	
	sprintf(max_mod_buf, "%10.3e", max_mod);
	painter->drawText(10, line += 20, (std::string("|Pf|_max = ") + std::string(max_mod_buf)).c_str());
	printf("|Pf|_max =%10.3e\n", max_mod);
	
    painter->drawText(10, line += 20, (std::string("nx = ") + std::to_string(nx)).c_str());
	painter->drawText(10, line += 20, (std::string("ny = ") + std::to_string(ny)).c_str());
    painter->drawText(10, line += 20, (std::string("mx = ") + std::to_string(mx)).c_str());
	painter->drawText(10, line += 20, (std::string("my = ") + std::to_string(my)).c_str());
    painter->drawText(10, line += 20, (std::string("s = ") + std::to_string(s)).c_str());
    painter->drawText(10, line += 20, (std::string("p = ") + std::to_string(p)).c_str());
}

void Window::DrawResidual(QPainter* painter) {
	double hx = (b - a) / mx;
    double hy = (d - c) / my;
	
	double res_max = std::abs(f(a + hx * (1. / 3), c + hy * (2. / 3)) - ApproximationPointValue(a + hx * (1. / 3), c + hy * (2. / 3)));
	double res_min = std::abs(f(a + hx * (1. / 3), c + hy * (2. / 3)) - ApproximationPointValue(a + hx * (1. / 3), c + hy * (2. / 3)));
	
	QColor color;
	double temp;
	
	for (int i = 0; i < mx; ++i) {
		for (int j = 0; j < my; ++j) {
			temp = std::abs(f(a + hx * (i + 1. / 3), c + hy * (j + 2. / 3)) - ApproximationPointValue(a + hx * (i + 1. / 3), c + hy * (j + 2. / 3)));
			res_max = std::max(res_max, temp);
			res_min = std::min(res_min, temp);
			
			temp = std::abs(f(a + hx * (i + 2. / 3), c + hy * (j + 1. / 3)) - ApproximationPointValue(a + hx * (i + 2. / 3), c + hy * (j + 1. / 3)));
			res_max = std::max(res_max, temp);
			res_min = std::min(res_min, temp);
		}
	}
	
	for (int i = 0; i < mx; ++i) {
        for (int j = 0; j < my; ++j) {
			color = GetColor(std::abs(f(a + hx * (i + 1. / 3), c + hy * (j + 2. / 3)) - ApproximationPointValue(a + hx * (i + 1. / 3), c + hy * (j + 2. / 3))), res_max, res_min);
			DrawTriangle(l2g(a + hx * i, c + hy * j), l2g(a + hx * i, c + hy * (j + 1)), l2g(a + hx * (i + 1), c + hy * (j + 1)), painter, color);
			
			color = GetColor(std::abs(f(a + hx * (i + 2. / 3), c + hy * (j + 1. / 3)) - ApproximationPointValue(a + hx * (i + 2. / 3), c + hy * (j + 1. / 3))), res_max, res_min);
			DrawTriangle(l2g(a + hx * i, c + hy * j), l2g(a + hx * (i + 1), c + hy * j), l2g(a + hx * (i + 1), c + hy * (j + 1)), painter, color);
        }
    }
    
	int line = 0;
	char max_mod_buf[20];
	double max_mod = std::max(std::abs(res_max), std::abs(res_min));
	
    painter->setPen("white");
    painter->drawText(10, line += 20, (std::string("k = ") + std::to_string(func_id)).c_str());
    painter->drawText(10, line += 20, f_name);
	
	sprintf(max_mod_buf, "%10.3e", max_mod);
	painter->drawText(10, line += 20, (std::string("|F - Pf|_max = ") + std::string(max_mod_buf)).c_str());
	printf("|F - Pf|_max =%10.3e\n", max_mod);
	
    painter->drawText(10, line += 20, (std::string("nx = ") + std::to_string(nx)).c_str());
	painter->drawText(10, line += 20, (std::string("ny = ") + std::to_string(ny)).c_str());
    painter->drawText(10, line += 20, (std::string("mx = ") + std::to_string(mx)).c_str());
	painter->drawText(10, line += 20, (std::string("my = ") + std::to_string(my)).c_str());
    painter->drawText(10, line += 20, (std::string("s = ") + std::to_string(s)).c_str());
    painter->drawText(10, line += 20, (std::string("p = ") + std::to_string(p)).c_str());
}

void Window::paintEvent(QPaintEvent * /* event */) {
    QPainter painter(this);
	
	if (type_of_graph == 0) {
		DrawFunc(&painter);
	}
	if (type_of_graph == 1) {
		DrawApproximation(&painter);
	}
	if (type_of_graph == 2) {
		DrawResidual(&painter);
	}
}
