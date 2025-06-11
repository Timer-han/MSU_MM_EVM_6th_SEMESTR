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

void Window::update_image_bounds() {	
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

	// if (I_print) delete[] I_print;
	// if (A_print) delete[] A_print;
	if (x_print) delete[] x_print;
	// if (B_print) delete[] B_print;
	// if (r_print) delete[] r_print;
	// if (u_print) delete[] u_print;
	// if (v_print) delete[] v_print;

	if (thrd_data) delete[] thrd_data;
	if (tid) delete[] tid;
}

void Window::close() {
	if (is_threads_ready()) {
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

int Window::read_args(char* argv[]) {
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
	a_print = a;
	b_print = b;
	c_print = c;
	d_print = d;
	nx_print = nx;
	
	return 0;
}

bool Window::is_threads_ready() {
	for (int i = 0; i < p_thread; ++i) {
		if (!thrd_data[i].ready) {
			return false;
		}
	}
	
	return true;
}

void Window::waiting_threads() {
    if (is_threads_ready()) {
        printf ("%s : Task = %d R1 = %e R2 = %e R3 = %e R4 = %e T1 = %.2f T2 = %.2f It = %d E = %e K = %d Nx = %d Ny = %d P = %d\n", 
				file_name, 1, thrd_data->r1, thrd_data->r2, thrd_data->r3, thrd_data->r4, thrd_data->t1, thrd_data->t2, thrd_data->it, eps, func_id, nx, ny, p);
        
		save_prev_results();
		update();
    } else {
        QTimer::singleShot(200, this, &Window::waiting_threads);
    }
}

int Window::memory_realloc() {
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

	// if (!I_print) {
	// 	int size = n + get_len_msr_all_diag(nx, ny) + 1;
	// 	I_print = new int[size];
	// 	for (int i = 0; i < size; ++i) I_print[i] = I[i];
	// }
	// if (!A_print) {
	// 	int size = n + get_len_msr_all_diag(nx, ny) + 1;
	// 	A_print = new double[size];
	// 	for (int i = 0; i < size; ++i) A_print[i] = A[i];
	// }
	if (!x_print) {
		x_print = new double[n];
		for (int i = 0; i < n; ++i) x_print[i] = x[i];
	}
	// if (!B_print) {
	// 	B_print = new double[n];
	// 	for (int i = 0; i < n; ++i) B_print[i] = B[i];
	// }
	// if (!r_print) {
	// 	r_print = new double[n];
	// 	for (int i = 0; i < n; ++i) r_print[i] = r[i];
	// }
	// if (!u_print) {
	// 	u_print = new double[n];
	// 	for (int i = 0; i < n; ++i) u_print[i] = u[i];
	// }
	// if (!v_print) {
	// 	v_print = new double[n];
	// 	for (int i = 0; i < n; ++i) v_print[i] = v[i];
	// }

	if (nx_print == -1) {
		nx_print = nx;
		ny_print = ny;
		a_print = a;
		b_print = b;
		c_print = c;
		d_print = d;
	}
	// else {
	// 	if (nx != nx_print || ny != ny_print || a != a_print || b != b_print || c != c_print || d != d_print) {
	// 		printf("Error: the size of the grid has changed!\n");
	// 		return -2;
	// 	}
	// }

	// update_thread_data();
	
	// for (int i = 0; i < p_thread; ++i) {
	// 	thrd_data[i].ready = true;
	// }
	
	return 0;
}

void Window::save_prev_results() {
	// printf("I'm here! %s:%d\n", __FILE__, __LINE__);
	// if (I_print) delete[] I_print;
	// if (A_print) delete[] A_print;
	if (x_print) delete[] x_print;
	// if (B_print) delete[] B_print;
	// if (r_print) delete[] r_print;
	// if (u_print) delete[] u_print;
	// if (v_print) delete[] v_print;

	// printf("I'm here! %s:%d\n", __FILE__, __LINE__);

	// I_print = I;
	// A_print = A;
	// x_print = x;
	// B_print = B;
	// r_print = r;
	// u_print = u;
	// v_print = v;

	int n = (nx + 1) * (ny + 1);
	// int size = n + get_len_msr_all_diag(nx, ny) + 1;

	// printf("I'm here! %s:%d\n", __FILE__, __LINE__);
	// I_print = new int[size];
	// for (int i = 0; i < size; ++i) I_print[i] = I[i];
	// printf("I'm here! %s:%d\n", __FILE__, __LINE__);

	// A_print = new double[size];
	// for (int i = 0; i < size; ++i) A_print[i] = A[i];

	x_print = new double[n];
	for (int i = 0; i < n; ++i) x_print[i] = x[i];

	// B_print = new double[n];
	// for (int i = 0; i < n; ++i) B_print[i] = B[i];

	// r_print = new double[n];
	// for (int i = 0; i < n; ++i) r_print[i] = r[i];

	// u_print = new double[n];
	// for (int i = 0; i < n; ++i) u_print[i] = u[i];

	// v_print = new double[n];
	// for (int i = 0; i < n; ++i) v_print[i] = v[i];

	// printf("I'm here! %s:%d\n", __FILE__, __LINE__);
	// I = nullptr;
	// A = nullptr;
	// x = nullptr;
	// B = nullptr;
	// r = nullptr;
	// u = nullptr;
	// v = nullptr;

	nx_print = nx;
	ny_print = ny;
	a_print = a;
	b_print = b;
	c_print = c;
	d_print = d;

}

void Window::update_thread_data() {
	for (int i = 0; i < p_thread; ++i) {
        thrd_data[i].a = a;
        thrd_data[i].b = b;
        thrd_data[i].c = c;
        thrd_data[i].d = d;
		
		thrd_data[i].nx = nx;
		thrd_data[i].ny = ny;
		
		thrd_data[i].f = f;
		thrd_data[i].p_mid = p;
		thrd_data[i].f_abs = f_abs;
		
		thrd_data[i].A = A;
		thrd_data[i].I = I;
		thrd_data[i].x = x;
		thrd_data[i].B = B;
		
		thrd_data[i].r = r;
		thrd_data[i].u = u;
		thrd_data[i].v = v;
    }
}

int Window::parse_command_line(int argc, char* argv[]) {
    if (argc != 13) {
        return -1;
    }

    if (read_args(argv)) {
        return -2;
	}
	
	init_reduce_sum(p_thread);
	
	memory_realloc();
	
	select_func(func_id);
	type_of_graph = 0;
	update_image_bounds();
	
	thrd_data = new thread_data[p_thread];
    tid = new pthread_t[p_thread];
    for (int i = 0; i < p_thread; ++i) {
        thrd_data[i].a = a;
        thrd_data[i].b = b;
        thrd_data[i].c = c;
        thrd_data[i].d = d;
		
		thrd_data[i].nx = nx;
		thrd_data[i].ny = ny;
		
		thrd_data[i].f = f;
		thrd_data[i].p_mid = p;
		thrd_data[i].f_abs = f_abs;
		
		thrd_data[i].eps = eps;
		thrd_data[i].max_it = maxit;
		thrd_data[i].p = p_thread;
		thrd_data[i].pi = i;
		
		thrd_data[i].A = A;
		thrd_data[i].I = I;
		thrd_data[i].x = x;
		thrd_data[i].B = B;
		
		thrd_data[i].r = r;
		thrd_data[i].u = u;
		thrd_data[i].v = v;
		
        thrd_data[i].mutex = &mutex;
        thrd_data[i].cond = &cond;
		
        thrd_data[i].ready = false;
        thrd_data[i].finish = false;
    }
    
	for (int i = 0; i < p_thread; ++i) {
		pthread_create(&tid[i], 0, thread_func, thrd_data + i);
	}
	
	QTimer::singleShot(200, this, &Window::waiting_threads);
	
	return 0;
}

void Window::change_f() {
	if (is_threads_ready()) {
        for (int i = 0; i < p_thread; ++i) {
            thrd_data[i].ready = false;
        }
		save_prev_results();
		func_id = (func_id + 1) % 8;
		select_func(func_id);
		memory_realloc();
		update_image_bounds();
		update_thread_data();

		for (int i = 0; i < p_thread; ++i) {
            thrd_data[i].ready = false;
        }

		pthread_cond_broadcast(&cond);

        QTimer::singleShot(200, this, &Window::waiting_threads);

    } else {
        QMessageBox::warning(0, "Warning!", "The calculations are not finished yet!");
    }
}

void Window::change_show_mode() {
    type_of_graph  = (type_of_graph + 1) % 3;
    update();
}

void Window::increase_visible_area() {
	if (is_threads_ready()) {
		for (int i = 0; i < p_thread; ++i) {
            thrd_data[i].ready = false;
        }
		save_prev_results();
		double mid = (a + b) / 2;
		double len = (b - a);
		a = mid - len;
		b = mid + len;
		
		mid = (c + d) / 2;
		len = (d - c);
		c = mid - len;
		d = mid + len;
		
		++s;
		
		// update();

		memory_realloc();
		update_image_bounds();
		update_thread_data();
		
        for (int i = 0; i < p_thread; ++i) {
            thrd_data[i].ready = false;
        }

		pthread_cond_broadcast(&cond);

        QTimer::singleShot(200, this, &Window::waiting_threads);

    } else {
        QMessageBox::warning(0, "Warning!", "The calculations are not finished yet!");
    }
}

void Window::decrease_visible_area() {
	if (is_threads_ready()) {
		for (int i = 0; i < p_thread; ++i) {
            thrd_data[i].ready = false;
        }
		save_prev_results();
		double mid = (a + b) / 2;
		double len = (b - a);
		a = mid - len / 4;
		b = mid + len / 4;
		
		mid = (c + d) / 2;
		len = (d - c);
		c = mid - len / 4;
		d = mid + len / 4;
		
		--s;
		
		// update();

		memory_realloc();
		update_image_bounds();
		update_thread_data();
		
        for (int i = 0; i < p_thread; ++i) {
            thrd_data[i].ready = false;
        }

		pthread_cond_broadcast(&cond);

        QTimer::singleShot(200, this, &Window::waiting_threads);

    } else {
        QMessageBox::warning(0, "Warning!", "The calculations are not finished yet!");
    }
}

void Window::increase_triangulation() {
	if (is_threads_ready()) {
		for (int i = 0; i < p_thread; ++i) {
            thrd_data[i].ready = false;
        }
		save_prev_results();
		nx = 2 * nx;
		ny = 2 * ny;
		
		// update();

		memory_realloc();
		update_image_bounds();
		update_thread_data();
		
        for (int i = 0; i < p_thread; ++i) {
            thrd_data[i].ready = false;
        }

		pthread_cond_broadcast(&cond);

        QTimer::singleShot(200, this, &Window::waiting_threads);
		
    } else {
        QMessageBox::warning(0, "Warning!", "The calculations in process!");
    } 
}

void Window::decrease_triangulation() {
	if (is_threads_ready()) {
		for (int i = 0; i < p_thread; ++i) {
            thrd_data[i].ready = false;
        }
		save_prev_results();
		nx = nx / 2;
		ny = ny / 2;
		
		nx = std::max(nx, 1);
		ny = std::max(ny, 1);
		
		// update();

		memory_realloc();
		update_image_bounds();
		update_thread_data();
		
		for (int i = 0; i < p_thread; ++i) {
			thrd_data[i].nx = nx;
			thrd_data[i].ny = ny;
			thrd_data[i].ready = false;
		}

		pthread_cond_broadcast(&cond);

        QTimer::singleShot(200, this, &Window::waiting_threads);

    } else {
        QMessageBox::warning(0, "Warning!", "The calculations are not finished yet!");
    }
}

void Window::increase_protrusion() {
    if (is_threads_ready()) {
		for (int i = 0; i < p_thread; ++i) {
            thrd_data[i].ready = false;
        }
		save_prev_results();
		++p;
		
		// update();
		
		memory_realloc();
		update_image_bounds();
		update_thread_data();
		
		for (int i = 0; i < p_thread; ++i) {
			thrd_data[i].ready = false;
		}

		pthread_cond_broadcast(&cond);

        QTimer::singleShot(200, this, &Window::waiting_threads);

    } else {
        QMessageBox::warning(0, "Warning!", "The calculations are not finished yet!");
    }
}

void Window::decrease_protrusion() {
	if (is_threads_ready()) {
		for (int i = 0; i < p_thread; ++i) {
            thrd_data[i].ready = false;
        }
		save_prev_results();
		--p;
		
		// update();

		memory_realloc();
		update_image_bounds();
		update_thread_data();
		
		for (int i = 0; i < p_thread; ++i) {
			thrd_data[i].ready = false;
		}

		pthread_cond_broadcast(&cond);

        QTimer::singleShot(200, this, &Window::waiting_threads);

    } else {
        QMessageBox::warning(0, "Warning!", "The calculations are not finished yet!");
    }
}

void Window::increase_m() {
	mx = 2 * mx;
	my = 2 * my;
	
	mx = std::min(1000, mx);
	my = std::min(1000, my);
	
	update_image_bounds();
	
	update();
}

void Window::decrease_m() {
	mx = mx / 2;
	my = my / 2;
	
	mx = std::max(1, mx);
	my = std::max(1, my);
	
	update_image_bounds();
	
	update();
}

QPointF Window::l2g(double x_loc, double y_loc) {
	double x_gl = (x_loc - a) / (b - a) * width();
	double y_gl = (d - y_loc) / (d - c) * height();
	return QPointF(x_gl, y_gl);
}

void Window::draw_triangle(QPointF p_1, QPointF p_2, QPointF p_3, QPainter* painter, QColor color) {
	QPointF points[3] = {p_1, p_2, p_3};
	painter->setPen(Qt::NoPen);
	painter->setBrush(QBrush(color));
	painter->drawPolygon(points, 3);
	// QPainterPath path;
	// path.moveTo(p_1);
	// path.lineTo(p_2);
	// path.lineTo(p_3);
	// path.closeSubpath();
	
	// painter->fillPath(path, QBrush{color});

}

QColor Window::get_graph_color(double value, double max_value, double min_value) {
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
		return QColor(0, 0, 200 * coeff);
	}
	if (type_of_graph == 1) {
		return QColor(0, 150 * coeff, 0);
	}
	return QColor(200 * coeff, 0, 0);
}

void Window::draw_f(QPainter* painter) {
	double hx = (b - a) / mx;
    double hy = (d - c) / my;
    
	QColor color;
	
    for (int i = 0; i < mx; ++i) {
        for (int j = 0; j < my; ++j) {
			color = get_graph_color(f(a + hx * (i + 1. / 3), c + hy * (j + 2. / 3)), f_max, f_min);
			draw_triangle(l2g(a + hx * i, c + hy * j), l2g(a + hx * i, c + hy * (j + 1)), l2g(a + hx * (i + 1), c + hy * (j + 1)), painter, color);
			
			color = get_graph_color(f(a + hx * (i + 2. / 3), c + hy * (j + 1. / 3)), f_max, f_min);
			draw_triangle(l2g(a + hx * i, c + hy * j), l2g(a + hx * (i + 1), c + hy * j), l2g(a + hx * (i + 1), c + hy * (j + 1)), painter, color);
        }
    }
    
    int line = 0;
	char max_mod_buf[20];
	double max_mod = std::max(std::abs(f_max), std::abs(f_min));
	
	QPointF points[4] = {
		QPointF(4, 4),
		QPointF(4, 220),
		QPointF(240, 220),
		QPointF(240, 4)
	};
	painter->setPen(QColor(255, 255, 255, 60));
	painter->setBrush(QColor(0, 0, 0, 100));
	painter->drawPolygon(points, 4);

    painter->setPen("white");
    painter->drawText(15, line += 20, (std::string("func_id = ") + std::to_string(func_id)).c_str());
    painter->drawText(15, line += 20, f_name);
	
	sprintf(max_mod_buf, "%10.3e", max_mod);
	printf("[+] ||F||_max =%10.3e\n", max_mod);
	painter->drawText(15, line += 20, (std::string("||F||_max = ") + std::string(max_mod_buf)).c_str());
	
    painter->drawText(15, line += 20, (std::string("nx = ") + std::to_string(nx)).c_str());
	painter->drawText(15, line += 20, (std::string("ny = ") + std::to_string(ny)).c_str());
    painter->drawText(15, line += 20, (std::string("mx = ") + std::to_string(mx)).c_str());
	painter->drawText(15, line += 20, (std::string("my = ") + std::to_string(my)).c_str());
    painter->drawText(15, line += 20, (std::string("s = ")  + std::to_string(s)).c_str());
    painter->drawText(15, line += 20, (std::string("p = ")  + std::to_string(p)).c_str());
	painter->drawText(15, line += 20, (std::string("Function")).c_str());
}

double Window::get_aprx_value(double px, double py) {
	int i = (int)((px - a_print) / (b_print - a_print) * nx_print);
	int j = (int)((py - c_print) / (d_print - c_print) * ny_print);
	
	i = std::min(std::max(i, 0), nx_print - 1);
    j = std::min(std::max(j, 0), ny_print - 1);
	
	int l;
    ij2l(nx_print, ny_print, i, j, l);
	
	if (px - a_print - ((double)i) / nx_print * (b_print - a_print) <= 0 && py - c_print - ((double)j) / ny_print * (d_print - c_print) <= 0) {
		return x_print[l];
	}
	
	if (px - a_print - ((double)i) / nx_print * (b_print - a_print) <= 0) {
		double point_0_0 = x_print[l];
		double point_0_1 = x_print[l + nx_print + 1];
		return point_0_0 + (point_0_1 - point_0_0) * ((py - c_print - ((double)j) / ny_print * (d_print - c_print)) / ((d_print - c_print) / ny_print));
	}
	
	if (py - c_print - ((double)j) / ny_print * (d_print - c_print) <= 0) {
		double point_0_0 = x_print[l];
		double point_1_0 = x_print[l + 1];
		return point_0_0 + (point_1_0 - point_0_0) * ((px - a_print - ((double)i) / nx_print * (b_print - a_print)) / ((b_print - a_print) / nx_print));
	}
	
	double coeff = (py - c_print - ((double)j) / ny_print * (d_print - c_print)) / (px - a_print - ((double)i) / nx_print * (b_print - a_print));
	
	if (coeff > ((d_print - c_print) / ny_print) / ((b_print - a_print) / nx_print)) {
		double point_0_0 = x_print[l];
		double point_0_1 = x_print[l + nx_print + 1];
		double point_1_1 = x_print[l + 1 + nx_print + 1];
		double point_1_0 = point_0_0 + point_1_1 - point_0_1;
		
		double buf_1 = point_0_0 + (point_1_0 - point_0_0) * ((px - a_print - ((double)i) / nx_print * (b_print - a_print)) / ((b_print - a_print) / nx_print));
		double buf_2 = point_0_0 + (point_0_1 - point_0_0) * ((py - c_print - ((double)j) / ny_print * (d_print - c_print)) / ((d_print - c_print) / ny_print));
		
		return buf_1 + buf_2 - point_0_0;
	}
	
	double point_0_0 = x_print[l];
	double point_1_0 = x_print[l + 1];
	double point_1_1 = x_print[l + 1 + nx_print + 1];
	double point_0_1 = point_0_0 + point_1_1 - point_1_0;
	
	double buf_1 = point_0_0 + (point_1_0 - point_0_0) * ((px - a_print - ((double)i) / nx_print * (b_print - a_print)) / ((b_print - a_print) / nx_print));
	double buf_2 = point_0_0 + (point_0_1 - point_0_0) * ((py - c_print - ((double)j) / ny_print * (d_print - c_print)) / ((d_print - c_print) / ny_print));
	
	return buf_1 + buf_2 - point_0_0;
}

void Window::draw_aprx(QPainter* painter) {
	double hx = (b - a) / mx;
    double hy = (d - c) / my;
	
	double app_max = get_aprx_value(a + hx * (1. / 3), c + hy * (2. / 3));
	double app_min = get_aprx_value(a + hx * (1. / 3), c + hy * (2. / 3));
	
	QColor color;
	double temp;
	
	for (int i = 0; i < mx; ++i) {
		for (int j = 0; j < my; ++j) {
			temp = get_aprx_value(a + hx * (i + 1. / 3), c + hy * (j + 2. / 3));
			app_max = std::max(app_max, temp);
			app_min = std::min(app_min, temp);
			
			temp = get_aprx_value(a + hx * (i + 2. / 3), c + hy * (j + 1. / 3));
			app_max = std::max(app_max, temp);
			app_min = std::min(app_min, temp);
		}
	}
	
	for (int i = 0; i < mx; ++i) {
        for (int j = 0; j < my; ++j) {
			color = get_graph_color(get_aprx_value(a + hx * (i + 1. / 3), c + hy * (j + 2. / 3)), app_max, app_min);
			draw_triangle(l2g(a + hx * i, c + hy * j), l2g(a + hx * i, c + hy * (j + 1)), l2g(a + hx * (i + 1), c + hy * (j + 1)), painter, color);
			
			color = get_graph_color(get_aprx_value(a + hx * (i + 2. / 3), c + hy * (j + 1. / 3)), app_max, app_min);
			draw_triangle(l2g(a + hx * i, c + hy * j), l2g(a + hx * (i + 1), c + hy * j), l2g(a + hx * (i + 1), c + hy * (j + 1)), painter, color);
        }
    }
    
	int line = 0;
	char max_mod_buf[20];
	double max_mod = std::max(std::abs(app_max), std::abs(app_min));

	QPointF points[4] = {
		QPointF(4, 4),
		QPointF(4, 220),
		QPointF(240, 220),
		QPointF(240, 4)
	};
	painter->setPen(QColor(255, 255, 255, 60));
	painter->setBrush(QColor(0, 0, 0, 100));
	painter->drawPolygon(points, 4);

    painter->setPen("white");
    painter->drawText(15, line += 20, (std::string("func_id = ") + std::to_string(func_id)).c_str());
    painter->drawText(15, line += 20, f_name);
	
	sprintf(max_mod_buf, "%10.3e", max_mod);
	printf("[+] ||P_f||_max =%10.3e\n", max_mod);
	painter->drawText(15, line += 20, (std::string("||P_f||_max = ") + std::string(max_mod_buf)).c_str());
	
    painter->drawText(15, line += 20, (std::string("nx = ") + std::to_string(nx)).c_str());
	painter->drawText(15, line += 20, (std::string("ny = ") + std::to_string(ny)).c_str());
    painter->drawText(15, line += 20, (std::string("mx = ") + std::to_string(mx)).c_str());
	painter->drawText(15, line += 20, (std::string("my = ") + std::to_string(my)).c_str());
    painter->drawText(15, line += 20, (std::string("s = ")  + std::to_string(s)).c_str());
    painter->drawText(15, line += 20, (std::string("p = ")  + std::to_string(p)).c_str());
	painter->drawText(15, line += 20, (std::string("Approximation")).c_str());
}

void Window::draw_res(QPainter* painter) {
	double hx = (b - a) / mx;
    double hy = (d - c) / my;
	
	double res_max = std::abs(f(a + hx * (1. / 3), c + hy * (2. / 3)) - get_aprx_value(a + hx * (1. / 3), c + hy * (2. / 3)));
	double res_min = std::abs(f(a + hx * (1. / 3), c + hy * (2. / 3)) - get_aprx_value(a + hx * (1. / 3), c + hy * (2. / 3)));
	
	QColor color;
	double temp;
	
	for (int i = 0; i < mx; ++i) {
		for (int j = 0; j < my; ++j) {
			temp = std::abs(f(a + hx * (i + 1. / 3), c + hy * (j + 2. / 3)) - get_aprx_value(a + hx * (i + 1. / 3), c + hy * (j + 2. / 3)));
			res_max = std::max(res_max, temp);
			res_min = std::min(res_min, temp);
			
			temp = std::abs(f(a + hx * (i + 2. / 3), c + hy * (j + 1. / 3)) - get_aprx_value(a + hx * (i + 2. / 3), c + hy * (j + 1. / 3)));
			res_max = std::max(res_max, temp);
			res_min = std::min(res_min, temp);
		}
	}
	
	for (int i = 0; i < mx; ++i) {
        for (int j = 0; j < my; ++j) {
			color = get_graph_color(std::abs(f(a + hx * (i + 1. / 3), c + hy * (j + 2. / 3)) - get_aprx_value(a + hx * (i + 1. / 3), c + hy * (j + 2. / 3))), res_max, res_min);
			draw_triangle(l2g(a + hx * i, c + hy * j), l2g(a + hx * i, c + hy * (j + 1)), l2g(a + hx * (i + 1), c + hy * (j + 1)), painter, color);
			
			color = get_graph_color(std::abs(f(a + hx * (i + 2. / 3), c + hy * (j + 1. / 3)) - get_aprx_value(a + hx * (i + 2. / 3), c + hy * (j + 1. / 3))), res_max, res_min);
			draw_triangle(l2g(a + hx * i, c + hy * j), l2g(a + hx * (i + 1), c + hy * j), l2g(a + hx * (i + 1), c + hy * (j + 1)), painter, color);
        }
    }
    
	int line = 0;
	char max_mod_buf[20];
	double max_mod = std::max(std::abs(res_max), std::abs(res_min));

	QPointF points[4] = {
		QPointF(4, 4),
		QPointF(4, 220),
		QPointF(240, 220),
		QPointF(240, 4)
	};
	painter->setPen(QColor(255, 255, 255, 60));
	painter->setBrush(QColor(0, 0, 0, 150));
	painter->drawPolygon(points, 4);
	
    painter->setPen("white");
    painter->drawText(15, line += 20, (std::string("func_id = ") + std::to_string(func_id)).c_str());
    painter->drawText(15, line += 20, f_name);
	
	sprintf(max_mod_buf, "%10.3e", max_mod);
	printf("[+] ||F - P_f||_max =%10.3e\n", max_mod);
	painter->drawText(15, line += 20, (std::string("||F - P_f||_max = ") + std::string(max_mod_buf)).c_str());
	
    painter->drawText(15, line += 20, (std::string("nx = ") + std::to_string(nx)).c_str());
	painter->drawText(15, line += 20, (std::string("ny = ") + std::to_string(ny)).c_str());
    painter->drawText(15, line += 20, (std::string("mx = ") + std::to_string(mx)).c_str());
	painter->drawText(15, line += 20, (std::string("my = ") + std::to_string(my)).c_str());
    painter->drawText(15, line += 20, (std::string("s = ")  + std::to_string(s)).c_str());
    painter->drawText(15, line += 20, (std::string("p = ")  + std::to_string(p)).c_str());
	painter->drawText(15, line += 20, (std::string("Residual")).c_str());
}

void Window::paintEvent(QPaintEvent * /* event */) {
    QPainter painter(this);
	
	switch (type_of_graph)
	{
		case 0:
			draw_f(&painter);
			break;
		case 1:
			draw_aprx(&painter);
			break;
		case 2:
			draw_res(&painter);
			break;
		default:
			break;
	}
}
