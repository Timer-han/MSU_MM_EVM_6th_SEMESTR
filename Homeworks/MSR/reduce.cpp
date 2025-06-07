#include <pthread.h>  
#include <sys/resource.h>
#include <sys/sysinfo.h>

#include "reduce.h"

static double* result = nullptr;

void synchronize(int p, double* a, int n) {
    static pthread_mutex_t m = PTHREAD_MUTEX_INITIALIZER;
    static pthread_cond_t c_in = PTHREAD_COND_INITIALIZER;
    static pthread_cond_t c_out = PTHREAD_COND_INITIALIZER;
    static int t_in = 0;
    static int t_out = 0;
    static double* r = nullptr;
    int i;
    if (p <= 1) return;
    pthread_mutex_lock(&m);
    if (r == nullptr) r = a;
    else for(i = 0; i < n; ++i) r[i] += a[i];
    ++t_in;
    if (t_in >= p) {
        t_out = 0;
        pthread_cond_broadcast(&c_in);
    } else {
        while (t_in < p) {
            pthread_cond_wait(&c_in, &m);
        }
    }
    if (r != a) {
        for (i = 0; i < n; ++i) a[i] = r[i];
    }
    ++t_out;
    if (t_out >= p) {
        t_in = 0;
        r = nullptr;
        pthread_cond_broadcast(&c_out);
    } else {
        while (t_out < p) {
            pthread_cond_wait(&c_out, &m);
        }
    }
    pthread_mutex_unlock(&m);
}

int init_reduce_sum(int p) {
	result = new double[p];
	if (result == nullptr)
		return -1;
	return 0;
}

void delete_reduce_sum() {
	delete[] result;
}

double reduce_sum_det(int p, int k, double s) {
	double sum = 0;
	int l;
	result[k] = s;
	synchronize(p);
	for(l = 0; l < p; ++l) {
		sum += result[l];
	}
	synchronize(p);
	return sum;
}

void sync_max(int p, double* a, int n) {
    static pthread_mutex_t m = PTHREAD_MUTEX_INITIALIZER;
    static pthread_cond_t c_in = PTHREAD_COND_INITIALIZER;
    static pthread_cond_t c_out = PTHREAD_COND_INITIALIZER;
    static int t_in = 0;
    static int t_out = 0;
    static double* r = nullptr;
    int i;
    if (p <= 1) return;
    pthread_mutex_lock(&m);
    if (r == nullptr) r = a;
    else for(i = 0; i < n; ++i) r[i] = (r[i] > a[i] ? r[i] : a[i]);
    ++t_in;
    if (t_in >= p) {
        t_out = 0;
        pthread_cond_broadcast(&c_in);
    } else {
        while (t_in < p) {
            pthread_cond_wait(&c_in, &m);
        }
    }
    if (r != a) {
        for (i = 0; i < n; ++i) a[i] = r[i];
    }
    ++t_out;
    if (t_out >= p) {
        t_in = 0;
        r = nullptr;
        pthread_cond_broadcast(&c_out);
    } else {
        while (t_out < p) {
            pthread_cond_wait(&c_out, &m);
        }
    }
    pthread_mutex_unlock(&m);
}
