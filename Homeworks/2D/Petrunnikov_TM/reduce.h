int init_reduce_sum (int p);
void delete_reduce_sum ();
void sync_max    (int p, double* a = nullptr, int n = 0);
void synchronize (int p, double* a = nullptr, int n = 0);
double reduce_sum_det (int p, int pi, double s);
