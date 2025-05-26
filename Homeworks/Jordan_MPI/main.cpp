#include "mpi.h"
#include "functions.h"


int main(int argc, char *argv[])
{
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Init(&argc, &argv);


    if (argc < 5)
    {
        fprintf(stderr, "[-] Not enough arguments.\n");
        MPI_Finalize();
        return 1;
    }

    int n = 0, m = 0, r = 0, s = 0;
    if (
        sscanf(argv[1], "%d", &n) != 1 ||
        sscanf(argv[2], "%d", &m) != 1 ||
        sscanf(argv[3], "%d", &r) != 1 ||
        sscanf(argv[4], "%d", &s) != 1)
    {
        printf("[-] Mistake in args!\n");
        MPI_Finalize();
        return 2;
    }

    int rank, world_size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &world_size);

    int p = world_size;
    int l = n % m;
    int k = n / m;
    int sum = 0, reduce_sum = 0;

    if (l > 0 && k + 1 < p)
        p = k + 1;
    if (l == 0 && k < p)
        p = k;

    if (world_size > p)
    {
        fprintf(stderr, "[-] Too many processes.\n");
        MPI_Finalize();
        return 3;
    }

    int max_cols = get_max_cols(n, m, p);

    double *matrix = new double[n * max_cols];
    double *inversed_matrix = new double[n * max_cols];
    double *buffer = new double[m * n];
    double *norm = new double[n];

    if (!matrix || !inversed_matrix || !norm || !buffer)
    {
        sum = 1;
    }

    MPI_Allreduce(&sum, &reduce_sum, 1, MPI_INT, MPI_SUM, comm);
    if (reduce_sum > 0)
    {
        if (rank == 0)
        {
            printf("[-] Not enough memory!\n");
        }
        if (matrix)
            delete[] matrix;
        if (inversed_matrix)
            delete[] inversed_matrix;
        if (buffer)
            delete[] buffer;
        if (norm)
            delete[] norm;
        MPI_Finalize();
        return -1;
    }

    // int error;
    // int i;


    char *filename = nullptr;
    if (s == 0) {
        if (argc < 5) {
            fprintf(stderr, "[-] File name do not defined.\n");
            delete[] matrix;
            delete[] inversed_matrix;
            delete[] buffer;
            delete[] norm;
            MPI_Finalize();
            return 1;
        }
        filename = argv[5];
        if (read_matrix(matrix, n, m, p, rank, filename, buffer, comm)) {
            delete[] matrix;
            delete[] inversed_matrix;
            delete[] buffer;
            delete[] norm;
            MPI_Finalize();
            return -2;
        }

    } else {
        init_matrix(matrix, n, m, p, rank, s);
    }
    MPI_Barrier(comm);
    print_matrix_mpi(matrix, n, m, p, rank, buffer, 4, comm);
    MPI_Barrier(comm);


    unit_matrix_mpi(inversed_matrix, n, m, p, rank);
    MPI_Barrier(comm);
    print_matrix_mpi(inversed_matrix, n, m, p, rank, buffer, 4, comm);
    MPI_Barrier(comm);

    mpi_calculate(
        matrix, inversed_matrix, buffer, n, m, p, rank, comm
    );
    MPI_Barrier(comm);


    if (s == 0) {
        if (argc < 5) {
            fprintf(stderr, "[-] File name do not defined.\n");
            delete[] matrix;
            delete[] inversed_matrix;
            delete[] buffer;
            delete[] norm;
            MPI_Finalize();
            return 1;
        }
        filename = argv[5];
        if (read_matrix(matrix, n, m, p, rank, filename, buffer, comm)) {
            delete[] matrix;
            delete[] inversed_matrix;
            delete[] buffer;
            delete[] norm;
            MPI_Finalize();
            return -2;
        }

    } else {
        init_matrix(matrix, n, m, p, rank, s);
    }
    

    if (rank == 0) printf("---------------- MATRIX ----------------\n");
    print_matrix_mpi(matrix, n, m, p, rank, buffer, 4, comm);
    if (rank == 0) printf("----------------------------------------\n");
    if (rank == 0) printf("------------ INVERSED_MATRIX ------------\n");
    print_matrix_mpi(inversed_matrix, n, m, p, rank, buffer, 4, comm);
    if (rank == 0) printf("-----------------------------------------\n");
    double r1, r2 = 0;

    r1 = residual_calculate_mpi(matrix, inversed_matrix, n, m, p, rank, comm);
    // r2 = residual_calculate_mpi(inversed_matrix, matrix, n, m, p, rank, comm);
    
    if (rank == 0) 
        printf("%s : Task = %d Res1 = %e Res2 = %e T1 = %.2f T2 = %.2f S = %d N = "
            "%d M = %d P = %d\n",
            argv[0], 18, r1, r2, -0., 0., s, n, m, p);


    delete[] matrix;
    delete[] inversed_matrix;
    delete[] buffer;
    delete[] norm;
    MPI_Finalize();
    return 0;
}
