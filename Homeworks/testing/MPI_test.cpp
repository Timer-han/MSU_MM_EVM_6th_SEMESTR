#include "/usr/local/intel/impi/5.1.3.223/intel64/include/mpi.h"

#include <stdio.h>

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // std::cout << "Привет от процесса " << rank << " из " << size << " процессов!" << std::endl;
    printf("%d %d\n", rank, size);

    MPI_Finalize();
    return 0;
}
