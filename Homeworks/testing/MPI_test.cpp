#include <mpi.h>
#include <stdio.h>

int main(int argc, char** argv) {
    // Initialize the MPI environment
    MPI_Init(NULL, NULL);

    // Get the number of processes
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    // Get the name of the processor
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);

    double * data = new double[100];

    for (int i = 0; i < 100; i++) {
        data[i] = world_rank + i;
    }

    if (world_rank == 1) {
        MPI_Bcast(data, 100, MPI_DOUBLE, world_rank, MPI_COMM_WORLD);
    }
    else {
        MPI_Bcast(data, 100, MPI_DOUBLE, 1, MPI_COMM_WORLD);
    }



    // Print off a hello world message
    printf("Hello world from processor %s, rank %d out of %d processors\n"
           "data[0]: %f\n",
           processor_name, world_rank, world_size, data[0]);
    

    // Finalize the MPI environment.
    MPI_Finalize();
}