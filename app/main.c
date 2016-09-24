#include <stdio.h>
#include "mpi.h"

#include "equations.c"

#define PARAMETER_A 2
#define PARAMETER_C 3
#define PARAMETER_D 4
#define PARAMETER_K 1
#define PARAMETER_ALPHA 1
#define PARAMETER_BETA 2
#define PARAMETER_GAMMA 1
#define PARAMETER_THETA 2

int main(int argc, char *argv[]) {

	MPI_Init(&argc, &argv);

	int totalProcesses;
	MPI_Comm_size(MPI_COMM_WORLD, &totalProcesses);

	int currentProcess;
	MPI_Comm_rank(MPI_COMM_WORLD, &currentProcess);

	printf("current: %d out of %d", currentProcess, totalProcesses);

	MPI_Finalize();
	return 0;
}
