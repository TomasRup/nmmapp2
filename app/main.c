#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <mpi.h>
#include <complex.h>
#include <math.h>
#include <string.h>
#include <sys/time.h>

#include "configuration.c"
#include "helpers.c"
#include "equations.c"
#include "algorithms.c"

int main(int argc, char *argv[]) {

    struct timeval appExecutionStarted;
	gettimeofday(&appExecutionStarted, NULL);

	// Multi core implementation parameters
	MPI_Init(&argc, &argv);

	int currentProcess;
	MPI_Comm_rank(MPI_COMM_WORLD, &currentProcess);

	int totalProcesses;
	MPI_Comm_size(MPI_COMM_WORLD, &totalProcesses);

	// Reading configuration by processor id
	struct configuration cfg;
	cfg = getConfigurationBy(currentProcess);

	// Executing calculations
	solveNonLinearEquationsSystem(cfg.N, cfg.h, cfg.Tau, cfg.T, cfg.alpha, cfg.delta, currentProcess);

	// Finalizing MPI process
	MPI_Finalize();

	// Stats
	struct timeval appExecutionEnded;
	gettimeofday(&appExecutionEnded, NULL);

	printStats(appExecutionStarted, appExecutionEnded, currentProcess, totalProcesses);

	return 0;
}
