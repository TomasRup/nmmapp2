#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <mpi.h>
#include <complex.h>
#include <math.h>
#include <string.h>
#include <sys/time.h>
#include <assert.h>

#include "configuration.c"
#include "helpers.c"
#include "equations.c"
#include "algorithms.c"
#include "tests.c"

int main(int argc, char *argv[]) {

	// Multi core implementation parameters
	MPI_Init(&argc, &argv);

	int currentProcess;
	MPI_Comm_rank(MPI_COMM_WORLD, &currentProcess);

	int totalProcesses;
	MPI_Comm_size(MPI_COMM_WORLD, &totalProcesses);

	// Initializing stats
    struct timeval appExecutionStarted;
	gettimeofday(&appExecutionStarted, NULL);

	// Run tests
	assertTest1();
	assertTest2();
	assertTest3();
	assertTest4();

	printf("\n\nTests have been completed, beginning initial execution...\n");

	// Reading configuration by processor id
	struct configuration cfg;
	cfg = getConfigurationBy(currentProcess);

	// Initializing results matrix
	const int amountOfIterations = (int) (cfg.T / cfg.Tau);
	long double complex **finalResultsMatrix = malloc((amountOfIterations + 1) * sizeof(long double complex*));

	for (int i = 0 ; i <= amountOfIterations ; i++) {
		finalResultsMatrix[i] = malloc((cfg.N + 1) * sizeof(*finalResultsMatrix[i]));
	}

	// Actually solving equation
	solve(finalResultsMatrix, cfg);

	// Printint out the result
	for (int i = 0 ; i < amountOfIterations ; i++) {
		printListFileForSage(finalResultsMatrix[i], cfg.N, i, currentProcess);
	}

	// Finalizing
	free(finalResultsMatrix);
	MPI_Finalize();

	// Printing stats
	struct timeval appExecutionEnded;
	gettimeofday(&appExecutionEnded, NULL);
	printStats(appExecutionStarted, appExecutionEnded, currentProcess, totalProcesses);

	// Success
	return 0;
}
