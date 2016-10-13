#define true 1
#define false 0

// C = 2 - ((2ih^2) / Tau)
complex double cFromC1(
        const double h, 
        const double Tau) {

    return 2 - ((2 * I * pow(h, 2)) / Tau);
}

// F = (u(j-1) - 2u(j) + u(j+1))
//      - alpha * i * h^2 * | (uNext(j) + u(j)) / 2 | * (((uNext(j+1) - uNext(j-1)) / 2h) + ((u(j+1) - u(j-1)) / 2h)
//      - i * h^2 * (fNext(j) + f(j))
complex double fFromC1(
        const int j,
        const double t,
        const double tNext,
        const double h, 
        const double alpha) {

    return (uAccurate((j - 1) * h, t) - 2 * uAccurate(j * h, t) + uAccurate((j + 1) * h, t))
        - alpha * I * pow(h, 2) 
            * pow(cabs((uAccurate(j * h, tNext) + uAccurate(j * h, t)) / 2), 2)
            * (((uAccurate((j + 1) * h, tNext) - uAccurate((j - 1) * h, tNext)) / (2 * h)) 
                + ((uAccurate((j + 1) * h, t) - uAccurate((j - 1) * h, t)) / (2 * h)))
        - (I * pow(h, 2) * (f(j * h, tNext, alpha) + f(j * h, t, alpha)));
}

void thomasAlgorithm(
        complex double *oldResults, 
        complex double *results, 
        const double h, 
        const double t, 
        const double tNext,
        const int N, 
        const double alpha,
        const double Tau) {

    // Upper tridiagonal matrix diagonal
    complex double *upperDiagonal = NULL;
    upperDiagonal = malloc((N + 1) * sizeof(complex double));
    for (int i = 0 ; i <= N ; i++) upperDiagonal[i] = 1;

    // Middle tridiagonal matrix diagonal
    complex double *middleDiagonal = NULL;
    middleDiagonal = malloc((N + 1) * sizeof(complex double));
    for (int i = 0 ; i <= N ; i++) middleDiagonal[i] = (-1) * cFromC1(h, Tau);

    // Lower tridiagonal matrix diagonal
    complex double *lowerDiagonal = NULL;
    lowerDiagonal = malloc((N + 1) * sizeof(complex double));
    for (int i = 0 ; i <= N ; i++) lowerDiagonal[i] = 1;

    // Initializing results vector
    for (int i = 1 ; i < N ; i++) results[i] = (-1) * fFromC1(i, t, tNext, h, alpha);
    results[0] = results[1];
    results[N] = results[N - 1];

    // Forward elimination phase
    for (int k = 1 ; k <= N ; k++) {
        const complex double m = 1.0 / (middleDiagonal[k] - (lowerDiagonal[k] * upperDiagonal[k - 1])); 
        upperDiagonal[k] = m * upperDiagonal[k];
        results[k] = (results[k] - lowerDiagonal[k] * results[k - 1]) * m;
    }

    // Backward substitution phase
    for (int k = N - 1 ; k > 0 ; k--) {
        results[k] = results[k] - (results[k + 1] * upperDiagonal[k]);
    }

    // Freeing memory
    free(upperDiagonal);
    free(middleDiagonal);
    free(lowerDiagonal);
}

// Solving using Thomas Algorithm (Tridiagonal Matrix Algorithm). The supported form is x(j-1) - C*x(j) + x(j+1) = -F
void solveNonLinearEquationsSystem(
        const int N, 
        const double h, 
        const double Tau,
        const double T,
        const double alpha,
        const double delta,
        const int currentProcess) { 
    
    // Initializing time moment 0
    complex double *oldResults = NULL;
    oldResults = malloc((N + 1) * sizeof(complex double));
    for (int i = 0 ; i <= N ; i++) oldResults[i] = uAccurate(i * h, 0);

    // Iterating over seconds starting at time moment 1 * Tau
    complex double *results = NULL;
	results = malloc((N + 1) * sizeof(complex double));
    const int maxIterationsPossible = (int) (T / Tau);

    for (int iteration = 1 ; iteration <= maxIterationsPossible ; iteration++) {

        const double t = iteration * Tau;
        const double tNext = (1 + iteration) * Tau;

        int lowerThanDelta = false;
        while (lowerThanDelta == false) {

            // Executing TDMA
            thomasAlgorithm(oldResults, results, h, t, tNext, N, alpha, Tau);
            
            // Checking max|uNew - uOld| < delta
            double *subtractions = NULL;
            subtractions = malloc((N + 1) * sizeof(double));

            for (int i = 0 ; i <= N ; i++) {
                subtractions[i] = fabs(cabs(results[i]) - cabs(oldResults[i]));
            }

            if (max(subtractions, N + 1) < delta) {
                lowerThanDelta = true;
            }

            free(subtractions);

            // If it's not over caching current results
            if (lowerThanDelta == false) {
                memcpy(oldResults, results, (N + 1) * sizeof(complex double));
            }
        }

        // Result output
        printListFileForSage(results, N, iteration, currentProcess);
    }
    
    // Freeing memory
    free(oldResults);
    free(results);
}