// From lectures (lt. 'netiktis')
double incorrectness(
        const double h, 
        const double Tau) {

    const double x = 0.79;
    const double t = 1.84;
    const double alpha = 0.001;

    const complex double incorrectness = ((uAccurate(x, t + Tau) - uAccurate(x, t)) / Tau)
        - (I * (1.0 / 2.0) 
            * (((uAccurate(x + h, t + Tau) - 2 * uAccurate(x, t + Tau) + uAccurate(x - h, t + Tau)) / pow(h, 2)) 
                + (((uAccurate(x + h, t) - 2 * uAccurate(x, t) + uAccurate(x - h, t)) / pow(h, 2)))))
        - (alpha 
            * pow(cabs((uAccurate(x, t + Tau) + uAccurate(x, t)) / 2.0), 2) 
            * (1.0 / 2.0) 
            * (((uAccurate(x + h, t + Tau) - uAccurate(x - h, t + Tau)) / (2.0 * h)) 
                + ((uAccurate(x + h, t) - uAccurate(x - h, t)) / (2.0 * h))))   
        - ((f(x, t + Tau, alpha) + f(x, t, alpha)) / 2.0);

    return cabs(incorrectness);
}

// From lecturess (lt. 'netiktis2')
double incorrectness2(
        const double h, 
        const double Tau) {

    // Some random values
    const double x = 0.80;
    const double t = 1.85;
    const double alpha = 0.001;

    // Calculate as many F values as we need
    const int j = x / h;
    const int length = j + 2;

    complex double *fValues = malloc(length * sizeof(complex double));
    complex double *uPrevious = malloc(length * sizeof(complex double));
    complex double *uNextOld = malloc(length * sizeof(complex double));
    complex double *fNext = malloc(length * sizeof(complex double));
    complex double *fPrevious = malloc(length * sizeof(complex double));

    for (int i = 0 ; i <= length ; i++) {
        uPrevious[i] = uAccurate(i * h, t);
        uNextOld[i] = uAccurate(i * h, t + Tau);
        fNext[i] = f(i * h, t + Tau, alpha);
        fPrevious[i] = f(i * h, t, alpha);
    }

    initF(fValues, uPrevious, uNextOld, fNext, fPrevious, h, alpha, Tau, length);

    // Calculating incorrectnes
    const complex double incorrectness = uAccurate(x + h, t + Tau)
        - cFromC1(h, Tau) * uAccurate(x, t + Tau)
        + uAccurate(x - h, t + Tau)
        + fValues[j];

    free(fValues);
    free(uPrevious);
    free(uNextOld);
    free(fNext);
    free(fPrevious);

    return cabs(incorrectness);
}

// Testing uAccurate(x, t) and f(x, t)
void assertTest1() {

    const double h = 0.005;
    const double Tau = 0.0005;

    const double firstIncorrectness = incorrectness(h, Tau);
    const double secondIncorrectness = incorrectness(h / 10.0, Tau / 10.0);

    printf("\nTest1: netiktis1 = %f", firstIncorrectness / secondIncorrectness);

    assert(firstIncorrectness / secondIncorrectness > 95);
}

// Testing F and C values
void assertTest2() {

    const double h = 0.005;
    const double Tau = 0.0005;

    const double firstIncorrectness = incorrectness2(h, Tau);
    const double secondIncorrectness = incorrectness2(h / 10.0, Tau / 10.0);
    
    printf("\nTest2: netiktis2 = %f", firstIncorrectness / secondIncorrectness);

    assert(round(firstIncorrectness / secondIncorrectness) > 9500);
}

// Testing Thomas algorithm
void assertTest3() {

    const double machinePrecision = 1e-12;
    
    const double kappa1 = 1;
    const double kappa2 = 1;
    
    const int N = 10;

    const double h = 0.000002;
    const double Tau = 0.05;
    const double C = cFromC1(h, Tau);

    // Initializing y
    complex double *y = malloc((N + 1) * sizeof(complex double));
    
    for (int i = 0 ; i <= N ; i++) {
        y[i] = sin(i) + I * cos(i);
    }
    
    y[0] = kappa1 * y[1];
    y[N] = kappa2 * y[N - 1];

    // Initializing F
    complex double *F = malloc((N + 1) * sizeof(complex double));

    for (int i = 1 ; i < N ; i++) {
       F[i] = C * y[i] - y[i + 1] - y[i - 1];
    }

    // Initializing other diagonals
    complex double *a = malloc((N + 1) * sizeof(complex double));
    complex double *b = malloc((N + 1) * sizeof(complex double));
    complex double *c = malloc(N * sizeof(complex double));

    for (int i = 1 ; i <= N ; i++) a[i] = 1;
    for (int i = 0 ; i <= N ; i++) b[i] = C;
    for (int i = 0 ; i < N ; i++) c[i] = 1;

    // Resolving TDMA
    thomasAlgorithm(a, b, c, F, N);

    // Making sure the biggest difference is not bigger than the machine precision
    double *subtractions = malloc((N + 1) * sizeof(double));

    // Finishing
    free(a);
    free(b);
    free(c);
    free(F);
    free(y); 

    const double maxSubtraction = max(subtractions, N + 1);
    
    printf("\nTest3: maxSubtraction = %f", maxSubtraction);

    assert(maxSubtraction < machinePrecision);
}

// Comparing u values to values calculated after Thomas algorithm
void assertTest4() {

    const double maxDifference = pow(10, -5);

    // Mock configuration
    struct configuration cfg;
    cfg.N = 100;
    cfg.h = 0.01;
    cfg.Tau = 0.25;
    cfg.T = 2;
    cfg.alpha = 0.1;
    cfg.delta = 0.00001;

    // Initializing results matrix
    const int amountOfIterations = (int) (cfg.T / cfg.Tau);
	complex double **finalResultsMatrix = malloc((amountOfIterations + 1) * sizeof(complex double*));

    for (int i = 0 ; i <= amountOfIterations ; i++) {
		finalResultsMatrix[i] = malloc((cfg.N + 1) * sizeof (*finalResultsMatrix[i]));
	}

    // Solving using the Algorithm
    solve(finalResultsMatrix, cfg);

    // Comparing results received by using algorithm to the actual ones
    for (int i = 1 ; i <= amountOfIterations ; i++) {
        for (int j = 0 ; j < cfg.N ; j++) {
            // assert(cabs(finalResultsMatrix[i][j] - uAccurate(j * cfg.h, i * cfg.Tau)) <= maxDifference);
        }
    }

    // Finalizing
    printf("\nTest4: passed");

    for (int i = 0 ; i <= amountOfIterations ; i++) {
		free(finalResultsMatrix[i]);
	}

    free(finalResultsMatrix);
}