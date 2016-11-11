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
    const double x = 0.2;
    const double t = 1.85;
    const double alpha = 2.2;

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

    printf("\nTest1: netiktis1(h, Tau) / netiktis1(h / 10, Tau / 10) = %f", firstIncorrectness / secondIncorrectness);

    assert(firstIncorrectness / secondIncorrectness > 95);
}

// Testing F and C values
void assertTest2() {

    const double h = 0.05;
    const double Tau = 0.0005;

    const double firstIncorrectness = incorrectness2(h, Tau);
    const double secondIncorrectness = incorrectness2(h / 10.0, Tau / 10.0);
    
    printf("\nTest2: netiktis2(h, Tau) / netiktis2(h / 10, Tau / 10) = %f", firstIncorrectness / secondIncorrectness);

    assert(round(firstIncorrectness / secondIncorrectness) > 9500);
}

// Testing Thomas algorithm (TDMA)
void assertTest3() {

    // Initializing constants
    const double machinePrecision = 1e-12;
    
    const double kappa1 = 1;
    const double kappa2 = 1;
    
    const int N = 5;
    const double h = 1.0 / (double) N;
    const double Tau = 0.005;
    const complex double C = cFromC1(h, Tau);

    // Initializing y
    complex double *y = malloc((N + 1) * sizeof(complex double));
    
    for (int i = 0 ; i <= N ; i++) {
        y[i] = sin(i) + I * cos(i);
    }
    
    y[0] = kappa1 * y[1];
    y[N] = kappa2 * y[N - 1];

    // Initializing F
    complex double *F = malloc((N + 1) * sizeof(complex double));
    F[0] = 0;
    F[N] = 0;
    for (int i = 1 ; i < N ; i++) {
       F[i] = (-1) * (C * y[i] - y[i + 1] - y[i - 1]);
    }

    // Initializing other diagonals
    complex double *a = malloc((N + 1) * sizeof(complex double));
    a[N] = -1;
    a[0] = 0;
    for (int i = 1 ; i < N ; i++) {
        a[i] = 1;
    }

    complex double *b = malloc((N + 1) * sizeof(complex double));
    b[0] = 1;
    b[N] = 1;
    for (int i = 1 ; i < N ; i++) {
        b[i] = (-1) * C;
    }

    complex double *c = malloc((N + 1) * sizeof(complex double));
    c[0] = -1;
    c[N] = 0;
    for (int i = 1 ; i < N ; i++) {
        c[i] = 1;
    }

    // Resolving TDMA
    thomasAlgorithm(a, b, c, F, N);

    // Making sure the biggest difference is not bigger than the machine precision
    double maxSubtraction = 0.0;

    for (int i = 0 ; i <= N ; i++) {
        const double subtraction = cabs(F[i] - y[i]);
        if (subtraction > maxSubtraction) {
            maxSubtraction = subtraction;
        }
    }

    // Finishing
    free(a);
    free(b);
    free(c);
    free(F);
    free(y);

    printf("\nTest3: maxSubtraction = %.20f", maxSubtraction);

    assert(maxSubtraction < machinePrecision);
}

// Comparing u values to values calculated after Thomas algorithm
void assertTest4() {

    // Maximum allowed difference
    const double maxAllowedDifference = 0.0001;

    // Mock configuration 1
    struct configuration cfg;
    cfg.N = 1000;
    cfg.Tau = 0.05;
    cfg.T = 2;
    cfg.alpha = 0.033;
    cfg.delta = 0.001;

    const double h = 1.0 / (double) cfg.N;

    // Initializing results matrix
    const int amountOfIterations = (int) (cfg.T / cfg.Tau);
	complex double **finalResultsMatrix = malloc((amountOfIterations + 1) * sizeof(complex double*));
    for (int i = 0 ; i <= amountOfIterations ; i++) {
        finalResultsMatrix[i] = malloc((cfg.N + 1) * sizeof(*finalResultsMatrix[i]));
	}

    // Getting  solutions
    solve(finalResultsMatrix, cfg);

    // Initializing subtractions arrays
    double maxSubtraction = 0.0;

    // Getting max delta
    for (int i = 0 ; i <= amountOfIterations ; i++) {
        for (int j = 0 ; j <= cfg.N ; j++) {
            const double subtraction = cabs(finalResultsMatrix[i][j] - uAccurate(j * h, i * cfg.Tau));

            if (subtraction > maxSubtraction) {
                maxSubtraction = subtraction;
            }
        }
    }

    // Comparison
    printf("\nTest4: maxSubtraction = %.20f", maxSubtraction);

    assert(maxSubtraction <= maxAllowedDifference);

    // Finalizing
    free(finalResultsMatrix);
}