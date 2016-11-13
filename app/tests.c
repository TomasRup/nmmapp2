// From lectures (lt. 'netiktis')
long double incorrectness(
        const long double h, 
        const long double Tau) {

    const long double x = 0.79;
    const long double t = 1.84;
    const long double alpha = 0.001;

    const long double complex incorrectness = ((uAccurate(x, t + Tau) - uAccurate(x, t)) / Tau)
        - (I * (1.0 / 2.0) 
            * (((uAccurate(x + h, t + Tau) - 2 * uAccurate(x, t + Tau) + uAccurate(x - h, t + Tau)) / pow(h, 2)) 
                + (((uAccurate(x + h, t) - 2 * uAccurate(x, t) + uAccurate(x - h, t)) / pow(h, 2)))))
        - (alpha 
            * pow(cabsl((uAccurate(x, t + Tau) + uAccurate(x, t)) / 2.0), 2) 
            * (1.0 / 2.0) 
            * (((uAccurate(x + h, t + Tau) - uAccurate(x - h, t + Tau)) / (2.0 * h)) 
                + ((uAccurate(x + h, t) - uAccurate(x - h, t)) / (2.0 * h))))   
        - ((f(x, t + Tau, alpha) + f(x, t, alpha)) / 2.0);

    return cabsl(incorrectness);
}

// From lecturess (lt. 'netiktis2')
long double incorrectness2(
        const long double h, 
        const long double Tau) {

    // Some random values
    const long double x = 0.2;
    const long double t = 1.85;
    const long double alpha = 2.2;

    // Calculate as many F values as we need
    const int j = x / h;
    const int length = j + 2;

    long double complex *uNextNew = malloc((length + 1) * sizeof(long double complex));
    long double complex *uPrevious = malloc((length + 1) * sizeof(long double complex));
    long double complex *uNextOld = malloc((length + 1) * sizeof(long double complex));
    long double complex *fNext = malloc((length + 1) * sizeof(long double complex));
    long double complex *fPrevious = malloc((length + 1) * sizeof(long double complex));

    for (int i = 0 ; i <= length ; i++) {
        uPrevious[i] = uAccurate(i * h, t);
        uNextOld[i] = uAccurate(i * h, t + Tau);
        fNext[i] = f(i * h, t + Tau, alpha);
        fPrevious[i] = f(i * h, t, alpha);
    }

    initF(uNextNew, uPrevious, uNextOld, fNext, fPrevious, h, alpha, Tau, length);

    // Calculating incorrectnes
    const long double complex incorrectness = uAccurate(x + h, t + Tau)
        - cFromC1(h, Tau) * uAccurate(x, t + Tau)
        + uAccurate(x - h, t + Tau)
        + uNextNew[j];

    free(uNextNew);
    free(uPrevious);
    free(uNextOld);
    free(fNext);
    free(fPrevious);

    return cabsl(incorrectness);
}

// Testing uAccurate(x, t) and f(x, t)
void assertTest1() {

    const long double h = 0.005;
    const long double Tau = 0.0005;

    const long double firstIncorrectness = incorrectness(h, Tau);
    const long double secondIncorrectness = incorrectness(h / 10.0, Tau / 10.0);

    printf("\nTest1: netiktis1(h, Tau) / netiktis1(h / 10, Tau / 10) = %Lf", firstIncorrectness / secondIncorrectness);

    assert(firstIncorrectness / secondIncorrectness > 95);
}

// Testing F and C values
void assertTest2() {

    const long double h = 0.05;
    const long double Tau = 0.0005;

    const long double firstIncorrectness = incorrectness2(h, Tau);
    const long double secondIncorrectness = incorrectness2(h / 10.0, Tau / 10.0);
    
    printf("\nTest2: netiktis2(h, Tau) / netiktis2(h / 10, Tau / 10) = %Lf", firstIncorrectness / secondIncorrectness);

    assert(round(firstIncorrectness / secondIncorrectness) > 9500);
}

// Testing Thomas algorithm (TDMA)
void assertTest3() {

    // Initializing constants
    const long double machinePrecision = 1e-12;
    
    const long double kappa1 = 1;
    const long double kappa2 = 1;
    
    const int N = 1000;
    const long double h = 1.0 / (long double) N;
    const long double Tau = 0.005;
    const long double complex C = cFromC1(h, Tau);

    // Initializing y
    long double complex *y = malloc((N + 1) * sizeof(long double complex));
    
    for (int i = 0 ; i <= N ; i++) {
        y[i] = sin(i) + I * cos(i);
    }
    
    y[0] = kappa1 * y[1];
    y[N] = kappa2 * y[N - 1];

    // Initializing F
    long double complex *F = malloc((N + 1) * sizeof(long double complex));
    F[0] = 0;
    F[N] = 0;
    for (int i = 1 ; i < N ; i++) {
       F[i] = (-1) * (C * y[i] - y[i + 1] - y[i - 1]);
    }

    // Initializing other diagonals
    long double complex *a = malloc((N + 1) * sizeof(long double complex));
    a[N] = -1;
    a[0] = 0;
    for (int i = 1 ; i < N ; i++) {
        a[i] = 1;
    }

    long double complex *b = malloc((N + 1) * sizeof(long double complex));
    b[0] = 1;
    b[N] = 1;
    for (int i = 1 ; i < N ; i++) {
        b[i] = (-1) * C;
    }

    long double complex *c = malloc((N + 1) * sizeof(long double complex));
    c[0] = -1;
    c[N] = 0;
    for (int i = 1 ; i < N ; i++) {
        c[i] = 1;
    }

    // Resolving TDMA
    thomasAlgorithm(a, b, c, F, N);

    // Making sure the biggest difference is not bigger than the machine precision
    long double maxSubtraction = 0.0;

    for (int i = 0 ; i <= N ; i++) {
        const long double subtraction = cabsl(F[i] - y[i]);
        if (subtraction > maxSubtraction) {
            maxSubtraction = subtraction;
        }
    }

    printf("\nTest3: maxSubtraction = %.25Lf", maxSubtraction);

    free(y);
    free(F);
    free(a);
    free(b);
    free(c);

    assert(maxSubtraction < machinePrecision);
}

// Comparing u values to values calculated after Thomas algorithm
void assertTest4() {

    // Maximum allowed difference
    const long double maxAllowedDifference = 0.0001;

    // Mock configuration 1
    struct configuration cfg;
    cfg.N = 1000;
    cfg.Tau = 0.05;
    cfg.T = 0.8;
    cfg.alpha = 0.23;
    cfg.delta = 1e-10;

    const long double h = 1.0 / (long double) cfg.N;

    // Initializing results matrix
    const int amountOfIterations = (int) (cfg.T / cfg.Tau);
	long double complex **finalResultsMatrix = malloc((amountOfIterations + 1) * sizeof(long double complex*));
    for (int i = 0 ; i <= amountOfIterations ; i++) {
        finalResultsMatrix[i] = malloc((cfg.N + 1) * sizeof(*finalResultsMatrix[i]));
	}

    // Getting  solutions
    solve(finalResultsMatrix, cfg);

    // Initializing subtractions arrays
    long double maxSubtraction = 0.0;

    // Getting max delta
    for (int i = 0 ; i <= amountOfIterations ; i++) {
        for (int j = 0 ; j <= cfg.N ; j++) {
            const long double subtraction = cabsl(finalResultsMatrix[i][j] - uAccurate(j * h, i * cfg.Tau));

            if (subtraction > maxSubtraction) {
                maxSubtraction = subtraction;
            }
        }
    }

    // Comparison
    printf("\nTest4: maxSubtraction = %.25Lf", maxSubtraction);
    free(finalResultsMatrix);
    assert(maxSubtraction <= maxAllowedDifference);
}