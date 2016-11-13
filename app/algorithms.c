#define true 1
#define false 0

// C
long double complex cFromC1(
        const long double h, 
        const long double Tau) {
            
    return 2 - ((I * 2.0 * pow(h, 2)) / Tau);
}

// Executing function F from C1 algorithm
void initF(
        long double complex *result,
        long double complex *uPrevious,
        long double complex *uNextOld,
        long double complex *fNext,
        long double complex *fPrevious,
        const long double h,
        const long double alpha,
        const long double Tau,
        const int length) {

    for (int j = 1 ; j < length ; j++) {
        result[j] = (uPrevious[j + 1]
            - 2.0 * uPrevious[j]
            + uPrevious[j - 1]
            - I
                * pow(h, 2)
                * alpha
                * pow(
                    cabsl(
                        (uNextOld[j] + uPrevious[j])
                            / 2.0),
                    2)
                * (((uNextOld[j + 1] - uNextOld[j - 1]) 
                        / (2 * h)) 
                    + ((uPrevious[j + 1] - uPrevious[j - 1]) 
                        / (2 * h)))
            - I
                * pow(h, 2)
                * (fNext[j] + fPrevious[j])
            - (uPrevious[j]
                * ((I * 2 * pow(h, 2)) 
                    / Tau)));
    }

    result[0] = 0;
    result[length] = 0;
}

// |bc0|   |d|
// |abc| = |d| 
// |0ab|   |d|
void thomasAlgorithm(
        long double complex *a,
        long double complex *b,
        long double complex *c,
        long double complex *d,
        int n) {

    c[0] /= b[0];
    d[0] /= b[0];

    for (int i = 1 ; i < n ; i++) {
        c[i] /= b[i] - a[i] * c[i - 1];
        d[i] = (d[i] - a[i] * d[i - 1]) / (b[i] - a[i] * c[i - 1]);
    }

    d[n] = (d[n] - a[n] * d[n - 1]) / (b[n] - a[n] * c[n - 1]);

    for (int i = n - 1 ; i >= 0 ; i--) {
        d[i] -= c[i] * d[i + 1];
    }
}

// Solving using Thomas Algorithm (Tridiagonal Matrix Algorithm). 
// The supported form is x(j-1) - C * x(j) + x(j+1) = -F
void solve(
        long double complex **finalResultsMatrix, 
        struct configuration cfg) { 
    
    // Initializing constants
    const long double h = 1.0 / (long double) cfg.N;

    // Initializing arrays
    long double complex *uPrevious = malloc((cfg.N + 1) * sizeof(long double complex));
    long double complex *uNextNew = malloc((cfg.N + 1) * sizeof(long double complex));
    long double complex *uNextOld = malloc((cfg.N + 1) * sizeof(long double complex));

    long double complex *fNext = malloc((cfg.N + 1) * sizeof(long double complex));
    long double complex *fPrevious = malloc((cfg.N + 1) * sizeof(long double complex));

    // 'u' at time moment 0
    for (int i = 0 ; i <= cfg.N ; i++) {
        uPrevious[i] = uAccurate(i * h, 0);
    }

    // 'u' next 'old' at time moment 0
    memcpy(uNextOld, uPrevious, (cfg.N + 1) * sizeof(long double complex));

    // The first result array is the one from time 0
    memcpy(finalResultsMatrix[0], uNextOld, (cfg.N + 1) * sizeof(long double complex));

    // 'f' at time moment 0
    for (int i = 0 ; i <= cfg.N ; i++) {
        fPrevious[i] = f(i * h, 0, cfg.alpha);
    }

    // Iterating over seconds starting at time moment 1 * Tau
    for (int timeIteration = 1 ; timeIteration <= (int) (cfg.T / cfg.Tau) ; timeIteration++) {

        const long double t = timeIteration * cfg.Tau;

        // f function values at times t and t + Tau
        for (int i = 0 ; i <= cfg.N ; i++) {
            fNext[i] = f(i * h, t, cfg.alpha);
        }
        
        // Loop until max|uNextNew - uNextOld| < delta
        int lowerThanDelta = false;
        while (lowerThanDelta == false) {

            // Calculating results (F) vector
            initF(uNextNew, uPrevious, uNextOld, fNext, fPrevious, h, cfg.alpha, cfg.Tau, cfg.N);
            
            // We need to negate F vector for solving of tridiagonal matrix
            for (int i = 1 ; i < cfg.N ; i++) {
                uNextNew[i] = -uNextNew[i];
            }

            // Initializing Matrix diagonal values
            long double complex *lowerDiagonal = malloc((cfg.N + 1) * sizeof(long double complex));
            lowerDiagonal[0] = 0;
            lowerDiagonal[cfg.N] = -1;
            for (int i = 1 ; i < cfg.N ; i++) {
                lowerDiagonal[i] = 1;
            }

            long double complex *middleDiagonal = malloc((cfg.N + 1) * sizeof(long double complex));
            middleDiagonal[0] = 1;
            middleDiagonal[cfg.N] = 1; 
            const long double complex C = cFromC1(h, cfg.Tau);
            for (int i = 1 ; i < cfg.N ; i++) {
                middleDiagonal[i] = -C;
            }

            long double complex *upperDiagonal = malloc((cfg.N + 1) * sizeof(long double complex));
            upperDiagonal[0] = -1;
            upperDiagonal[cfg.N] = 0;
            for (int i = 1 ; i < cfg.N ; i++) {
                upperDiagonal[i] = 1;
            }

            // Executing TDMA
            thomasAlgorithm(lowerDiagonal, middleDiagonal, upperDiagonal, uNextNew, cfg.N);

            // Checking max|uNextNew - uNextOld| < delta
            long double maxSubtraction = 0.0;
            for (int i = 0 ; i <= cfg.N ; i++) {
                const long double subtraction = cabsl(uNextNew[i] - uNextOld[i]);

                if (subtraction > maxSubtraction) {
                    maxSubtraction = subtraction;
                }
            }

            if (maxSubtraction < cfg.delta) {
                lowerThanDelta = true;
            }

            // The uNext new becomes old
            memcpy(uNextOld, uNextNew, (cfg.N + 1) * sizeof(long double complex));

            free(lowerDiagonal);
            free(middleDiagonal);
            free(upperDiagonal);
        }

        // Adding result
        memcpy(finalResultsMatrix[timeIteration], uNextNew, (cfg.N + 1) * sizeof(long double complex));

        // Swapping previous u values with next new values
        long double complex *tempList = malloc((cfg.N + 1) * sizeof(long double complex));
        
        memcpy(tempList, uNextNew, (cfg.N + 1) * sizeof(long double complex));
        memcpy(uNextNew, uPrevious, (cfg.N + 1) * sizeof(long double complex));
        memcpy(uPrevious, tempList, (cfg.N + 1) * sizeof(long double complex));
        
        free(tempList);

        // Current 'f' becomes previous
        memcpy(fPrevious, fNext, (cfg.N + 1) * sizeof(long double complex));
    }
}