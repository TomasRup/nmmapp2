#define true 1
#define false 0

// C
complex double cFromC1(
        const double h, 
        const double Tau) {

    return 2 - ((I * 2 * pow(h, 2)) / Tau);
}

// Executing function F from C1 algorithm
void initF(
        complex double *result,
        complex double *uPrevious,
        complex double *uNextOld,
        complex double *fNext,
        complex double *fPrevious,
        const double h,
        const double alpha,
        const double Tau,
        const int length) {

    for (int j = 1 ; j < length ; j++) {
        result[j] = (uPrevious[j + 1]
            - 2 * uPrevious[j]
            + uPrevious[j - 1]
            - I
                * pow(h, 2)
                * alpha
                * pow(
                    cabs(
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

    result[0] = result[1];
    result[length] = result[length - 1];
}

// |bc0|   |d|
// |abc| = |d| 
// |0ab|   |d|
void thomasAlgorithm(
        complex double *a,
        complex double *b,
        complex double *c,
        complex double *d,
        const int length) {

    c[0] /= b[0];
    d[0] /= b[0];

    for (int i = 1; i < length ; i++) {
        c[i] /= b[i] - a[i] * c[i - 1];
        d[i] = (d[i] - a[i] * d[i - 1]) / (b[i] - a[i] * c[i - 1]);
    }

    d[length] = (d[length] - a[length] * d[length - 1]) / (b[length] - a[length] * c[length - 1]);

    for (int i = length - 1 ; i >= 0 ; i--) {
        d[i] -= c[i] * d[i + 1];
    }
}

// Solving using Thomas Algorithm (Tridiagonal Matrix Algorithm). 
// The supported form is x(j-1) - C * x(j) + x(j+1) = -F
void solve(
        complex double **finalResultsMatrix, 
        struct configuration cfg) { 
    
    // Initializing arrays
    complex double *lowerDiagonal = malloc((cfg.N + 1) * sizeof(complex double));
    complex double *middleDiagonal = malloc((cfg.N + 1) * sizeof(complex double));
    complex double *upperDiagonal = malloc((cfg.N + 1) * sizeof(complex double));
    
    complex double *uPrevious = malloc((cfg.N + 1) * sizeof(complex double));
    complex double *uNextNew = malloc((cfg.N + 1) * sizeof(complex double));
    complex double *uNextOld = malloc((cfg.N + 1) * sizeof(complex double));

    complex double *fNext = malloc((cfg.N + 1) * sizeof(complex double));
    complex double *fPrevious = malloc((cfg.N + 1) * sizeof(complex double));

    // Matrix diagonal values
    lowerDiagonal[0] = 0;
    for (int i = 1 ; i <= cfg.N ; i++) {
        lowerDiagonal[i] = 1;
    }
    
    for (int i = 0 ; i <= cfg.N ; i++) {
        middleDiagonal[i] = (-1) * cFromC1(cfg.h, cfg.Tau);
    }

    upperDiagonal[cfg.N] = 0;    
    for (int i = 0 ; i < cfg.N ; i++) {
        upperDiagonal[i] = 1;
    }

    // 'u' at time moment 0
    for (int i = 1 ; i < cfg.N ; i++) {
        uPrevious[i] = uAccurate(i * cfg.h, 0);
    }

    uPrevious[0] = uPrevious[1];
    uPrevious[cfg.N] = uPrevious[cfg.N - 1];

    // 'u' next 'old' at time moment 0
    memcpy(uNextOld, uPrevious, (cfg.N + 1) * sizeof(complex double));

    // 'f' at time moment 0
    for (int i = 0 ; i <= cfg.N ; i++) {
        fPrevious[i] = f(i * cfg.h, 0, cfg.alpha);
    }

    // Iterating over seconds starting at time moment 1 * Tau
    for (int timeIteration = 1 ; timeIteration <= (int) (cfg.T / cfg.Tau) ; timeIteration++) {

        const double t = timeIteration * cfg.Tau;
        
        // f function values at times t and t + Tau
        for (int i = 0 ; i <= cfg.N ; i++) {
            fNext[i] = f(i * cfg.h, t, cfg.alpha);
        }
        
        // Loop until max|uNextNew - uNextOld| < delta
        int lowerThanDelta = false;
        while (lowerThanDelta == false) {

            // Calculating results (F) vector
            initF(uNextNew, uPrevious, uNextOld, fNext, fPrevious, cfg.h, cfg.alpha, cfg.Tau, cfg.N);
            
            // We need to negate F vector for solving of tridiagonal matrix
            for (int i = 0 ; i <= cfg.N ; i++) {
                uNextNew[i] = (-1) * uNextNew[i];
            }
            
            // Executing TDMA
            thomasAlgorithm(lowerDiagonal, middleDiagonal, upperDiagonal, uNextNew, cfg.N);

            // Checking max|uNextNew - uNextOld| < delta
            double *subtractions = malloc((cfg.N + 1) * sizeof(double));
            for (int i = 0 ; i <= cfg.N ; i++) {
                subtractions[i] = cabs(uNextNew[i] - uNextOld[i]);
            }

            if (max(subtractions, cfg.N + 1) < cfg.delta) {
                lowerThanDelta = true;
            }

            free(subtractions);

            // The uNext new becomes old
            memcpy(uNextOld, uNextNew, (cfg.N + 1) * sizeof(complex double));
        }

        // Adding result
        memcpy(finalResultsMatrix[timeIteration - 1], uNextNew, (cfg.N + 1) * sizeof(complex double));

        // Swapping current u values with new values
        complex double *tempList = malloc((cfg.N + 1) * sizeof(complex double));
        memcpy(tempList, uNextNew, (cfg.N + 1) * sizeof(complex double));
        memcpy(uNextNew, uPrevious, (cfg.N + 1) * sizeof(complex double));
        memcpy(uPrevious, tempList, (cfg.N + 1) * sizeof(complex double));
        free(tempList);

        // Current 'f' becomes previous
        memcpy(fPrevious, fNext, (cfg.N + 1) * sizeof(complex double));
    }
    
    // Freeing memory
    free(lowerDiagonal);
    free(middleDiagonal);
    free(upperDiagonal);

    free(uNextNew);
    free(uPrevious);
    free(uNextOld);

    free(fNext);
    free(fPrevious);
}