// uAccurate(x, t) - some manually typed function
long double complex uAccurate(
        const long double x, 
        const long double t) {
            
    assert(x >= 0);
    assert(x <= 1);
    
    return (2.0L + I * t) * cos(2.0L * M_PI * x);
}

// ∂u / ∂x
long double complex uDerivativeX(
        const long double x,
        const long double t) {

    assert(x >= 0);
    assert(x <= 1);

    return (-2.0L) * M_PI * (2.0L + I * t) * sin(2.0L * M_PI * x);
}

// ∂u / ∂t
long double complex uDerivativeT(
        const long double x, 
        const long double t) {

    assert(x >= 0);
    assert(x <= 1);

    return I * cos(2.0L * M_PI * x);
}

// ∂^2u / ∂x^2
long double complex uSquareDerivativeXSquare(
        const long double x, 
        const long double t) {

    assert(x >= 0);
    assert(x <= 1);

    return (-4.0L) * pow(M_PI, 2) * (2.0L + I * t) * cos(2.0L * M_PI * x);
}

// f(x, t) from Schrodinger's equation (C)
long double complex f(
        const long double x, 
        const long double t, 
        const long double alpha) {

    assert(x >= 0);
    assert(x <= 1);

    return uDerivativeT(x, t)
        - (I * uSquareDerivativeXSquare(x, t))
        - (alpha * pow(cabsl(uAccurate(x, t)), 2) * uDerivativeX(x, t));
}