// uAccurate(x, t) - some manually typed function
complex double uAccurate(
        const double x, 
        const double t) {
            
    assert(x >= 0);
    assert(x <= 1);
    
    return (2.0 + I * t) * cos(2.0 * M_PI * x);
}

// ∂u / ∂x
complex double uDerivativeX(
        const double x,
        const double t) {

    assert(x >= 0);
    assert(x <= 1);

    return (-2.0) * M_PI * (2.0 + I * t) * sin(2.0 * M_PI * x);
}

// ∂u / ∂t
complex double uDerivativeT(
        const double x, 
        const double t) {

    assert(x >= 0);
    assert(x <= 1);

    return I * cos(2.0 * M_PI * x);
}

// ∂^2u / ∂x^2
complex double uSquareDerivativeXSquare(
        const double x, 
        const double t) {

    assert(x >= 0);
    assert(x <= 1);

    return (-4.0) * pow(M_PI, 2) * (2.0 + I * t) * cos(2.0 * M_PI * x);
}

// f(x, t) from Schrodinger's equation (C)
complex double f(
        const double x, 
        const double t, 
        const double alpha) {

    assert(x >= 0);
    assert(x <= 1);

    return uDerivativeT(x, t)
        - (I * uSquareDerivativeXSquare(x, t))
        - (alpha * pow(cabs(uAccurate(x, t)), 2) * uDerivativeX(x, t));
}