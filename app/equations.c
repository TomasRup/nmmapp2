// u(x, t) - some manually typed function
complex double uAccurate(
        const double x, 
        const double t) {

    return (1 - 2 * x) * (t -  I * cos(t));
}

// ∂u / ∂x
complex double uAccurateDerivativeX(
        const double x,
        const double t) {

    return (-2 * t) + (2 * I * cos(t));
}

// ∂u / ∂t
complex double uAccurateDerivativeT(
        const double x, 
        const double t) {

    return (1 - 2 * x) * (1 + I * sin(t));
}

// ∂^2u / ∂x^2
complex double uAccurateSquareDerivativeXSquare(
        const double x, 
        const double t) {

    return pow((4 * t * x) - (4 * I * x * cos(t)), 2);
}

// f(x, t) from Schrodinger's equation (C)
complex double f(
        const double x, 
        const double t, 
        const double alpha) {

    return uAccurateDerivativeT(x, t)
        - (I * uAccurateSquareDerivativeXSquare(x, t))
        - (alpha * pow(cabs(uAccurate(x, t)), 2) * uAccurateDerivativeX(x, t));
}