//functions.cpp
complex<double> EI(double x);
complex<double> HI(double x);
complex<double> EP(double x);
complex<double> HP(double x);
complex<double> Y(double x);
complex<double> Z(double x);
complex<double> R(double x);
complex<double> T(double x);
complex<double> RR(double x);
complex<double> TT(double x);
complex<double> E1(double x);
complex<double> E2(double x);
complex<double> null(double x);
complex<double> unity(double x);
complex<double> delta(double x,double x1);

//fourier.cpp
complex<double> FT_delta(double u,double a, double T0);
complex<double> fourier(double u, complex<double> (*func2FT)(double x));
complex<double> sampledFT(double u, complex<double> (*func2FT)(double x));
complex<double> sampledFT(double u, cmplxpnt nearfield[]);
complex<double> convolute(double uprime, complex<double> (*func2FT)(double));
