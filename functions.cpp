#include "Header.h"

complex<double> EI(double x)
{
	complex<double> EI=-cos(theta)*exp(-j*k*sin(theta)*x);
	return EI;
}

complex<double> HI(double x)
{
	complex<double> HI=-1/eta*exp(-j*k*sin(theta)*x);
	return HI;
}

complex<double> EP(double x)
{
	complex<double> EP=-amp*cos(phi)*exp(-j*k*sin(phi)*x);
	return EP;
}

complex <double> HP(double x)
{
	complex<double> HP=-amp/eta*exp(-j*k*sin(phi)*x);
	return HP;
}

complex<double> Y(double x)
{
	complex<double> Y=2e0*(HI(x)-HP(x))/(EI(x)+EP(x));
	return Y;
}

complex<double> Z(double x)
{
	complex<double> Z=2e0*(EI(x)-EP(x))/(HI(x)+HP(x));
	return Z;
}

complex<double> R(double x)
{
	complex<double> R=-0.5/cos(theta)*(Y(x)*eta*pow(cos(theta),2)-Z(x)/eta)/(1e0+0.25*Z(x)*Y(x)+0.5*(Y(x)*eta*pow(cos(theta),2)+Z(x)/eta));
	return R;
}

complex<double> T(double x)
{
	complex<double> T=(1.e0-0.25*Z(x)*Y(x))/(1e0+0.25*Z(x)*Y(x)+0.5*(Y(x)*eta*pow(cos(theta),2)+Z(x)/eta));
	return T;
}

complex<double> RR(double x){
	complex<double> RR=cos(phi)/cos(theta)*pow(R(x),2);
	return RR;
}

complex<double> TT(double x){
	complex<double> TT=1e0-RR(x);
	return TT;
}

complex<double> E1(double x){
	complex<double> E1=EI(x)*(1e0+R(x));
	return E1;
}

complex<double> E2(double x)
{
	complex<double> E2=EI(x)*T(x);
	return E2;
}

complex<double> null(double x)
{
	complex<double> null=0;
	return null;
}

complex<double> delta(double x,double x1)
{
	complex<double> delta;
	if (x==x1)
		delta=1;
	else 
		delta=0;
	return delta;
}

complex<double> unity(double x)
{
	complex<double> unity=1e0;
	return unity;
}

