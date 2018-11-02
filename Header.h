#include <iostream>
#include <fstream>
#include <sstream>
#include <complex>
#include <cmath>

using namespace std;

namespace constants
{
	//wave properites
	const double pi=3.14159265358979;
	const complex<double> j(0,1);
	const double c=2.99792458e8f;
	const double f=10e9;
	const double lambda=c/f;
	const double k=2*pi/lambda;
	const double mu=4*pi*1e-7f;
	const double eps=1/mu/c/c;
	const double eta=119.9169832f*pi;
	//incident and deflection angles
	const double m=0, n=3;
	const double theta=pi*m/12;
	const double phi=pi*n/12;
	//additional parameters
	const double period=lambda/abs(sin(phi)-sin(theta));
	//const double period=3*10*lambda;
	const double amp = sqrt(cos(theta)/cos(phi));
	//sampling parameters
	const int NUMperiod=5;
	const int NUMsample=5;
	const int NUMpnts=10000;
	const int NUMtrials=20;
}

struct cmplxpnt{double x; complex<double> y;};
using namespace constants;

