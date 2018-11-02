#include "Header.h"

complex<double> delta(double x,double x1);

complex<double> FT_delta(double u,double a, double T0)
{
	double T=a+T0;
	complex<double> FT_delta=exp(j*2e0*pi*u*T);
	return FT_delta;
}

complex<double> fourier(double u, complex<double> (*func2FT)(double x))
{
	complex<double> sum=0;
	int NUMsum=500;
	double lwrlim=0;
	double upprlim=period;
	double dx=abs(upprlim-lwrlim)/NUMsum;
	double x=lwrlim+dx/2;
	for (int i=0;i<NUMsum;++i)
	{
		sum+=func2FT(x)*exp(j*2e0*pi*u*x)*dx;
		x+=dx;
	}
	return sum;
}

complex<double> sampledFT(double u, complex<double> (*func2FT)(double x))
{
	complex<double> sum=0;
	
	//double x=period/NUMsample/2;
	double x=0;
	for(int i=0;i<NUMsample;++i)
	{
		sum+=func2FT(x)*exp(j*2e0*pi*u*x);
		x+=period/NUMsample;
	}
	return sum;
}

complex<double> sampledFT(double u, cmplxpnt nearfield[])
{
	complex<double> sum=0;

	for(int i=0;i<NUMsample;++i)
	{
		sum+=nearfield[i].y*exp(j*2e0*pi*u*nearfield[i].x);
	}
	return sum;
}


complex<double> convolute(double uprime, complex<double> (*func2FT)(double))
{
	complex<double> sum=0;
	int NUMsum=100;

	double lwrlim=-1/lambda, upprlim=1/lambda;
	double du=abs(upprlim-lwrlim)/NUMsum;
	double u=lwrlim+du/2;

	cmplxpnt* FT=new cmplxpnt[NUMsum];

	//create FT
	for(int i=0; i<NUMsum; ++i)
	{
		FT[i].x=u;
		FT[i].y=fourier(u,(*func2FT));
		u+=du;
	}

	//parameters for FT_delta
	double a=period/NUMsample;
	double T0=0;

	//convolution
	for (int ii=0;ii<NUMsample;++ii)
	{
		u=lwrlim+du/2;
		for (int i=0;i<NUMsum;++i)
		{
			sum+=FT[i].y*FT_delta(uprime-u,ii*a,T0)*du;
			u+=du;
		}
	}
	delete[] FT;
	return sum;
}
