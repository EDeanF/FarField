#include "Header.h"
#include "functions.h"

void MakePlot(cmplxpnt plot[], int plotsize, double start, double stop, complex<double> (*func)(double))
{
	//plot a function

	//starting point
	double x=start;
	double range=abs(stop-start);

	//create plot
	for(int i=0; i<plotsize; ++i){
		plot[i].x=x;
		plot[i].y=(*func)(x);
		x+=range/plotsize;
	}
}

void MakeFarField(cmplxpnt farfield[], int plotsize, double start, double stop, complex<double> (*func)(double, complex<double> (*)(double x)),complex<double>(*func2FT)(double))
{
	//plot the farfield using fourier() for continuous nearfield or sampledFT() for discretized nearfield
	//nearfield given as a function

	//starting point
	double u=start;
	double range=abs(stop-start);

	//create plot
	for(int i=0; i<plotsize; ++i){
		farfield[i].x=u;
		farfield[i].y=(*func)(sin(u)/lambda,(*func2FT));
		u+=range/plotsize;
	}
}

void MakeFarField(cmplxpnt farfield[], int plotsize, double start, double stop, cmplxpnt nearfield[])
{
	//plot the farfield using fourier() for continuous nearfield or sampledFT() for discretized nearfield
	//nearfield given as an array

	//starting point
	double u=start;
	double range=abs(stop-start);

	//create plot
	for(int i=0;i<plotsize;++i){
		farfield[i].x=u;
		farfield[i].y=sampledFT(sin(u)/lambda,nearfield);
		u+=range/plotsize;
	}
}

void PrntPlot( cmplxpnt plot[], int plotsize, string title, string Y)
{
	ofstream file(title+".txt");
	file<<"#X\t"<<Y<<"real\t"<<Y<<"imag\n";		//title
	for(int i=0;i<plotsize; ++i){
		file<< plot[i].x<<"\t"<<plot[i].y.real()<<"\t"<<plot[i].y.imag()<<"\n";
	}
	file.close();
}

string MakeTitle()
{
	string title;
	ostringstream phistream;
	phistream<<n;
	ostringstream periodstring;
	periodstring<<NUMperiod;
	ostringstream samplestring;
	samplestring<<NUMsample;
	title="n"+phistream.str()+"_P"+periodstring.str()+"_S"+samplestring.str();
	return title;
}
