#include "Header.h"
#include "functions.h"
#include "arrays.h"

//plot and printing
void MakePlot(cmplxpnt plot[], int plotsize, double start, double stop, complex<double> (*func)(double));
void MakeFarField(cmplxpnt plot[], int plotsize, double start, double stop, complex<double> (*func)(double, complex<double> (*)(double x)),complex<double>(*func2FT)(double));
void MakeFarField(cmplxpnt farfield[], int plotsize, double start, double stop, cmplxpnt nearfield[]);
void PrntPlot( cmplxpnt plot[], int plotsize, string title, string Y);

void MakeIdealYZ(cmplxpnt Yplot[],cmplxpnt Zplot[])
{
	MakePlot(Yplot,NUMsample,0,period,Y);
	MakePlot(Zplot,NUMsample,0,period,Z);

}

void MakeLimitYZ(cmplxpnt Yplot[],cmplxpnt Zplot[])
{
	MakeIdealYZ(Yplot,Zplot);

	//apply limitations
	normalizeY(Yplot, NUMsample);
	normalizeZ(Zplot, NUMsample);

	limitYZ(Yplot, NUMsample);
	limitYZ(Zplot, NUMsample);

	revertY(Yplot, NUMsample);
	revertZ(Zplot, NUMsample);
}

void YZ2FarField(cmplxpnt farfield[],cmplxpnt Yplot[],cmplxpnt Zplot[])
{
	//Make Rplot and Tplot

	cmplxpnt Rplot[NUMsample]; cmplxpnt Tplot[NUMsample]; 

	MakePlot(Rplot, NUMsample,0,period,null);
	MakePlot(Tplot, NUMsample,0,period,null);

	MakeRplot(Rplot, NUMsample, Zplot, Yplot);
	MakeTplot(Tplot, NUMsample, Zplot, Yplot);

	//Make E2
		//Make EIplot
		//required for E2plot
	cmplxpnt EIplot[NUMsample];
	MakePlot(EIplot,NUMsample,0,period,EI);

	cmplxpnt E2plot[NUMsample]; 
	MakePlot(E2plot, NUMsample,0,period,null);
	MakeE2plot(E2plot, NUMsample,EIplot,Tplot);

	//PrntPlot(E2plot,NUMsample,"E2_test","E");

	//Calculate Farfield

	cmplxpnt temp[NUMpnts];
	MakePlot(farfield, NUMpnts, -pi/2, pi/2, null);

	MakeFarField(temp, NUMpnts, -pi/2, pi/2, E2plot);
	MakePlot(farfield, NUMpnts, -pi/2, pi/2, null);

	for(int i=0;i<NUMperiod;++i)
	{
		for(int ii=0;ii<NUMpnts;++ii)
		{
			farfield[ii].y+=temp[ii].y*FT_delta(sin(temp[ii].x)/lambda,i*period,-period*NUMperiod/2);
		}
	}
}

void FindLimits(int limits[2], cmplxpnt plot[], int plotsize)
{
	//find peak
	int peak = floor(plotsize/2*(1+2*phi/pi));

	//double check peak; this is necessary for small s
	//for(int find=(peak-10);find<(peak+10);++find)
	for(int find=0;find<plotsize;++find)
	{
		if(plot[find].y.real()>plot[peak].y.real())
			{peak=find;}
	}
		//find lowerlimit
	int i=peak;
	while((plot[i].y.real()>plot[i-1].y.real())&&(i!=0))
		{--i;}
	limits[0]=i;
	//find upperlimit
	i=peak;
	while((plot[i].y.real()>plot[i+1].y.real())&&(i!=(plotsize-1)))
		{++i;}
	limits[1]=i;
}

double integrate(cmplxpnt plot[], int lwrlim, int upprlim)
{
	double sum=0;
	//left Riemann Sum
	double delta = abs((plot[1].x-plot[0].x));
	for(int i=lwrlim; i<upprlim; ++i)
		{sum+=plot[i].y.real()*delta;}

	return sum;
}

double directivity(cmplxpnt plot[], int plotsize)
{
	double directivity;
	
	int limits[2];
	FindLimits(limits,plot,plotsize);

	directivity=integrate(plot,limits[0],limits[1])/integrate(plot,0,plotsize);

	return directivity;
}

double directivity(cmplxpnt plot[], int plotsize, int limits[2])
{
	double directivity;
	directivity=integrate(plot,limits[0],limits[1])/integrate(plot,0,plotsize);
	return directivity;
}

void copypaste(cmplxpnt paste[],cmplxpnt copy[], int plotsize)
{
	for(int i=0;i<plotsize;++i)
		{paste[i]=copy[i];}
};

void optimize(cmplxpnt farfield[NUMpnts],cmplxpnt Yplot[NUMsample],cmplxpnt Zplot[NUMsample],int limits[2])
{
	//goes through each sampled point and applies a +/-10% change
	//seeks to increase intensity in the region specified by limits

	cmplxpnt plusYplot[NUMsample];
	cmplxpnt minusYplot[NUMsample];
	cmplxpnt plusZplot[NUMsample];
	cmplxpnt minusZplot[NUMsample];

	cmplxpnt tempfarfield[NUMpnts];
	
	//for each sample
	for(int s=0;s<NUMsample;++s)
	{
		//initialize
		copypaste(plusYplot,Yplot,NUMsample);
		copypaste(minusYplot,Yplot,NUMsample);
		copypaste(plusZplot,Zplot,NUMsample);
		copypaste(minusZplot,Zplot,NUMsample);

		if(!(abs(Yplot[s].y.imag())>=3/eta))
		{
			plusYplot[s].y.imag(1.1*Yplot[s].y.imag());
			normalizeY(plusYplot,NUMsample);
			limitYZ(plusYplot,NUMsample);
			revertY(plusYplot,NUMsample);
		}
		if(!(abs(Yplot[s].y.imag())<=1/(3*eta)))
		{
			minusYplot[s].y.imag(0.9*Yplot[s].y.imag());
			normalizeY(minusYplot,NUMsample);
			limitYZ(minusYplot,NUMsample);
			revertY(minusYplot,NUMsample);
		}
		if(!(abs(Zplot[s].y.imag())>=3*eta))
		{
			plusZplot[s].y.imag(1.1*Zplot[s].y.imag());
			normalizeZ(plusZplot,NUMsample);
			limitYZ(plusZplot,NUMsample);
			revertZ(plusZplot,NUMsample);
		}
		if(!(abs(Zplot[s].y.imag())<=eta/3))
		{
			minusZplot[s].y.imag(0.9*Zplot[s].y.imag());
			normalizeZ(minusZplot,NUMsample);
			limitYZ(minusZplot,NUMsample);
			revertZ(minusZplot,NUMsample);
		}

		//print directivity
		cout<<directivity(farfield,NUMpnts,limits)<<"\n";

		//+10%Y+10%Z
		YZ2FarField(tempfarfield,plusYplot,plusZplot);
		intensity(tempfarfield,NUMpnts);
		//cout<<directivity(tempfarfield,NUMpnts,limits)<<"\n";
		if(directivity(tempfarfield,NUMpnts,limits)>directivity(farfield,NUMpnts,limits))
		{
			copypaste(Yplot,plusYplot,NUMsample);
			copypaste(Zplot,plusZplot,NUMsample);
			copypaste(farfield,tempfarfield,NUMpnts);
		}

		//+10%Y-10%Z
		YZ2FarField(tempfarfield,plusYplot,minusZplot);
		intensity(tempfarfield,NUMpnts);
		//cout<<directivity(tempfarfield,NUMpnts,limits)<<"\n";
		if(directivity(tempfarfield,NUMpnts,limits)>directivity(farfield,NUMpnts,limits))
		{
			copypaste(Yplot,plusYplot,NUMsample);
			copypaste(Zplot,minusZplot,NUMsample);
			copypaste(farfield,tempfarfield,NUMpnts);
		}
		
		//-10%Y+10%Z
		YZ2FarField(tempfarfield,minusYplot,plusZplot);
		intensity(tempfarfield,NUMpnts);
		//cout<<directivity(tempfarfield,NUMpnts,limits)<<"\n";
		if(directivity(tempfarfield,NUMpnts,limits)>directivity(farfield,NUMpnts,limits))
		{	
			copypaste(Yplot,minusYplot,NUMsample);
			copypaste(Zplot,plusZplot,NUMsample);
			copypaste(farfield,tempfarfield,NUMpnts);
		}

		//-10%Y-10%Z
		YZ2FarField(tempfarfield,minusYplot,minusZplot);
		intensity(tempfarfield,NUMpnts);
		//cout<<directivity(tempfarfield,NUMpnts,limits)<<"\n\n";
		if(directivity(tempfarfield,NUMpnts,limits)>directivity(farfield,NUMpnts,limits))
		{
			copypaste(Yplot,minusYplot,NUMsample);
			copypaste(Zplot,minusZplot,NUMsample);
			copypaste(farfield,tempfarfield,NUMpnts);
		}
	}
}
