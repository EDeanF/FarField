#include "Header.h"
#include "arrays.h"
#include "optimize.h"
#include "functions.h"

// declare plot and printing functions
void MakePlot(cmplxpnt plot[], int plotsize, double start, double stop, complex<double> (*func)(double));
void MakeFarField(cmplxpnt plot[], int plotsize, double start, double stop, complex<double> (*func)(double, complex<double> (*)(double x)),complex<double>(*func2FT)(double));
void MakeFarField(cmplxpnt farfield[], int plotsize, double start, double stop, cmplxpnt nearfield[]);
void PrntPlot( cmplxpnt plot[], int plotsize, string title, string Y);
string MakeTitle();

int main()
{
	string title=MakeTitle();
	//Make Zplot and Yplot
	cmplxpnt Yplot[NUMsample]; cmplxpnt Zplot[NUMsample]; 
	MakeLimitYZ(Yplot,Zplot);

	//Make FarField
	cmplxpnt farfield[NUMpnts];
	YZ2FarField(farfield,Yplot,Zplot);
	intensity(farfield, NUMpnts);

	//Print YZ and farfield
	normalizeY(Yplot,NUMsample);
	normalizeZ(Zplot,NUMsample);
	radians2degrees(farfield,NUMpnts);
	
	PrntPlot(Yplot,NUMsample,"Y_limit_"+title,"Y");
	PrntPlot(Zplot,NUMsample,"Z_limit_"+title,"Z");
	PrntPlot(farfield,NUMpnts,"FF_limit_"+title,"I");

	int limits[2];
	FindLimits(limits,farfield,NUMpnts);

	//show limits to make sure we aren't missing large portions
	cout<<farfield[limits[0]].x<<"\t"<<farfield[limits[1]].x<<"\n";

	revertY(Yplot,NUMsample);
	revertZ(Zplot,NUMsample);
	degrees2radians(farfield,NUMpnts);
	
	for(int trial=0;trial<NUMtrials;++trial){
		optimize(farfield,Yplot,Zplot,limits);
	}

	//Print YZ and farfield
	normalizeY(Yplot,NUMsample);
	normalizeZ(Zplot,NUMsample);
	radians2degrees(farfield,NUMpnts);
	
	PrntPlot(Yplot,NUMsample,"Y_optimized_"+title,"Y");
	PrntPlot(Zplot,NUMsample,"Z_optimized_"+title,"Z");
	PrntPlot(farfield,NUMpnts,"FF_optimized_"+title,"I");

	revertY(Yplot,NUMsample);
	revertZ(Zplot,NUMsample);
	degrees2radians(farfield,NUMpnts);

	return 0;

}
