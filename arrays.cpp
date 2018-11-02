#include "Header.h"

void normalizeY(cmplxpnt Yplot[], int plotsize)
{
	for(int i=0;i<plotsize;++i)
	{
		Yplot[i].y=Yplot[i].y*eta;
	}
}

void normalizeZ(cmplxpnt Zplot[], int plotsize)
{
	for(int i=0;i<plotsize;++i)
	{
		Zplot[i].y=Zplot[i].y/eta;
	}
}

void revertY(cmplxpnt Yplot[], int plotsize)
{
	for(int i=0;i<plotsize;++i)
	{
		Yplot[i].y=Yplot[i].y/eta;
	}
}

void revertZ(cmplxpnt Zplot[], int plotsize)
{
	for(int i=0;i<plotsize;++i)
	{
		Zplot[i].y=Zplot[i].y*eta;
	}
}

void limitYZ(cmplxpnt plot[], int plotsize)
{
	for(int i=0; i<plotsize; ++i)
	{
		//if imaginary part is greater than 3
		if(abs(plot[i].y.imag())>3)
		{
			if(plot[i].y.imag()>0){plot[i].y.imag(3);}
			else{plot[i].y.imag(-3);}
		}
		//if imaginary part is less than 1/3
		if(abs(plot[i].y.imag())<1e0/3)
		{
			if(plot[i].y.imag()>=0){plot[i].y.imag(1e0/3);}
			else{plot[i].y.imag(-1e0/3);}
		}
		//if real part is greater than 3
		int temp=i;
		while(abs(plot[i].y.real())>3)
		{
			if(temp==0)
			{
				temp=plotsize;
				plot[i].y.real(plot[temp].y.real());
			}
			else
			{
				temp-=1; 
				plot[i].y.real(plot[temp].y.real());
			}
		}
	}
}

void MakeRplot(cmplxpnt Rplot[], int plotsize, cmplxpnt Zplot[], cmplxpnt Yplot[])
{
	for(int i=0;i<plotsize;++i)
	{
		Rplot[i].y=-0.5/cos(theta)*(Yplot[i].y*eta*pow(cos(theta),2)-Zplot[i].y/eta)/(1e0+0.25*Zplot[i].y*Yplot[i].y+0.5*(Yplot[i].y*eta*pow(cos(theta),2)+Zplot[i].y/eta));
	}
}

void MakeTplot(cmplxpnt Tplot[], int plotsize, cmplxpnt Zplot[], cmplxpnt Yplot[])
{
	for(int i=0;i<plotsize;++i)
	{
		Tplot[i].y=(1.e0-0.25*Zplot[i].y*Yplot[i].y)/(1e0+0.25*Zplot[i].y*Yplot[i].y+0.5*(Yplot[i].y*eta*pow(cos(theta),2)+Zplot[i].y/eta));
	}
}

void MakeE2plot(cmplxpnt E2plot[], int plotsize, cmplxpnt EIplot[], cmplxpnt Tplot[])
{
	for(int i=0;i<plotsize;++i)
	{
		E2plot[i].y=EIplot[i].y*Tplot[i].y;
	}
}

void MakeH2plot(cmplxpnt H2plot[], int plotsize, cmplxpnt HIplot[], cmplxpnt Tplot[])
{
	for(int i=0;i<plotsize;++i)
	{
		H2plot[i].y=HIplot[i].y*Tplot[i].y;
	}
}

void u2radians(cmplxpnt plot[], int plotsize)
{
	for(int i=0;i<plotsize;++i)
	{
		plot[i].x=asin(plot[i].x*lambda);
	}
}

void radians2degrees(cmplxpnt plot[], int plotsize)
{
	for(int i=0;i<plotsize;++i)
	{
		plot[i].x=plot[i].x*180/pi;
	}
}

void intensity(cmplxpnt plot[], int plotsize)
{
	for(int i=0;i<plotsize;++i)
	{
		plot[i].y=c*eps/2*pow(abs(plot[i].y),2);
	}
}

void degrees2radians(cmplxpnt plot[],int plotsize)
{
	for(int i=0;i<plotsize;++i)
	{
		plot[i].x=plot[i].x*pi/180;
	}
}
