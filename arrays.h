void normalizeY(cmplxpnt Yplot[], int plotsize);
void normalizeZ(cmplxpnt Zplot[], int plotsize);
void revertY(cmplxpnt Yplot[], int plotsize);
void revertZ(cmplxpnt Zplot[], int plotsize);
void limitYZ(cmplxpnt plot[], int plotsize);
void MakeRplot(cmplxpnt Rplot[], int plotsize, cmplxpnt Zplot[], cmplxpnt Yplot[]);
void MakeTplot(cmplxpnt Tplot[], int plotsize, cmplxpnt Zplot[], cmplxpnt Yplot[]);
void MakeE2plot(cmplxpnt E2plot[], int plotsize, cmplxpnt EIplot[], cmplxpnt Tplot[]);
void MakeH2plot(cmplxpnt H2plot[], int plotsize, cmplxpnt HIplot[], cmplxpnt Tplot[]);

double integrate(cmplxpnt plot[], int upprlim, int lwrlim);
void FindLimits(int limits[2], cmplxpnt plot[], int plotsize);
double directivity(cmplxpnt plot[], int plotsize);

void u2radians(cmplxpnt plot[], int plotsize);
void radians2degrees(cmplxpnt plot[], int plotsize);
void intensity(cmplxpnt plot[], int plotsize);
void degrees2radians(cmplxpnt plot[],int plotsize);