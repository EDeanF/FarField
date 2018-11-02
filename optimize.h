void MakeIdealYZ(cmplxpnt Yplot[],cmplxpnt Zplot[]);
void MakeLimitYZ(cmplxpnt Yplot[],cmplxpnt Zplot[]);
void YZ2FarField(cmplxpnt farfield[],cmplxpnt Yplot[],cmplxpnt Zplot[]);
void FindLimits(int limits[2], cmplxpnt plot[], int plotsize);
double integrate(cmplxpnt plot[], int lwrlim, int upprlim);
double directivity(cmplxpnt plot[], int plotsize);
double directivity(cmplxpnt plot[], int plotsize, int limits[2]);
void copypaste(cmplxpnt paste[],cmplxpnt copy[], int plotsize);
void optimize(cmplxpnt farfield[NUMpnts],cmplxpnt Yplot[NUMsample],cmplxpnt Zplot[NUMsample],int limits[2]);
