/*
************************************************************************
 Herein are defined the canonical tripeptides coordinates used to build
 the combinatorial library. Current order of atoms is (N-CA-C-O)*3
 To obtain a new library, simpli define a new tripeptide with the
 desired phi-psi pair and replace the coordinates below.
************************************************************************
*/
#ifndef FI
#define FI 12
#endif

#ifndef COL3
#define COL3 3
#endif

double std25[8][3]={
{-3.610,-0.654,0.017},
{-2.432,0.308,0.088},
{-2.621,1.515,0.217},
{-1.213,-0.231,0.004},
};

/* canonical pair A phi = -150  psi = 150 */
double A[FI][COL3]={
{-4.859,0.073,0.122},
{-3.610,-0.654,0.017},
{-2.432,0.308,0.088},
{-2.621,1.515,0.217},
{-1.213,-0.231,0.004},
{-0.016,0.584,0.061},
{1.100,-0.061,-0.749},
{1.157,-1.283,-0.869},
{1.989,0.764,-1.306},
{3.095,0.267,-2.100},
{3.942,1.426,-2.607},
{3.646,2.585,-2.327},
};
/* canonical pair B phi = -70  psi = 140 */
double B[FI][COL3]={
{-4.859,0.073,0.122},
{-3.610,-0.654,0.017},
{-2.432,0.308,0.088},
{-2.621,1.515,0.217},
{-1.213,-0.231,0.004},
{-0.016,0.584,0.061},
{0.113,1.420,-1.204},
{-0.167,0.938,-2.300},
{0.540,2.675,-1.052},
{0.703,3.566,-2.184},
{1.199,4.927,-1.716},
{1.403,5.141,-0.522},
};
/* canonical pair C phi = -120  psi = 120 */
double C[FI][COL3]={
{-4.859,0.073,0.122},
{-3.610,-0.654,0.017},
{-2.432,0.308,0.088},
{-2.621,1.515,0.217},
{-1.213,-0.231,0.004},
{-0.016,0.584,0.061},
{0.794,0.422,-1.218},
{1.222,-0.680,-1.549},
{1.005,1.527,-1.937},
{1.760,1.497,-3.174},
{1.839,2.893,-3.777},
{1.303,3.846,-3.217},
};
/* canonical pair E phi = 60  psi = 30 */
double D[FI][COL3]={
{-4.859,0.073,0.122},
{-3.610,-0.654,0.017},
{-2.432,0.308,0.088},
{-2.621,1.515,0.217},
{-1.213,-0.231,0.004},
{-0.016,0.584,0.061},
{0.045,1.339,1.382},
{0.601,2.433,1.448},
{-0.527,0.751,2.435},
{-0.532,1.373,3.744},
{-1.250,0.483,4.748},
{-1.726,-0.595,4.398},
};
/* canonical pair F phi = -100  psi = 20 */
double E[FI][COL3]={
{-4.859,0.073,0.122},
{-3.610,-0.654,0.017},
{-2.432,0.308,0.088},
{-2.621,1.515,0.217},
{-1.213,-0.231,0.004},
{-0.016,0.584,0.061},
{0.527,0.822,-1.341},
{1.696,1.161,-1.506},
{-0.327,0.646,-2.351},
{0.075,0.843,-3.730},
{-1.096,0.572,-4.664},
{-2.189,0.237,-4.214},
};
/* canonical pair H phi = -65  psi = -40 */
double F[FI][COL3]={
{-4.859,0.073,0.122},
{-3.610,-0.654,0.017},
{-2.432,0.308,0.088},
{-2.621,1.515,0.217},
{-1.213,-0.231,0.004},
{-0.016,0.584,0.061},
{0.050,1.509,-1.146},
{0.436,2.669,-1.021},
{-0.330,0.994,-2.317},
{-0.311,1.778,-3.536},
{-0.792,0.939,-4.710},
{-1.128,-0.232,-4.545},
};

/********* THIS SET OF 9 DIHEDRALS IS NOT USED **************
// canonical pair A phi = -150  psi = 150 /
double A[FI][COL3]={
{-4.859,0.073,0.122},
{-3.610,-0.654,0.017},
{-2.432,0.308,0.088},
{-2.621,1.515,0.217},
{-1.213,-0.231,0.004},
{-0.016,0.584,0.061},
{1.100,-0.061,-0.749},
{1.157,-1.283,-0.869},
{1.989,0.764,-1.306},
{3.095,0.267,-2.100},
{3.942,1.426,-2.607},
{3.646,2.585,-2.327},
};
//canonical pair B phi = -70  psi = 140
double B[FI][COL3]={
{-4.859,0.073,0.122},
{-3.610,-0.654,0.017},
{-2.432,0.308,0.088},
{-2.621,1.515,0.217},
{-1.213,-0.231,0.004},
{-0.016,0.584,0.061},
{0.113,1.420,-1.204},
{-0.167,0.938,-2.300},
{0.540,2.675,-1.052},
{0.703,3.566,-2.184},
{1.199,4.927,-1.716},
{1.403,5.141,-0.522},
};
//canonical pair C phi = -120  psi = 120
double C[FI][COL3]={
{-4.859,0.073,0.122},
{-3.610,-0.654,0.017},
{-2.432,0.308,0.088},
{-2.621,1.515,0.217},
{-1.213,-0.231,0.004},
{-0.016,0.584,0.061},
{0.794,0.422,-1.218},
{1.222,-0.680,-1.549},
{1.005,1.527,-1.937},
{1.760,1.497,-3.174},
{1.839,2.893,-3.777},
{1.303,3.846,-3.217},
};
//canonical pair D phi = -80  psi = 120
double D[FI][COL3]={
{-4.859,0.073,0.122},
{-3.610,-0.654,0.017},
{-2.432,0.308,0.088},
{-2.621,1.515,0.217},
{-1.213,-0.231,0.004},
{-0.016,0.584,0.061},
{0.247,1.229,-1.292},
{0.429,0.533,-2.289},
{0.269,2.563,-1.325},
{0.510,3.290,-2.556},
{0.478,4.789,-2.298},
{0.268,5.224,-1.167},
};
//canonical pair E phi = 60  psi = 30
double E[FI][COL3]={
{-4.859,0.073,0.122},
{-3.610,-0.654,0.017},
{-2.432,0.308,0.088},
{-2.621,1.515,0.217},
{-1.213,-0.231,0.004},
{-0.016,0.584,0.061},
{0.045,1.339,1.382},
{0.601,2.433,1.448},
{-0.527,0.751,2.435},
{-0.532,1.373,3.744},
{-1.250,0.483,4.748},
{-1.726,-0.595,4.398},
};
//canonical pair F phi = -100  psi = 20
double F[FI][COL3]={
{-4.859,0.073,0.122},
{-3.610,-0.654,0.017},
{-2.432,0.308,0.088},
{-2.621,1.515,0.217},
{-1.213,-0.231,0.004},
{-0.016,0.584,0.061},
{0.527,0.822,-1.341},
{1.696,1.161,-1.506},
{-0.327,0.646,-2.351},
{0.075,0.843,-3.730},
{-1.096,0.572,-4.664},
{-2.189,0.237,-4.214},
};
//canonical pair G phi = -85  psi = -20
double G[FI][COL3]={
{-4.859,0.073,0.122},
{-3.610,-0.654,0.017},
{-2.432,0.308,0.088},
{-2.621,1.515,0.217},
{-1.213,-0.231,0.004},
{-0.016,0.584,0.061},
{0.316,1.129,-1.320},
{1.045,2.111,-1.443},
{-0.221,0.491,-2.362},
{0.024,0.917,-3.726},
{-0.710,0.010,-4.702},
{-1.391,-0.928,-4.291},
};
//canonical pair H phi = -65  psi = -40
double H[FI][COL3]={
{-4.859,0.073,0.122},
{-3.610,-0.654,0.017},
{-2.432,0.308,0.088},
{-2.621,1.515,0.217},
{-1.213,-0.231,0.004},
{-0.016,0.584,0.061},
{0.050,1.509,-1.146},
{0.436,2.669,-1.021},
{-0.330,0.994,-2.317},
{-0.311,1.778,-3.536},
{-0.792,0.939,-4.710},
{-1.128,-0.232,-4.545},
};
//canonical pair I phi = 60  psi = -150
double I[FI][COL3]={
{-4.859,0.073,0.122},
{-3.610,-0.654,0.017},
{-2.432,0.308,0.088},
{-2.621,1.515,0.217},
{-1.213,-0.231,0.004},
{-0.016,0.584,0.061},
{0.045,1.339,1.382},
{-0.460,0.862,2.396},
{0.666,2.520,1.368},
{0.788,3.330,2.564},
{1.539,4.617,2.255},
{1.945,4.845,1.117},
};
*/
