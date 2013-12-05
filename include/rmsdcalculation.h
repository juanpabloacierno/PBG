double rmsdcalculation(int from, int m){
int i=0, j=0;
double suma=0.0, tobesquared, rmsd;

	for (i=from;i<=(m-1);i++)
	{
	for (j=0;j<=(COL_XYZ-1);j++)
	{
		tobesquared=crd3[i][j]-crd2[i][j];
		suma+=(tobesquared*tobesquared);
	}
	
	}
	rmsd= sqrtf(suma)/sqrtf(m);
	
	return rmsd;
}

