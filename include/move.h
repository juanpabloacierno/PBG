int move(int jesimolibmember){

extern double tm[][];
extern double translation1[];
extern double translation2[];
double a[TOTAL_ATOM_NUMBER][COL_XYZ];
int i=0,j=0, k=0;


//memberout[jesimolibmember].coord[i][j]

	//library-------------------------------------------------
	for (i=0;i<=TOTAL_ATOM_NUMBER-1;i++){
		for (j=0;j<=COL_XYZ-1;j++)
		{
			a[i][j]=member[jesimolibmember].coord[i][j]-translation1[j];		
		}
	}
	//----------------------------------------------------------

	for (i=0;i<=(TOTAL_ATOM_NUMBER-1);i++)
	{
	for (j=0;j<=(COL_XYZ-1);j++)
		for (k=0;k<=(COL_XYZ-1);k++)memberout[jesimolibmember].coord[i][j]+=(a[i][k]*tm[j][k]);
	}
	
	for (i=0;i<=TOTAL_ATOM_NUMBER-1;i++){
		for (j=0;j<=COL_XYZ-1;j++)
		{
			memberout[jesimolibmember].coord[i][j]=memberout[jesimolibmember].coord[i][j]+translation2[j];		
		}
	}
return 0;
}
