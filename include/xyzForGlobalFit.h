//void coordforglobalfit(int, int)

//crea matrices para superponer las coordenadas de los nonapeptidos(residuos 1 a 9)

//ejemplo para iesimo9mer=1

//Residuo	atomo	x	y	z	tipo
//	1		1	x0	y0	z0	CA
//	1		2	x1	y1	z1	C
//	.		.	.	.	.
//	.		.	.	.	.
//	.		.	.	.	.
//	9		33	x9	y9	z9	CA

//ejemplo para jesimo library member=1

//Residuo	atomo	x	y	z	tipo
//	1		1	x0	y0	z0	CA
//	1		2	x1	y1	z1	C
//	.		.	.	.	.
//	.		.	.	.	.
//	.		.	.	.	.
//	9		33	x9	y9	z9	CA

void coordforglobalfit(int iesimo9mer, int jesimolibmember) {	
	
int n=COL_XYZ,k=0,l=0,j=0;
extern double crd1[][];
extern double crd2[][];



	//library-------------------------------------------------
	for (k=0,l=1;k<=32;k++,l++){
		for (j=0;j<=n-1;j++)
		{
			crd1[k][j]=member[jesimolibmember].coord[l][j];		
		}
	}
	
	//----------------------------------------------------------



	//target-------------------------------------------------
	for (k=0,l=-3;k<=32;k++,l++){
		for (j=0;j<=n-1;j++)
		{
			crd2[k][j]=targetdata[4*iesimo9mer+l].coord[j];		
		}
	}
	
	
}
