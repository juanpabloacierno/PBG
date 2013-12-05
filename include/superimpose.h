
int superimpose(long int m){

	int i,j,k;
	long int n=3;
	extern double crd1[][];
	extern double crd2[][];
	extern double crd3[][];
	extern double tm[][];
	extern double translation1[];
	extern double translation2[];
	
//Calcula las coordenadas del centro geométrico de ambas estruturas

//vector que contendrá las coordenadas del centro geométrico de la primer molécula
	double *mv1;
	mv1=vector(0,n-1);
	
	for (j=0;j<=(n-1);j++) mv1[j]=0;
	
	
	for (j=0;j<=(n-1);j++)
		{
			for (i=0;i<=(m-1); i++) mv1[j]+=crd1[i][j];
			mv1[j]=mv1[j]/m;
			
		}
	
//vector que contendrá las coordenadas del centro geométrico de la segunda molécula
	double *mv2;
	mv2=vector(0,n-1);
	
	for (j=0;j<=(n-1);j++) mv2[j]=0;
	
	for (j=0;j<=(n-1);j++)
		{
			for (i=0;i<=(m-1); i++) mv2[j]+=crd2[i][j];
			mv2[j]=mv2[j]/m;
		}


//Cálculos para Single Value Decomposition
	//Matrix with target coordinates x,y,z
	
	double **bplus;
	double **cplus;
	
	//Crea dos matrices de m-1 filas por n-1 columnas
	bplus=matrix(0,m-1,0,n-1);
	cplus=matrix(0,m-1,0,n-1);
	
	//llena matrices con las estructuras centradas en el mismo eje de coordenadas,
	//el centro geométrico de ambas está en 0,0,0*/
	
	
	
	
	for (i=0;i<=(m-1);i++)
	{
	for (j=0;j<=(n-1);j++)bplus[i][j]=crd1[i][j]-mv1[j];
	}
	
	for (i=0;i<=(m-1);i++)
	{
	for (j=0;j<=(n-1);j++)cplus[i][j]=crd2[i][j]-mv2[j];
	}
	
	
	
	
//calcula el producto XtY. Multiplicación de una matriz transpuesta por otra matriz*/
	
	double **xty;
	xty=matrix(0,n-1,0,n-1);
	
	for (i=0;i<=(n-1);i++)
	{
	for (j=0;j<=(n-1);j++)
		xty[i][j]=0;
	}
	
	for (i=0;i<=(n-1);i++)
	{
	for (j=0;j<=(n-1);j++)
		for (k=0;k<=(m-1);k++) xty[i][j]+=(bplus[k][i]*cplus[k][j]);
	}
	
	
	
//will be the diagonal vector from SVD
	
	double *cc;
	cc=vector(0,n-1);
	for (i=0;i<=(n-1);i++) cc[i]=0;
	
	
	
//will be the V matriz from SVD*/
	
	
	double **ddd;
	ddd=matrix(0,n-1,0,n-1);
	
	for (i=0;i<=(n-1);i++)
	{
	for (j=0;j<=(n-1);j++)
		ddd[i][j]=0;
	}
	
//will be the 1/w diagonal matriz from SVD ee[0...n-1][0...n-1]
	
	double **eee;
	eee=matrix(0,n-1,0,n-1);
	
	for (i=0;i<=(n-1);i++)
	{
	for (j=0;j<=(n-1);j++)
		eee[i][j]=0;
	}
	
	
SingValDecomp(xty, n, n, cc, ddd);
	
//crea la matriz diagonal inversa con los elementos de W*/
	
	for (i=0;i<=n-1;i++) eee[i][i]=1/cc[i];
	
	

//Multiplicación de una matriz por la transpuesta de otra: V por la transpuesta de U*/
	double **prd5;//prd5 contiene la matriz de superposicion calculada
	prd5=matrix(0,n-1,0,n-1);
	for (i=0;i<=(n-1);i++)
	{
	for (j=0;j<=(n-1);j++)
		prd5[i][j]=0;
	}
	for (i=0;i<=(n-1);i++)
	{
	for (j=0;j<=(n-1);j++)
		for (k=0;k<=(n-1);k++)prd5[i][j]+=(ddd[i][k]*xty[j][k]);
	}
//calcula el determinante
	
	float det=0;
	float detneg=0;
	
	det = prd5[0][0] * prd5[1][1] * prd5[2][2] +
	prd5[0][1] * prd5[1][2] * prd5[2][0] +
	prd5[0][2] * prd5[1][0] * prd5[2][1] -
	prd5[2][0] * prd5[1][1] * prd5[0][2] -
	prd5[2][1] * prd5[1][2] * prd5[0][0] -
	prd5[2][2] * prd5[1][0] * prd5[0][1];
	
	
//lo que sigue está tomado de MapCoord.c del MolMol, salvo el registro de eventos det negativos
	if (det < 0.0f) {
    if (cc[0] < cc[1] && cc[0] < cc[2])
      for (i = 0; i < 3; i++)
	for (k = 0; k < 3; k++)
	  prd5[i][k] =
	      - xty[k][0] * ddd[i][0] + xty[k][1] * ddd[i][1] + xty[k][2] * ddd[i][2];
    else if (cc[1] < cc[2])
      for (i = 0; i < 3; i++)
	for (k = 0; k < 3; k++)
	  prd5[i][k] =
	      xty[k][0] * ddd[i][0] - xty[k][1] * ddd[i][1] + xty[k][2] * ddd[i][2];
    else
      for (i = 0; i < 3; i++)
	for (k = 0; k < 3; k++)
	  prd5[i][k] =
	      xty[k][0] * ddd[i][0] + xty[k][1] * ddd[i][1] - xty[k][2] * ddd[i][2];
	detneg+=detneg;
  }
//aquí termina lo tomado de MapCoord.c del MolMol

	
	
//Multiplicación de una matriz por la transpuesta de otra, Obtiene las coordenadas de la superposición de crd1 sobre crd2 */
	
	
	
	for (i=0;i<=(m-1);i++)
	{
	for (j=0;j<=(n-1);j++)
		crd3[i][j]=0;
	}
	
	for (i=0;i<=(m-1);i++)
	{
	for (j=0;j<=(n-1);j++)
		for (k=0;k<=(n-1);k++)crd3[i][j]+=(bplus[i][k]*prd5[j][k]);
	}
	
	
	
//vuelve el centro geométrico de las coordenadas para la estructura rotada al lugar inicial de la estructura dos
	
	for (i=0;i<=(m-1);i++)
	{
	for (j=0;j<=(n-1);j++) crd3[i][j]=crd3[i][j]+mv2[j];
	}
	
//guarda la matriz de transformacion calculada prd5 como una matriz externa tm
	for (i=0;i<=(n-1);i++)
	{
	for (j=0;j<=(n-1);j++) tm[i][j]=prd5[i][j];
	}
	
	for (j=0;j<=(n-1);j++) translation1[j]=mv1[j];
	for (j=0;j<=(n-1);j++) translation2[j]=mv2[j];
	
	/*------------------------------------------
	cierra archivos y libera memoria*/ 
	
	free_vector(cc,0,n-1);
	free_matrix(bplus,0,m-1,0,n-1);
	free_matrix(cplus,0,m-1,0,n-1);
	free_matrix(ddd,0,n-1,0,n-1);
	free_matrix(xty,0,n-1,0,n-1);
	free_matrix(eee,0,n-1,0,n-1);
	free_matrix(prd5,0,n-1,0,n-1);
	
	return 0;
}


