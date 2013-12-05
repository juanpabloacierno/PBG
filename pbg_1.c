/*
ParametricBackboneGenerator 
Authors: Juan Pablo Acierno and Mario R. Ermacora


Input file example:

AAAAAAA
AAAAAAB
AAAAAAC
AAAAAAD
AAAAAAE
AAAAAAF
AAAAABA
AAAAABB
AAAAABC
	.
	.
	.
FFFFFFF

*/


#define maxlibsize 300000
#define MAXLINEALONG 120
#define TOTAL_ATOM_NUMBER 36
#define COL_XYZ 3
#define SI_ATOM_NUMBER 10
#define PDBSIZE 6000

double crd1[36][3];//contiene las coordenadas del miembro de la biblioteca que se van a superponer al target
double crd2[36][3];//contiene las coordenadas del target
double crd3[36][3];//contiene las coordenadas del miembro de la biblioteca ya superpuestas al target
double tm[3][3];//contiene la última matriz de superposicion calculada
double translation1[3];
double translation2[3];
//---------------------------------------------------------------
//parameters initialized by reading the file RGBparameters.txt
long int libsize;
char library[50],pdbslist[50];
double rmsdcutof1;
int wpdb;
//----------------------------------------------------------------
struct lib6to7 {
 int numModel;
 char combination[8];
 double coord[TOTAL_ATOM_NUMBER][COL_XYZ];
} member[maxlibsize];

struct libout {
 int numModel;
 char combination[8];
 double targetC1libraryC2dis;
 double targetC9libraryC8dis;
 double rmsd_ends;
 double rmsd_N;
 double rmsd_C;
 double rmsd_global;
 double coord[TOTAL_ATOM_NUMBER][COL_XYZ];
} memberout[maxlibsize];


struct BackboneCoordinates {
 int resn;
 char atomcode[4];
 double coord[COL_XYZ];
} targetdata[PDBSIZE];

struct TargetInfo {
 char aPDB[MAXLINEALONG];
 int targettotalbbatoms;
 int numberofres;
 int numberof9Mers;
} targetinfo;


#include "include/somefunctiondeclarations.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include "include/bool.h"
#include "include/nrutil_FitCoord.h"
#include "include/parametersinput.h"
//#include "include/canonical_dihedrals.h"
#include "include/alignment.h"
#include "include/write_pdb.h"
#include "include/GetTargetBB.h"
#include "include/CheckMyRecord.h"
#include "include/xyzForFitOfEnds.h"
#include "include/superimpose.h"
#include "include/rmsdcalculation.h"
#include "include/move.h"
#include "include/xyzForGlobalFit.h"
#include "include/WriteFitLibrary.h"
#include "include/C1C2andC8C9dis.h"
#include "include/ReadLibraryFromFile.h"


int main()
{

int a1=0, h=0, i=0, j=0, k=0;

char scrap[MAXLINEALONG],pdbname[MAXLINEALONG];

//collects run parameters from file
parametersin();

//open a file for data output
FILE *mg;
mg = fopen("checkout.dat","w");
fprintf(mg,
		"pdbname\toutputindex\t9mer\tmember\trmsdends\trmsdN\trmsdC\trmsdglobal\tgap12\tgap89\n");

/*Read peptide library from File*/
readlibraryfromfile();


//find the list of pdb targets and handle one file name at the time for processing
FILE *fy;
	fy = fopen(pdbslist,"r");
	if (fy==NULL){
	
		printf("Can find pdbnames file!\n");
		printf("Enter any character to exit the program:\t");
		scanf("%s",scrap);
		exit(1);
    }

while (fscanf(fy,"%s", pdbname)!=EOF) {
//load and check target bb coordinates and info
	GetTargetBB(pdbname);//in GetTargetBB.h


for (i=1;i<=targetinfo.numberof9Mers; i++) {//primary loop i index is over all the 9mers in a protein chain
	a1=i+8;
	k++;
	
	
	for (j=0;j<libsize;j++) {//secondary loop j index enumerates all the library members
	
	//fit the library member ends to the target ends and fils some data in the library out array 
	memberout[j].numModel=member[j].numModel;
	strcpy(memberout[j].combination,member[j].combination);
	coordforfitofends(i, j);/*i=iesimo9mer; j=jesimolibmember*/
	
	superimpose(SI_ATOM_NUMBER);//superimpose ends
	
	//Calcula RMSDs de los extremos
	memberout[j].rmsd_ends=rmsdcalculation(0, SI_ATOM_NUMBER);//parameters are the limits of the atom index included in rmsd calculations
	
	memberout[j].rmsd_N=rmsdcalculation(0, 4);//parameters are the limits of the atom index included in rmsd calculations
	
	memberout[j].rmsd_C=rmsdcalculation(5, SI_ATOM_NUMBER);//parameters are the limits of the atom index included in rmsd calculations
	
	
	/*Aplica a las coordenadas del jesimo miembro de la biblioteca la translacion calculada para
	llevar el centro geometrico de los extremos a 0,0,0 (translation1).
	Aplica a las coordenadas del jesimo miembro la matriz de rotacion
	calculada por superimpose para los extremos (tm).
	Vuelve las coordenadas rotadas al centro geometrico original de los extremos (translation2).
	El resultado es la superposición global del jesimo miembro a iesimo target 9mer que minimiza
	el rmsd de los extremos*/
	move(j);
	
	

	/*select atoms for global fit (el ajuste para minimzar rmsd de todos
	los atomos del jesimo miembro al target*/
	coordforglobalfit(i,j);
	
	//global fit
	superimpose(TOTAL_ATOM_NUMBER-3);
	memberout[j].rmsd_global=rmsdcalculation(0,TOTAL_ATOM_NUMBER-4);
	
	/*calculates closure gap (las distancias entre CA1 y CA9 del blanco y CA2 y CA8
	del miembro de la bibioteca respectivamente)*/
	C1C2andC8C9dis(i,j);
	
	}//ends secondary loop (library members)
	
	//This section organize and select the results of library members
	
		h=0;
	for (j=0;j<libsize;j++) {
	if (memberout[j].rmsd_N<=rmsdcutof1&&memberout[j].rmsd_C<=rmsdcutof1&&memberout[j].targetC1libraryC2dis<=4.0&&
		memberout[j].targetC1libraryC2dis>=3.6&&memberout[j].targetC9libraryC8dis<=4.0&&
		memberout[j].targetC9libraryC8dis>=3.6) {
	/*if (memberout[j].rmsd_global<=rmsdcutof1) {*/
		if (wpdb==1) WriteFitLibrary(pdbname, i,j);
		fprintf(mg,
		"%s\t%i\t%i\t%i\t%8.3f\t%8.3f\t%8.3f\t%8.3f\t%8.3f\t%8.3f\n",
		pdbname,h,i,j,memberout[j].rmsd_ends, memberout[j].rmsd_N, memberout[j].rmsd_C, memberout[j].rmsd_global,
		memberout[j].targetC1libraryC2dis,memberout[j].targetC9libraryC8dis);
		h++;
	}
	}
	
}//ends primary loop (target protein 9mers)
}
//printf("Enter whatever\n");
//scanf("%s",scrap);
fclose(mg);
return 0;
}

