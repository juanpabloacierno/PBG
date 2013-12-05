#define MAXRESNUM 9999
#define MINRESNUM -999
#define MAXCOORDINATE 999.999
#define MINCOORDINATE -999.999
#define MAXOCCUPANCY 999.99
#define MINOCCUPANCY 0.00
#define MAXTEMP 999.99
#define MINTEMP 0.00
#define CHECKRECORDOK 9
#define LINES_READ 35000
#define ALA "ALA"
#define CYS "CYS"
#define ASP "ASP"
#define GLU "GLU"
#define PHE "PHE"
#define GLY "GLY"
#define HIS "HIS"
#define ILE "ILE"
#define LYS "LYS"
#define LEU "LEU"
#define MET "MET"
#define ASN "ASN"
#define PRO "PRO"
#define GLN "GLN"
#define ARG "ARG"
#define SER "SER"
#define THR "THR"
#define VAL "VAL"
#define TRP "TRP"
#define TYR "TYR"
#define MSE "MSE"
#define AXN "AXN"
#define GXN "GXN"
#define PCA "PCA"



int GetTargetBB(char *pdbname)
{
char line[MAXLINEALONG],atomname[4],alternatelocator;	
int nrecords=0,j=0,k,l=0,residue_number=0,isatomrecord,isresrecord,checkreturn;
double xval,yval,zval;



FILE *fz;
	fz = fopen("PBG.log","a");

	
//=====================================================
					//GET PDB DATA
//=====================================================


	
	FILE *fp;
	fp = fopen(pdbname,"r");	
	
	printf("%s\n",pdbname);	
strcpy(targetinfo.aPDB,pdbname);	
k=0;	
while((fgets(line,MAXLINEALONG,fp)!=NULL) && nrecords<LINES_READ)
	{
		checkreturn=CheckMyRecord(line, &isatomrecord, &isresrecord,&xval,&yval,&zval,atomname,
		&residue_number, &alternatelocator);
		if (isatomrecord==1 && checkreturn==10 && isresrecord==1 &&
			(strncmp(atomname," N  ",3)==0 || strncmp(atomname," CA ",3)==0 ||
			strncmp(atomname," C  ",3)==0 || strncmp(atomname," O  ",3)==0) &&
			(alternatelocator==' '|| alternatelocator=='A')) {
			
			
			
			targetdata[k].coord[0]=xval;
			targetdata[k].coord[1]=yval;
			targetdata[k].coord[2]=zval;
			
			strcpy(targetdata[k].atomcode,atomname);
			
			targetdata[k].resn=residue_number;
			//fprintf(fz,"ATOM  %i\t%8.3f\t%8.3f\t%8.3f\n",k,targetdata[k].coord[0],targetdata[k].coord[1],targetdata[k].coord[2]);
			
			k++;
			
		}
	}

while (j<k) {
	
if (strncmp(targetdata[j].atomcode, " N  ",3)==0 &&
strncmp(targetdata[j+1].atomcode, " CA ",3)==0 &&
strncmp(targetdata[j+2].atomcode, " C  ",3)==0 &&
strncmp(targetdata[j+3].atomcode, " O  ",3)==0) l++;
else fprintf(fz,"%s: error in target backbone",pdbname);
j+=4;
}
if (l!=k/4 || l<9)fprintf(fz,"%s: error in target backbone",pdbname);

targetinfo.targettotalbbatoms=k;
targetinfo.numberofres=l;
targetinfo.numberof9Mers=l-9+1;

fclose(fp);

fprintf(fz,"%s\t: target backbone OK, residues=%i, 9Mers=%i, total atoms=%i\n",
targetinfo.aPDB,targetinfo.numberofres,targetinfo.numberof9Mers,
targetinfo.targettotalbbatoms);

fclose (fz);
return 0;
}

