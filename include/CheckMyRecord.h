/*int CheckMyRecord(char *, int *, int *, double *, double *,
double *, char *, int *, char *);

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
*/
/* Verifies PDB file register
COLUMNS DATA TYPE FIELD DEFINITION
-------------------------------------------------------------
1 - 6 			Record name "ATOM "
7 - 11		 	Integer serial Atom serial number.
13 - 16 		Atom name Atom name.
17 				Character altLoc Alternate location indicator.
18 - 20 		Residue name resName Residue name.
22 				Character chainID Chain identifier.
23 - 26 		Integer resSeq Residue sequence number.
27 				AChar iCode Code for insertion of residues.
31 - 38 		Real(8.3) x Orthogonal coordinates for X in Angstroms
39 - 46			Real(8.3) y Orthogonal coordinates for Y in Angstroms
47 - 54         Real(8.3) z Orthogonal coordinates for Z in Angstroms
55 - 60         Real(6.2) occupancy Occupancy.
61 - 66         Real(6.2) tempFactor Temperature factor.
77 - 78         LString(2) element Element symbol, right-justified.
79 - 80         LString(2) charge Charge on the atom.
Example
         1         2         3         4         5         6         7         8
12345678901234567890123456789012345678901234567890123456789012345678901234567890
ATOM    145  N   VAL A  25      32.433  16.336  57.540  1.00 11.92           N1+
ATOM    146  CA  VAL A  25      31.132  16.439  58.160  1.00 11.85           C
ATOM    147  C   VAL A  25      30.447  15.105  58.363  1.00 12.34           C
ATOM    148  O   VAL A  25      29.520  15.059  59.174  1.00 15.65           O
ATOM    149  CB AVAL A  25      30.385  17.437  57.230  0.28 13.88           C
ATOM    150  CB BVAL A  25      30.166  17.399  57.373  0.72 15.41           C
ATOM    151  CG1AVAL A  25      28.870  17.401  57.336  0.28 12.64           C
ATOM    152  CG1BVAL A  25      30.805  18.788  57.449  0.72 15.11           C
ATOM    153  CG2AVAL A  25      30.835  18.826  57.661  0.28 13.58           C
ATOM    154  CG2BVAL A  25      29.909  16.996  55.922  0.72 13.25           C
...
...
ATOM    588 1HG  GLU    18     -13.363  -4.163  -2.372  1.00  0.00           H
ATOM    589 2HG  GLU    18     -12.634  -3.023  -3.475  1.00  0.00           H
*/
int CheckMyRecord(char *record, int *isatom, int *isres, double *xxx,
	double *yyy, double *zzz, char *atmnameatm, int *resseqnumres, char *altlocalt)
{
int b=0,i,j;
char recordname[7];
int serial;
char serialstring[6];
char atmname[5];
char altloc;
char resname[4];
char chainID; 
int resseqnum;
char resseqnumstring[5];
char iCode;
double xx;
char xxstring[9];
double yy;
char yystring[9];
double zz;
char zzstring[9];
float occupancy;
char occupancystring[7];
float tempfact;
char tempfactstring[7];
char element[3];
char charge[3];
char **endp;


/*-----------------------------------------------------*/
	for (i=0,j=0;i<=5;i++,j++) recordname[j]=record[i];
	recordname[6]='\0';	
	
if (strncmp(recordname,"ATOM  ",6)==0 || strncmp(recordname,"HETATM",6)==0) {
	*isatom=1;}
	 else {
		 *isatom=0;
	 }
	
/*-----------------------------------------------------*/
	for (i=6,j=0;i<=10;i++,j++) serialstring[j]=record[i];
	serialstring[5]='\0';
	serial=atoi(serialstring);
	

/*-----------------------------------------------------*/
	for (i=12,j=0;i<=15;i++,j++) atmname[j]=record[i];
	atmname[4]='\0';
	strcpy(atmnameatm,atmname);
	

/*-----------------------------------------------------*/
	altloc=record[16];
	*altlocalt=record[16];

/*-----------------------------------------------------*/
	for (i=17,j=0;i<=19;i++,j++) resname[j]=record[i];
	resname[3]='\0';	
	

/*-----------------------------------------------------*/
	chainID=record[21];

/*-----------------------------------------------------*/
	for (i=22,j=0;i<=25;i++,j++) resseqnumstring[j]=record[i];
	resseqnumstring[4]='\0';	
	resseqnum=atoi(resseqnumstring);
	*resseqnumres=resseqnum;

/*-----------------------------------------------------*/
	iCode=record[26];

/*-----------------------------------------------------*/
	for (i=30,j=0;i<=37;i++,j++) xxstring[j]=record[i];
	xxstring[8]='\0';	
	xx=strtod(xxstring,endp);
	*xxx=xx;

/*-----------------------------------------------------*/
	for (i=38,j=0;i<=45;i++,j++) yystring[j]=record[i];
	yystring[8]='\0';
	yy=strtod(yystring,endp);
	*yyy=yy;

/*-----------------------------------------------------*/
	for (i=46,j=0;i<=53;i++,j++) zzstring[j]=record[i];
	zzstring[8]='\0';	
	zz=strtod(zzstring,endp);
	*zzz=zz;

/*-----------------------------------------------------*/
	for (i=54,j=0;i<=59;i++,j++) occupancystring[j]=record[i];
	occupancystring[6]='\0';
	occupancy=atof(occupancystring);
	

/*-----------------------------------------------------*/
	for (i=60,j=0;i<=65;i++,j++) tempfactstring[j]=record[i];
	tempfactstring[6]='\0';
	tempfact=atof(tempfactstring);
	

/*-----------------------------------------------------*/
	for (i=76,j=0;i<=77;i++,j++) element[j]=record[i];
	element[2]='\0';	
	

/*-----------------------------------------------------*/
	for (i=78,j=0;i<=79;i++,j++) charge[j]=record[i];
	charge[2]='\0';	
	

/*-----------------------------------------------------*/


if ((isdigit(atmname[0]) || atmname[0]==' ' || atmname[0]=='S'|| atmname[0]=='H') && (atmname[1]=='C'||
	atmname[1]=='N'|| atmname[1]=='O' || atmname[1]=='H'|| atmname[1]=='S'|| atmname[1]=='E'||
	atmname[1]=='G' || atmname[1]=='D' || isdigit(atmname[1])) &&
	(atmname[2]=='A' || atmname[2]=='B' || atmname[2]=='G' || atmname[2]=='D' ||
	 atmname[2]=='E' || atmname[2]=='Z' || atmname[2]=='H' || atmname[2]==' ' ||
	isdigit(atmname[2]) || atmname[2]=='T' || atmname[2]=='X') 
	&& (isdigit(atmname[3]) || atmname[3]==' ' || atmname[3]=='T')) ++b;
	
	
if (isupper(altloc)||altloc==' ') ++b;
	

if (strncmp(resname,ALA,3)==0||strncmp(resname,CYS,3)==0||strncmp(resname,ASP,3)==0||strncmp(resname,GLU,3)==0||strncmp(resname,PHE,3)==0||
	strncmp(resname,GLY,3)==0||strncmp(resname,HIS,3)==0||strncmp(resname,ILE,3)==0||strncmp(resname,LYS,3)==0||strncmp(resname,LEU,3)==0||
	strncmp(resname,MET,3)==0||strncmp(resname,ASN,3)==0||strncmp(resname,PRO,3)==0||strncmp(resname,GLN,3)==0||strncmp(resname,ARG,3)==0||
	strncmp(resname,SER,3)==0||strncmp(resname,THR,3)==0||strncmp(resname,VAL,3)==0||strncmp(resname,TRP,3)==0||strncmp(resname,TYR,3)==0||
	strncmp(resname,MSE,3)==0||strncmp(resname,AXN,3)==0||strncmp(resname,GXN,3)==0||strncmp(resname,PCA,3)==0)
	{
		*isres=1;
		++b;}
	else *isres=0;

if (isupper(chainID)||chainID==' ') ++b;
	

if (resseqnum<=MAXRESNUM && resseqnum>=MINRESNUM) ++b;
	

if (xxstring[4]=='.' && yystring[4]=='.' && zzstring[4]=='.' && occupancystring[3]=='.' &&
	tempfactstring[3]=='.')
	++b;

if (xx<=MAXCOORDINATE && xx>=MINCOORDINATE && yy<=MAXCOORDINATE &&
	yy>=MINCOORDINATE && zz<=MAXCOORDINATE && zz>=MINCOORDINATE)
	++b;

if (occupancy<=MAXOCCUPANCY && occupancy>=MINOCCUPANCY) ++b;
	
if (tempfact<=MAXTEMP && tempfact>=MINTEMP) ++b;

if ((element[0]==' ' || element[0]=='S') &&
	((element[1]=='C')||(element[1]=='N')||
	(element[1]=='O')||(element[1]=='H')||(element[1]=='S')||(element[1]=='E')))
	++b;

return b;
}
