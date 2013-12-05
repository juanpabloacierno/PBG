
int parametersin()
{

extern char library[];
extern char pdbslist[];
extern long int libsize;
extern double rmsdcutof1;
extern int wpdb;

FILE *fp;
fp = fopen("RGBparameters.txt","r");

FILE *fz;
	fz = fopen("PBG.log","a");

char scrap[50],line[MAXLINEALONG],s1[50],s2[50];
char **scrap2;

int nrecords=0,monitor=0;
char parameterinwatch[]="beginparameters";

if(fgets(line,MAXLINEALONG,fp)!=NULL){
sscanf(line, "%s",s1);
	if((monitor=strncmp(parameterinwatch,s1,15))!=0){
		printf("ill formatted parameters file/n");
		printf("Enter whatever\n");
		scanf("%s",scrap);
		exit(0);	
	}
}
while (monitor==0 && nrecords<=100 && strncmp(parameterinwatch,"beginparameters",15)==0)
	{
		if(fgets(line,MAXLINEALONG,fp)!=NULL){ 
		
		sscanf(line, "%s%s",s1,s2);
		if (strncmp(s1,"libraryfilename",15)==0) strcpy(library,s2);
		else if (strncmp(s1,"libsize",7)==0) libsize=atol(s2);
		else if (strncmp(s1,"listofpdbs",10)==0) strcpy(pdbslist,s2);
		else if (strncmp(s1,"writepdbfile",12)==0) wpdb=atoi(s2);
		else if (strncmp(s1,"rmsdcutof",9)==0) rmsdcutof1=strtod(s2,scrap2);
		else if (strncmp(s1,"endparameters",13)==0) strcpy(parameterinwatch,s1);
		else {
		printf("ill formatted parameters file/n");
		printf("Enter whatever\n");
		scanf("%s",scrap);
		exit(0);		
		}
		}
		nrecords++;
	}
fprintf(fz,"Run Parameters\n");
fprintf(fz,"libraryfilename=\t%s\n",library);
fprintf(fz,"libsize=\t%li\n",libsize);
fprintf(fz,"listofpdbs=\t%s\n",pdbslist);
fprintf(fz,"rmsdcutof=\t%f\n",rmsdcutof1);
fprintf(fz,"writepdb=\t%i\n",wpdb);

fclose(fp);
fclose(fz);

return 0;	
}
