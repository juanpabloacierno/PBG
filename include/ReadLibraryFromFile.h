int readlibraryfromfile()
{
extern char library[];	
long int i,j;
char modelcode[8];
float x,y,z;

FILE *fp;
	fp = fopen(library,"r");		
	
for (i=0;i<=libsize-1;i++){
	fscanf(fp,"%7s",modelcode);
	strcpy(member[i].combination,modelcode);
	member[i].numModel=i;
	for (j=0;j<=35;j++){
	fscanf(fp,"%8f%8f%8f",&x,&y,&z);
	member[i].coord[j][0]=x;
	member[i].coord[j][1]=y;
	member[i].coord[j][2]=z;
	}
}
fclose(fp);
return 0;	
}
