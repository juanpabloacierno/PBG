/*
Writes 9-Gly PDB files
*/

void writePDB(char *input, int numModel, char combination[8], double transformed_coord[TOTAL_ATOM_NUMBER][COL_XYZ], int k)
{
char atom[]={'N','C','A','C','O',' '}, pdb_out[50], output_pdb[50],kchar[50];
int i,j, numatom=1,numres=1;

FILE *outfile_pdb;

strcpy(output_pdb,input);
strncat(output_pdb,".",1);

if(k<10){
sprintf(kchar,"00%d",k);
strncat(output_pdb, kchar,3);
}
else if(k>=10 && k<100){
sprintf(kchar,"0%d",k);
strncat(output_pdb, kchar,3);
}
else if(k>=100){
sprintf(kchar,"%d",k);
strncat(output_pdb, kchar,3);
}

strcpy(pdb_out,output_pdb);
strncat(pdb_out,".pdb",4);

outfile_pdb = fopen(pdb_out,"a+");
	if (outfile_pdb==NULL)
		{
		printf("Can't write output PDB.\n");
		exit(666);
		}

fprintf(outfile_pdb,"MODEL%10i\nCOMPND    ",numModel);

for(i=0;i<7;i++)
    {
        fprintf(outfile_pdb,"%c",combination[i]);
    }

fprintf(outfile_pdb,"\n");
numatom=1;
numres=1;

for(i=0;i<TOTAL_ATOM_NUMBER;i++)
    {
		if(numatom<4)numres=1;
			else if(numatom>=5&&numatom<=8)numres=2;
			else if(numatom>=9&&numatom<=12)numres=3;
			else if(numatom>=13&&numatom<=16)numres=4;
			else if(numatom>=17&&numatom<=20)numres=5;
			else if(numatom>=21&&numatom<=24)numres=6;
			else if(numatom>=25&&numatom<=28)numres=7;
            else if(numatom>=29&&numatom<=32)numres=8;
			else if(numatom>=33&&numatom<=36)numres=9;

		if(((numatom+4)%4)==1 )
            {
				fprintf(outfile_pdb,"ATOM     %2i  %c%c  GLY    %2i    ",numatom++,atom[0],atom[5],numres);
			}
			else if(((numatom+4)%4)==2 )
				{
					fprintf(outfile_pdb,"ATOM     %2i  %c%c  GLY    %2i    ",numatom++,atom[1],atom[2],numres);
				}
			else if(((numatom+4)%4)==3 )
				{
					fprintf(outfile_pdb,"ATOM     %2i  %c%c  GLY    %2i    ",numatom++,atom[1],atom[5],numres);
				}
			else if(((numatom+4)%4)==0 )
				{
					fprintf(outfile_pdb,"ATOM     %2i  %c%c  GLY    %2i    ",numatom++,atom[4],atom[5],numres);
				}

        for(j=0;j<3;j++)
            {
                fprintf(outfile_pdb,"%8.3f",transformed_coord[i][j]);
            }

        fprintf(outfile_pdb,"  1.00  0.00");

        if(((numatom+3)%4)==1 )
            {
                fprintf(outfile_pdb,"%12c\n",atom[0]);
            }
        else if(((numatom+3)%4)==2 )
			{
                fprintf(outfile_pdb,"%12c\n",atom[1]);
            }
        else if(((numatom+3)%4)==3 )
			{
				fprintf(outfile_pdb,"%12c\n",atom[1]);
			}
		else if(((numatom+3)%4)==0 )
			{
                fprintf(outfile_pdb,"%12c\n",atom[4]);
            }
    }

fprintf(outfile_pdb,"TER      %2i      GLU     9\nENDMDL\n",numatom-1);

fclose(outfile_pdb);
}
