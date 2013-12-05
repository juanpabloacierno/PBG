/*
Writes 9-Gly PDB files
*/

int WriteFitLibrary(char name[20], int i, int j)

{

char outfilename[20];	
sprintf(outfilename,"%s.%03i.pdb",name,i);

FILE *mh;
mh = fopen(outfilename,"a");

fprintf(mh,"MODEL%8i\n",memberout[j].numModel);
fprintf(mh,"COMPND%8s\n",memberout[j].combination);
fprintf(mh,"ATOM      1  N   GLY     1%12.3f%8.3f%8.3f  1.00  0.00           N\n",memberout[j].coord[0][0],memberout[j].coord[0][1],memberout[j].coord[0][2]);
fprintf(mh,"ATOM      2  CA  GLY     1%12.3f%8.3f%8.3f  1.00  0.00           C\n",memberout[j].coord[1][0],memberout[j].coord[1][1],memberout[j].coord[1][2]);
fprintf(mh,"ATOM      3  C   GLY     1%12.3f%8.3f%8.3f  1.00  0.00           C\n",memberout[j].coord[2][0],memberout[j].coord[2][1],memberout[j].coord[2][2]);
fprintf(mh,"ATOM      4  O   GLY     1%12.3f%8.3f%8.3f  1.00  0.00           O\n",memberout[j].coord[3][0],memberout[j].coord[3][1],memberout[j].coord[3][2]);
fprintf(mh,"ATOM      5  N   GLY     2%12.3f%8.3f%8.3f  1.00  0.00           N\n",memberout[j].coord[4][0],memberout[j].coord[4][1],memberout[j].coord[4][2]);
fprintf(mh,"ATOM      6  CA  GLY     2%12.3f%8.3f%8.3f  1.00  0.00           C\n",memberout[j].coord[5][0],memberout[j].coord[5][1],memberout[j].coord[5][2]);
fprintf(mh,"ATOM      7  C   GLY     2%12.3f%8.3f%8.3f  1.00  0.00           C\n",memberout[j].coord[6][0],memberout[j].coord[6][1],memberout[j].coord[6][2]);
fprintf(mh,"ATOM      8  O   GLY     2%12.3f%8.3f%8.3f  1.00  0.00           O\n",memberout[j].coord[7][0],memberout[j].coord[7][1],memberout[j].coord[7][2]);
fprintf(mh,"ATOM      9  N   GLY     3%12.3f%8.3f%8.3f  1.00  0.00           N\n",memberout[j].coord[8][0],memberout[j].coord[8][1],memberout[j].coord[8][2]);
fprintf(mh,"ATOM     10  CA  GLY     3%12.3f%8.3f%8.3f  1.00  0.00           C\n",memberout[j].coord[9][0],memberout[j].coord[9][1],memberout[j].coord[9][2]);
fprintf(mh,"ATOM     11  C   GLY     3%12.3f%8.3f%8.3f  1.00  0.00           C\n",memberout[j].coord[10][0],memberout[j].coord[10][1],memberout[j].coord[10][2]);
fprintf(mh,"ATOM     12  O   GLY     3%12.3f%8.3f%8.3f  1.00  0.00           O\n",memberout[j].coord[11][0],memberout[j].coord[11][1],memberout[j].coord[11][2]);
fprintf(mh,"ATOM     13  N   GLY     4%12.3f%8.3f%8.3f  1.00  0.00           N\n",memberout[j].coord[12][0],memberout[j].coord[12][1],memberout[j].coord[12][2]);
fprintf(mh,"ATOM     14  CA  GLY     4%12.3f%8.3f%8.3f  1.00  0.00           C\n",memberout[j].coord[13][0],memberout[j].coord[13][1],memberout[j].coord[13][2]);
fprintf(mh,"ATOM     15  C   GLY     4%12.3f%8.3f%8.3f  1.00  0.00           C\n",memberout[j].coord[14][0],memberout[j].coord[14][1],memberout[j].coord[14][2]);
fprintf(mh,"ATOM     16  O   GLY     4%12.3f%8.3f%8.3f  1.00  0.00           O\n",memberout[j].coord[15][0],memberout[j].coord[15][1],memberout[j].coord[15][2]);
fprintf(mh,"ATOM     17  N   GLY     5%12.3f%8.3f%8.3f  1.00  0.00           N\n",memberout[j].coord[16][0],memberout[j].coord[16][1],memberout[j].coord[16][2]);
fprintf(mh,"ATOM     18  CA  GLY     5%12.3f%8.3f%8.3f  1.00  0.00           C\n",memberout[j].coord[17][0],memberout[j].coord[17][1],memberout[j].coord[17][2]);
fprintf(mh,"ATOM     19  C   GLY     5%12.3f%8.3f%8.3f  1.00  0.00           C\n",memberout[j].coord[18][0],memberout[j].coord[18][1],memberout[j].coord[18][2]);
fprintf(mh,"ATOM     20  O   GLY     5%12.3f%8.3f%8.3f  1.00  0.00           O\n",memberout[j].coord[19][0],memberout[j].coord[19][1],memberout[j].coord[19][2]);
fprintf(mh,"ATOM     21  N   GLY     6%12.3f%8.3f%8.3f  1.00  0.00           N\n",memberout[j].coord[20][0],memberout[j].coord[20][1],memberout[j].coord[20][2]);
fprintf(mh,"ATOM     22  CA  GLY     6%12.3f%8.3f%8.3f  1.00  0.00           C\n",memberout[j].coord[21][0],memberout[j].coord[21][1],memberout[j].coord[21][2]);
fprintf(mh,"ATOM     23  C   GLY     6%12.3f%8.3f%8.3f  1.00  0.00           C\n",memberout[j].coord[22][0],memberout[j].coord[22][1],memberout[j].coord[22][2]);
fprintf(mh,"ATOM     24  O   GLY     6%12.3f%8.3f%8.3f  1.00  0.00           O\n",memberout[j].coord[23][0],memberout[j].coord[23][1],memberout[j].coord[23][2]);
fprintf(mh,"ATOM     25  N   GLY     7%12.3f%8.3f%8.3f  1.00  0.00           N\n",memberout[j].coord[24][0],memberout[j].coord[24][1],memberout[j].coord[24][2]);
fprintf(mh,"ATOM     26  CA  GLY     7%12.3f%8.3f%8.3f  1.00  0.00           C\n",memberout[j].coord[25][0],memberout[j].coord[25][1],memberout[j].coord[25][2]);
fprintf(mh,"ATOM     27  C   GLY     7%12.3f%8.3f%8.3f  1.00  0.00           C\n",memberout[j].coord[26][0],memberout[j].coord[26][1],memberout[j].coord[26][2]);
fprintf(mh,"ATOM     28  O   GLY     7%12.3f%8.3f%8.3f  1.00  0.00           O\n",memberout[j].coord[27][0],memberout[j].coord[27][1],memberout[j].coord[27][2]);
fprintf(mh,"ATOM     29  N   GLY     8%12.3f%8.3f%8.3f  1.00  0.00           N\n",memberout[j].coord[28][0],memberout[j].coord[28][1],memberout[j].coord[28][2]);
fprintf(mh,"ATOM     30  CA  GLY     8%12.3f%8.3f%8.3f  1.00  0.00           C\n",memberout[j].coord[29][0],memberout[j].coord[29][1],memberout[j].coord[29][2]);
fprintf(mh,"ATOM     31  C   GLY     8%12.3f%8.3f%8.3f  1.00  0.00           C\n",memberout[j].coord[30][0],memberout[j].coord[30][1],memberout[j].coord[30][2]);
fprintf(mh,"ATOM     32  O   GLY     8%12.3f%8.3f%8.3f  1.00  0.00           O\n",memberout[j].coord[31][0],memberout[j].coord[31][1],memberout[j].coord[31][2]);
fprintf(mh,"ATOM     33  N   GLY     9%12.3f%8.3f%8.3f  1.00  0.00           N\n",memberout[j].coord[32][0],memberout[j].coord[32][1],memberout[j].coord[32][2]);
fprintf(mh,"ATOM     34  CA  GLY     9%12.3f%8.3f%8.3f  1.00  0.00           C\n",memberout[j].coord[33][0],memberout[j].coord[33][1],memberout[j].coord[33][2]);
fprintf(mh,"ATOM     35  C   GLY     9%12.3f%8.3f%8.3f  1.00  0.00           C\n",memberout[j].coord[34][0],memberout[j].coord[34][1],memberout[j].coord[34][2]);
fprintf(mh,"ATOM     36  O   GLY     9%12.3f%8.3f%8.3f  1.00  0.00           O\n",memberout[j].coord[35][0],memberout[j].coord[35][1],memberout[j].coord[35][2]);
fprintf(mh,"TER      36      GLU     9\n");
fprintf(mh,"ENDMDL\n");
fclose(mh);
return 0;
}
