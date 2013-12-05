/*
Writes 9-Gly PDB files
*/

int WriteFitLibrary(char filename[20], long numModel, char code[8], double a[TOTAL_ATOM_NUMBER][COL_XYZ], int k)

{


FILE *mh;
mh = fopen(filename,"a+");

fprintf(mh,"MODEL%6s\n",code);
fprintf(mh,"COMPND%8i\n",numModel);
fprintf(mh,"ATOM      1  N   GLY     1%12.3f%8.3f%8.3f  1.00  0.00           N\n",a[0][0],a[0][1],a[0][2]);
fprintf(mh,"ATOM      2  CA  GLY     1%12.3f%8.3f%8.3f  1.00  0.00           C\n",a[1][0],a[1][1],a[1][2]);
fprintf(mh,"ATOM      3  C   GLY     1%12.3f%8.3f%8.3f  1.00  0.00           C\n",a[2][0],a[2][1],a[2][2]);
fprintf(mh,"ATOM      4  O   GLY     1%12.3f%8.3f%8.3f  1.00  0.00           O\n",a[3][0],a[3][1],a[3][2]);
fprintf(mh,"ATOM      5  N   GLY     2%12.3f%8.3f%8.3f  1.00  0.00           N\n",a[4][0],a[4][1],a[4][2]);
fprintf(mh,"ATOM      6  CA  GLY     2%12.3f%8.3f%8.3f  1.00  0.00           C\n",a[5][0],a[5][1],a[5][2]);
fprintf(mh,"ATOM      7  C   GLY     2%12.3f%8.3f%8.3f  1.00  0.00           C\n",a[6][0],a[6][1],a[6][2]);
fprintf(mh,"ATOM      8  O   GLY     2%12.3f%8.3f%8.3f  1.00  0.00           O\n",a[7][0],a[7][1],a[7][2]);
fprintf(mh,"ATOM      9  N   GLY     3%12.3f%8.3f%8.3f  1.00  0.00           N\n",a[8][0],a[8][1],a[8][2]);
fprintf(mh,"ATOM     10  CA  GLY     3%12.3f%8.3f%8.3f  1.00  0.00           C\n",a[9][0],a[9][1],a[9][2]);
fprintf(mh,"ATOM     11  C   GLY     3%12.3f%8.3f%8.3f  1.00  0.00           C\n",a[10][0],a[10][1],a[10][2]);
fprintf(mh,"ATOM     12  O   GLY     3%12.3f%8.3f%8.3f  1.00  0.00           O\n",a[11][0],a[11][1],a[11][2]);
fprintf(mh,"ATOM     13  N   GLY     4%12.3f%8.3f%8.3f  1.00  0.00           N\n",a[12][0],a[12][1],a[12][2]);
fprintf(mh,"ATOM     14  CA  GLY     4%12.3f%8.3f%8.3f  1.00  0.00           C\n",a[13][0],a[13][1],a[13][2]);
fprintf(mh,"ATOM     15  C   GLY     4%12.3f%8.3f%8.3f  1.00  0.00           C\n",a[14][0],a[14][1],a[14][2]);
fprintf(mh,"ATOM     16  O   GLY     4%12.3f%8.3f%8.3f  1.00  0.00           O\n",a[15][0],a[15][1],a[15][2]);
fprintf(mh,"ATOM     17  N   GLY     5%12.3f%8.3f%8.3f  1.00  0.00           N\n",a[16][0],a[16][1],a[16][2]);
fprintf(mh,"ATOM     18  CA  GLY     5%12.3f%8.3f%8.3f  1.00  0.00           C\n",a[17][0],a[17][1],a[17][2]);
fprintf(mh,"ATOM     19  C   GLY     5%12.3f%8.3f%8.3f  1.00  0.00           C\n",a[18][0],a[18][1],a[18][2]);
fprintf(mh,"ATOM     20  O   GLY     5%12.3f%8.3f%8.3f  1.00  0.00           O\n",a[19][0],a[19][1],a[19][2]);
fprintf(mh,"ATOM     21  N   GLY     6%12.3f%8.3f%8.3f  1.00  0.00           N\n",a[20][0],a[20][1],a[20][2]);
fprintf(mh,"ATOM     22  CA  GLY     6%12.3f%8.3f%8.3f  1.00  0.00           C\n",a[21][0],a[21][1],a[21][2]);
fprintf(mh,"ATOM     23  C   GLY     6%12.3f%8.3f%8.3f  1.00  0.00           C\n",a[22][0],a[22][1],a[22][2]);
fprintf(mh,"ATOM     24  O   GLY     6%12.3f%8.3f%8.3f  1.00  0.00           O\n",a[23][0],a[23][1],a[23][2]);
fprintf(mh,"ATOM     25  N   GLY     7%12.3f%8.3f%8.3f  1.00  0.00           N\n",a[24][0],a[24][1],a[24][2]);
fprintf(mh,"ATOM     26  CA  GLY     7%12.3f%8.3f%8.3f  1.00  0.00           C\n",a[25][0],a[25][1],a[25][2]);
fprintf(mh,"ATOM     27  C   GLY     7%12.3f%8.3f%8.3f  1.00  0.00           C\n",a[26][0],a[26][1],a[26][2]);
fprintf(mh,"ATOM     28  O   GLY     7%12.3f%8.3f%8.3f  1.00  0.00           O\n",a[27][0],a[27][1],a[27][2]);
fprintf(mh,"ATOM     29  N   GLY     8%12.3f%8.3f%8.3f  1.00  0.00           N\n",a[28][0],a[28][1],a[28][2]);
fprintf(mh,"ATOM     30  CA  GLY     8%12.3f%8.3f%8.3f  1.00  0.00           C\n",a[29][0],a[29][1],a[29][2]);
fprintf(mh,"ATOM     31  C   GLY     8%12.3f%8.3f%8.3f  1.00  0.00           C\n",a[30][0],a[30][1],a[30][2]);
fprintf(mh,"ATOM     32  O   GLY     8%12.3f%8.3f%8.3f  1.00  0.00           O\n",a[31][0],a[31][1],a[31][2]);
fprintf(mh,"ATOM     33  N   GLY     9%12.3f%8.3f%8.3f  1.00  0.00           N\n",a[32][0],a[32][1],a[32][2]);
fprintf(mh,"ATOM     34  CA  GLY     9%12.3f%8.3f%8.3f  1.00  0.00           C\n",a[33][0],a[33][1],a[33][2]);
fprintf(mh,"ATOM     35  C   GLY     9%12.3f%8.3f%8.3f  1.00  0.00           C\n",a[34][0],a[34][1],a[34][2]);
fprintf(mh,"ATOM     36  O   GLY     9%12.3f%8.3f%8.3f  1.00  0.00           O\n",a[35][0],a[35][1],a[35][2]);
fprintf(mh,"TER      36      GLU     9\n");
fprintf(mh,"ENDMDL\n");
fclose(mh);
return 0;
}
