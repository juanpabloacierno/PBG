/*This is an example of the data in struc libout (the coordinates are before any fit)
lib memberout[0],numModel=1
		[atom]						[X		Y	    Z]
ATOM      1  N   GLY     1      -4.859   0.073   0.122  1.00  0.00           N
ATOM      2  CA  GLY     1      -3.610  -0.654   0.017  1.00  0.00           C
ATOM      3  C   GLY     1      -2.432   0.308   0.088  1.00  0.00           C
ATOM      4  O   GLY     1      -2.621   1.515   0.217  1.00  0.00           O
ATOM      5  N   GLY     2      -1.213  -0.231   0.004  1.00  0.00           N
ATOM      6  CA  GLY     2      -0.016   0.584   0.061  1.00  0.00           C
ATOM      7  C   GLY     2       1.100  -0.061  -0.749  1.00  0.00           C
ATOM      8  O   GLY     2       1.157  -1.283  -0.869  1.00  0.00           O
ATOM      9  N   GLY     3       1.989   0.764  -1.306  1.00  0.00           N
ATOM     10  CA  GLY     3       3.094   0.267  -2.101  1.00  0.00           C
ATOM     11  C   GLY     3       4.286   1.206  -1.991  1.00  0.00           C
ATOM     12  O   GLY     3       4.116   2.409  -1.804  1.00  0.00           O
ATOM     13  N   GLY     4       5.495   0.654  -2.109  1.00  0.00           N
ATOM     14  CA  GLY     4       6.704   1.448  -2.021  1.00  0.00           C
ATOM     15  C   GLY     4       7.804   0.825  -2.868  1.00  0.00           C
ATOM     16  O   GLY     4       7.842  -0.390  -3.047  1.00  0.00           O
ATOM     17  N   GLY     5       8.701   1.663  -3.391  1.00  0.00           N
ATOM     18  CA  GLY     5       9.793   1.188  -4.215  1.00  0.00           C
ATOM     19  C   GLY     5      10.999   2.104  -4.069  1.00  0.00           C
ATOM     20  O   GLY     5      10.848   3.299  -3.823  1.00  0.00           O
ATOM     21  N   GLY     6      12.199   1.540  -4.221  1.00  0.00           N
ATOM     22  CA  GLY     6      13.419   2.311  -4.104  1.00  0.00           C
ATOM     23  C   GLY     6      14.504   1.713  -4.987  1.00  0.00           C
ATOM     24  O   GLY     6      14.523   0.506  -5.223  1.00  0.00           O
ATOM     25  N   GLY     7      15.410   2.562  -5.474  1.00  0.00           N
ATOM     26  CA  GLY     7      16.489   2.110  -6.327  1.00  0.00           C
ATOM     27  C   GLY     7      17.709   3.000  -6.146  1.00  0.00           C
ATOM     28  O   GLY     7      17.578   4.185  -5.842  1.00  0.00           O
ATOM     29  N   GLY     8      18.898   2.426  -6.333  1.00  0.00           N
ATOM     30  CA  GLY     8      20.130   3.173  -6.187  1.00  0.00           C
ATOM     31  C   GLY     8      21.200   2.602  -7.104  1.00  0.00           C
ATOM     32  O   GLY     8      21.199   1.407  -7.398  1.00  0.00           O
ATOM     33  N   GLY     9      22.116   3.459  -7.556  1.00  0.00           N
ATOM     34  CA  GLY     9      23.184   3.033  -8.437  1.00  0.00           C
ATOM     35  C   GLY     9      24.073   4.216  -8.791  1.00  0.00           C
ATOM     36  O   GLY     9      23.836   5.335  -8.343  1.00  0.00           O

En targetdata el 9mer i=1
ATOM	1	1	N	-5.066	   0.058	  13.305

ATOM	34	9	CA	 9.454	  -7.193	   1.726

*/

void C1C2andC8C9dis(int iesimo9mer, int jesimo9mer)
{
double c1[COL_XYZ],c2[COL_XYZ],c8[COL_XYZ],c9[COL_XYZ],bond12, bond89;




		//calculates targetCA1-libraryCA2 distance
		
		c1[0]=targetdata[4*iesimo9mer-3].coord[0];//CA1,x
		c1[1]=targetdata[4*iesimo9mer-3].coord[1];//CA1,y
		c1[2]=targetdata[4*iesimo9mer-3].coord[2];//CA1,z	
		c2[0]=memberout[jesimo9mer].coord[5][0];//CA2,x
		c2[1]=memberout[jesimo9mer].coord[5][1];//CA2,y
		c2[2]=memberout[jesimo9mer].coord[5][2];//CA2,z	
		bond12=sqrtf((c1[0]-c2[0])*(c1[0]-c2[0]) + (c1[1]-c2[1])*(c1[1]-c2[1]) + (c1[2]-c2[2])*(c1[2]-c2[2]));
		memberout[jesimo9mer].targetC1libraryC2dis=bond12;
		
		
		//calculates targetCA9-libraryCA8 distance
		
		c9[0]=targetdata[4*iesimo9mer+29].coord[0];//CA9,x
		c9[1]=targetdata[4*iesimo9mer+29].coord[1];//CA9,y
		c9[2]=targetdata[4*iesimo9mer+29].coord[2];//CA9,z	
		c8[0]=memberout[jesimo9mer].coord[29][0];//CA8,x
		c8[1]=memberout[jesimo9mer].coord[29][1];//CA8,y
		c8[2]=memberout[jesimo9mer].coord[29][2];//CA8,z	
		bond89=sqrtf((c8[0]-c9[0])*(c8[0]-c9[0]) + (c8[1]-c9[1])*(c8[1]-c9[1]) + (c8[2]-c9[2])*(c8[2]-c9[2]));
		memberout[jesimo9mer].targetC9libraryC8dis=bond89;
		


}
