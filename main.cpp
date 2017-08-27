#include <stdio.h>
#include <stdlib.h>
#include <cstdlib>
#include <math.h>
#include <string.h>
#include <time.h>
#include <limits>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>
#include <vector>
#include "lattice_variables.h"
#include "lattice.h"


int main( int argc, char **argv )
{

   /*	Output File	*/
   FILE *lat;

   /*      Character array for directory pathname and filename     */
   char filepath[256];

   /*      Filepaths for output files      */
   sprintf(filepath, "lattice.dat");
   printf("Filename of Lattice Details: %s\n",filepath);
   lat = fopen(filepath, "w");
   if (lat == NULL)
   {
   	print_and_exit("Could Not Open File:lattice.dat");
   }


   /*      Byte size of structures         */
   size_t nBytes = sizeof(latticeStruct); 
   
   /*      Dynamic memory allocation for coordinates       */
   latticeStruct *h_coords = (latticeStruct *)malloc(nBytes);

   /*      Lattice positions	    */
   initialLatticeStruct(h_coords, LEN);
  
   /*      Generating the Bond Matrix	*/
   lattice_connectivity();

/*   for(int i=0;i<LEN;i++)
   {
	for(int j=0; j<LEN;j++)
	{
		printf("%d ",bond_mat[i][j]);
	}
	printf("\n");
   }
*/   
   /*	Checking if Bond Matrix is Symmetric	*/
   check_bond_mat(); 

   /*	Generating the Dihedrals	*/
   generate_dihedrals();

   /*	Particle Type ID	*/
   particle_typeid();

   /*	Printing lattice configuration	*/
   /*	Total Particles		*/
   fprintf(lat,"%d\n",LEN);
   /*	Initial Paerticle Position in Flat configuration	*/
   for(int i=0;i<LEN;i++)
   {
        fprintf(lat,"%.8f,%.8f,%.8f\n",h_coords->x[i],h_coords->y[i],h_coords->z[i]);
   }
   /*   Printing the Bond pairs */
   fprintf(lat,"%d\n",num_bonds);
   bonds(lat);
   /*   Printing the Dihedrals  */
   fprintf(lat,"%d\n",cnt_dihedrals);
   out_dihedrals(lat);
   /*   Printing particle type Ids	*/
   out_typeId(lat);

   /*	Close output file	*/
   fclose(lat);
   return 0;
}

