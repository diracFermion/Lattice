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
#include "streamCompaction.h"

int NX,NY,RUN,LEN;
double KAPPA;
float position[NMAX*3];

int main( int argc, char **argv )
{

   /*	Files	*/
   FILE *lat,*therm;

   /*      Character array for directory pathname and filename     */
   char filepath[256],thermalpos_file[256];

  switch (argc){
     case 5:
       sscanf(argv[1],"%d",&NX);
       sscanf(argv[2],"%d",&NY);
       sscanf(argv[3],"%lf",&KAPPA);
       sscanf(argv[4],"%d",&RUN);
       break;
     default:
       print_and_exit("Usage: %s NX NY KAPPA RUN\n",
           argv[0]);
   }
   LEN = NX * NY;

   /*      Filepaths for output files      */
   //sprintf(filepath, "../Sim_dump/lattice.dat");
   //printf("Enter path where the lattice.dat file will be written\n");
   //sprintf(filepath,argv[1]);
   for(int run=1;run<=RUN;run++)
   {
	   sprintf(filepath,"../Sim_dump_ribbon/L%d/W%d/k%.1f/r%d/lattice.dat",NX,NY,KAPPA,run);
	   printf("Filename of Lattice Details: %s\n",filepath);
	   lat = fopen(filepath, "w");
	   if (lat == NULL)
	   {
		print_and_exit("Could Not Open File:lattice.dat");
	   }

	   sprintf(thermalpos_file,"../Sim_dump_ribbon/L%d/W%d/k%.1f/r%d/thermalPosFrame.bin",NX,NY,KAPPA,run);
	   printf("thermalposition file : %s\n",thermalpos_file);
	   therm = fopen(thermalpos_file, "rb");
	   if (therm == NULL)
	   {
		print_and_exit("Could Not Open File to write thermalised position data");
	   }
	   /* Reading binary file with thermalized frame position data	*/
	   fread(position,sizeof(position),1,therm);


	   /*      Byte size of structures         */
	   size_t nBytes = sizeof(latticeStruct); 
		   
	   /*      Dynamic memory allocation for coordinates       */
	   latticeStruct *h_coords = (latticeStruct *)malloc(nBytes);

	   /*      Lattice positions	    */
	   //initialLatticeStruct(h_coords, LEN);
	   for(int i=0;i<LEN;i++)
	   {
		h_coords->x[i] = position[3*i];
		h_coords->y[i] = position[3*i+1];
		h_coords->z[i] = position[3*i+2];;
	   }
	  
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
	   /*	Sorting Dihedrals using 2nd particle as pivot	*/
	   insertionSortDihedrals(cnt_dihedrals);

	   /*	Particle Type ID	*/
	   particle_typeid();

	   /*	Printing lattice configuration	*/
	   /*	Total Particles		*/
	   fprintf(lat,"%d\n",LEN);
	   /*	Particle Position in THERMALIZED configuration	*/
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
	   fclose(therm);

	   bond_compaction(LEN);
   }
   return 0;
}

