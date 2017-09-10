#include <stdio.h>
#include <stdlib.h>
#include <cstdlib>
#include <math.h>
#include <limits>
#include <string.h>
#include <iostream>
#include <stdarg.h>
#include <unistd.h>
#include <errno.h>
#include "lattice_variables.h"
#include "lattice.h"
#include "stdint.h"


int bond_mat[LEN][LEN];
int num_bonds = 0;
int dihedrals[MAXDIHEDRALS][4];
int cnt_dihedrals=0;
int particle_id[LEN];
/*--------------------------------------------------------------------------------------*/
/*      Function for setting up Initial lattice coordinates                             */
/*--------------------------------------------------------------------------------------*/
int initialLatticeStruct(latticeStruct *ip,  int size)
{
  int k = 0;/*keeping track of coords along x-axis*/
  int j;
  for (int i = 0; i < size; i++)
  {
        j = i/NX;
        ip->x[i] = (j%2 == 0) ? (double)(k * a) : (double)(k * a + 0.5);
        k++;
        k = (k%NX == 0) ? 0 : k; /*resetting k to 0 once the lattice has been traversed in the x direction*/
        ip->y[i] = sqrt(3) * a * 0.5 * j;
        ip->z[i]=0; /*  Starting with flat configuration        */
        //ip->z[NX/2]=grandom(0,ZRAND);//Moving one node randomly
   }
   return 0;
}

/*----------------------------------------------------------------------------------------*/
/*Function for generating the Bond Matrix,Matrix[i][j] = 1,Bond between particles i and j */
/*----------------------------------------------------------------------------------------*/

int lattice_connectivity()
{

  /* Initializing elements of the Bond Matrix to 0	*/
  for(int i=0;i<LEN;i++)
  {
	for(int j=0;j<LEN;j++)
	{
		bond_mat[i][j]=0;
	}
  }
		
  for(int i=0;i<LEN;i++)
  {
     /*	EVEN Row	*/
     if((i/NX)%2==0)
     {
        /*	Left boundary lattice sites	*/
        if(i%NX==0)
        {
	   if(i+NX < LEN)
	      bond_mat[i][i+NX]=1;
	   if(i-NX >= 0)
	      bond_mat[i][i-NX]=1;
	   bond_mat[i][i+1]=1;
	}

	/*	Right boundary lattice sites     */
	else if(i%NX==NX-1)
	{
	   if(i+NX<LEN)
	      bond_mat[i][i+NX]=1;
	   if(i+NX-1<LEN)
	      bond_mat[i][i+NX-1]=1;
	   if(i-NX-1>=0)
	      bond_mat[i][i-NX-1]=1;
	   if(i-NX>=0)
	      bond_mat[i][i-NX]=1;
	   bond_mat[i][i-1]=1;
	}

        /*	Intermediate lattice sites	*/
        else
        {
           if(i+NX<LEN)
              bond_mat[i][i+NX]=1;
           if(i+NX-1<LEN)
              bond_mat[i][i+NX-1]=1;
           bond_mat[i][i-1]=1;
           bond_mat[i][i+1]=1;
           if(i-NX-1>=0)
	       bond_mat[i][i-NX-1]=1;
	   if(i-NX>=0)
	       bond_mat[i][i-NX]=1;
	}
      }

      /*  ODD Row   */
      else
      {
         /*	Left boundary lattice sites     */
	 if(i%NX==0)
	 {
	    if(i+NX+1<LEN)
	       bond_mat[i][i+NX+1]=1;
	    if(i+NX<LEN)
	       bond_mat[i][i+NX]=1;
	    bond_mat[i][i+1]=1;
	    if(i-NX>=0)
	       bond_mat[i][i-NX]=1;
	    if(i-NX+1>=0)
	       bond_mat[i][i-NX+1]=1;
	  }

	  /* Right boundary lattice sites     */
	  else if(i%NX==NX-1)
	  {
	     if(i+NX < LEN)
	        bond_mat[i][i+NX]=1;
             if(i-NX >= 0)
		bond_mat[i][i-NX]=1;
	     bond_mat[i][i-1]=1;
	  }

	  /*      Intermediate lattice sites      */
	  else
	  {
              if(i+NX+1<LEN)
                 bond_mat[i][i+NX+1]=1;
              if(i+NX<LEN)
                 bond_mat[i][i+NX]=1;
              bond_mat[i][i-1]=1;
              bond_mat[i][i+1]=1;
              if(i-NX>=0)
                 bond_mat[i][i-NX]=1;
              if(i-NX+1>=0)
                 bond_mat[i][i-NX+1]=1;
           }
        }							           
     }
  return 0;
}

/*	Bond Matrix is Symmetric	*/
/*	Function to check the same	*/

int check_bond_mat()
{
   for(int i=0;i<LEN;i++)
   {
	for(int j=0;j<i;j++)
	{
		if(bond_mat[i][j]!=bond_mat[j][i])
		{
			printf("At row = %d Col = %d\n",i,j);
			print_and_exit("ERROR: Bond Matrix is not Symmetric\n");
		}
		if(bond_mat[i][j] == 1)
		{
                        num_bonds++;
                }
	}
   }
   return 0;
}

/*	Culling bonds inside the hole	*/
int cull_bonds()
{
   for(int i=0;i<LEN;i++)
   {
        for(int j=0;j<i;j++)
        {
		if(particle_id[i]==2 || particle_id[j]==2)
		{
			bond_mat[i][j]=0;//Bond culled
			bond_mat[j][i]=0;	
		}
	}
   } 
   return 0;
}


/*	Printing the Bond Pairs, no repeats	*/
int bonds(FILE *fp)
{
   for(int i=0;i<LEN;i++)
   {
        for(int j=0;j<i;j++)
        {
	       if(bond_mat[i][j] == 1)
 	       {   
			fprintf(fp,"%d,%d\n",i,j);
	       }
        }
   }
   return 0;
}

/*----------------------------------------------*/
/*	Generate dihedrals HOOMD convention	*/
/*----------------------------------------------*/
int generate_dihedrals()
{
   int j=0;
  
   /*      Type I dihedral         */
   for(int i=0;i<LEN-NX;i++)
   {
	/*	EVEN Row	*/
	if((i/NX)%2==0 && i%NX!=NX-1)
	{
		if(particle_id[i]!=2 && particle_id[i+NX]!=2 && particle_id[i+1]!=2 && particle_id[i+NX+1]!=2)
		{	
			dihedrals[cnt_dihedrals][j]=i;j++;
			dihedrals[cnt_dihedrals][j]=i+NX;j++;
			dihedrals[cnt_dihedrals][j]=i+1;j++;
			dihedrals[cnt_dihedrals][j]=i+NX+1;
			cnt_dihedrals++;
			j=0;
		}
	}
	
	/*	ODD Row		*/
	if((i/NX)%2==1 && i%NX!=NX-1 && i%NX!=NX-2)
	{
		if(particle_id[i]!=2 && particle_id[i+NX+1]!=2 && particle_id[i+1]!=2 && particle_id[i+NX+2]!=2)
                {  
                	dihedrals[cnt_dihedrals][j]=i;j++;
               		dihedrals[cnt_dihedrals][j]=i+NX+1;j++;
                	dihedrals[cnt_dihedrals][j]=i+1;j++;
                	dihedrals[cnt_dihedrals][j]=i+NX+2;
                	cnt_dihedrals++;
                	j=0;
		}
        }
   }

   /*      Type II dihedral         */
   for(int i=NX;i<LEN;i++)
   {
        /*      EVEN Row        */
        if((i/NX)%2==0 && i%NX!=NX-1)
        {
		if(particle_id[i]!=2 && particle_id[i-NX]!=2 && particle_id[i+1]!=2 && particle_id[i-NX+1]!=2)
                {
                	dihedrals[cnt_dihedrals][j]=i;j++;
                	dihedrals[cnt_dihedrals][j]=i-NX;j++;
                	dihedrals[cnt_dihedrals][j]=i+1;j++;
                	dihedrals[cnt_dihedrals][j]=i-NX+1;
                	cnt_dihedrals++;
                	j=0;
		}
        }

        /*      ODD Row         */
        if((i/NX)%2==1 && i%NX!=NX-1 && i%NX!=NX-2)
        {
		if(particle_id[i]!=2 && particle_id[i-NX]!=2 && particle_id[i+1]!=2 && particle_id[i-NX+2]!=2)
                {
                	dihedrals[cnt_dihedrals][j]=i;j++;
                	dihedrals[cnt_dihedrals][j]=i-NX+1;j++;
                	dihedrals[cnt_dihedrals][j]=i+1;j++;
                	dihedrals[cnt_dihedrals][j]=i-NX+2;
                	cnt_dihedrals++;
                	j=0;
		}
        }
   }

   /*      Type III dihedral         */
   for(int i=2*NX;i<LEN;i++)
   {
        /*      EVEN Row        */
        if((i/NX)%2==0 && i%NX!=0)
        {
		if(particle_id[i]!=2 && particle_id[i-NX-1]!=2 && particle_id[i-NX]!=2 && particle_id[i-2*NX]!=2)
                {
                	dihedrals[cnt_dihedrals][j]=i;j++;
                	dihedrals[cnt_dihedrals][j]=i-NX-1;j++;
                	dihedrals[cnt_dihedrals][j]=i-NX;j++;
                	dihedrals[cnt_dihedrals][j]=i-2*NX;
                	cnt_dihedrals++;
                	j=0;
		}
        }

        /*      ODD Row         */
        if((i/NX)%2==1 && i%NX!=NX-1)
        {
		if(particle_id[i]!=2 && particle_id[i-NX]!=2 && particle_id[i-NX+1]!=2 && particle_id[i-2*NX]!=2)
                {
                	dihedrals[cnt_dihedrals][j]=i;j++;
                	dihedrals[cnt_dihedrals][j]=i-NX;j++;
                	dihedrals[cnt_dihedrals][j]=i-NX+1;j++;
                	dihedrals[cnt_dihedrals][j]=i-2*NX;
                	cnt_dihedrals++;
                	j=0;
		}
        }
   }

   return 0;
}	



/*		Print dihedrals		*/		
int out_dihedrals(FILE *fp)
{
   for(int i=0;i<cnt_dihedrals;i++)
   {
	fprintf(fp,"%d,%d,%d,%d\n",dihedrals[i][0],dihedrals[i][1],dihedrals[i][2],dihedrals[i][3]);
   }
   return 0;
}

/*	Particle Typeid		*/
/*	Clamped - Group Id 2	*/
/*	Otherwise - GroupId 1   */

int particle_typeid()
{
   for(int i=0;i<LEN;i++)
   {
//	if(i%NX==0)// || i%NX==NX-1)
	//if(i==0 || i==1 || i==NX || i==LEN-1 || i==LEN-2 || i==LEN-NX-1)
	if(i==0 || i==1 || i==NX || i==LEN-1 || i==LEN-2 || i==LEN-NX-1)
		particle_id[i]=1;
	else if(PIN4==1 && (i==NX-1 || i==NX-2 || i==2*NX-1 || i==2*NX-2))
		particle_id[i]=1;
	else if(PIN4==1 && (i==LEN-NX || i==LEN-NX+1 || i==LEN-2*NX || i==LEN-2*NX+1))
                particle_id[i]=1;
	else
		particle_id[i]=0;
   }
   
	if(STRIP_SIZE < NX/2)/* For no hole input stripsize greater than NX/2   */
        {
		printf("Generating Hole\n");
                /*      Generating the inner boundary of the hole, nodes on inner boundary ID is 5      */
                /*      TOP RIGHT CORNER INNER BOUNDARY         */
                int top_right = (LEN-1)-(int(2*STRIP_SIZE/sqrt(3))*NX)-STRIP_SIZE;
                /*      TOP LEFT CORNER INNER BOUNDARY          */
                int top_left = (LEN-NX)-(int(2*STRIP_SIZE/sqrt(3))*NX)+STRIP_SIZE;
                /*      BOTTOM RIGHT CORNER INNER BOUNDARY      */
                int bottom_right = NX-1+(int(2*STRIP_SIZE/sqrt(3))*NX)-STRIP_SIZE;
                /*      BOTTOM LEFT CORNER INNER BOUNDARY       */
                int bottom_left = (int(2*STRIP_SIZE/sqrt(3))*NX)+STRIP_SIZE;

                int i = top_right;
                int j = top_left;
                int k = 0;
                while(i>bottom_right)
                {
                        i = i-NX;
                        j = j-NX;
                        //nodeGroupId[i] = 5;
                        //nodeGroupId[j] = 5;
                        k = i-1;

                        /*      Nodes which are removed have ID 2       */
                        while(k>j && j >= bottom_left+NX)
                                {
                                        particle_id[k]=2;
					printf("%d %d\t",k,particle_id[k]);
                                	k--;
                        	}
                }
        }
   return 0;
}	

/*	Print particles TypeId		*/
int out_typeId(FILE *fp)
{
   for(int i=0;i<LEN;i++)
   {
//	if(particle_id[i]!=2)
	fprintf(fp,"%u\n",particle_id[i]);
   }
   return 0;
}

void print_and_exit(char *format, ...)
{
    va_list list;
    va_start(list,format);
    vfprintf(stderr,format,list);
    va_end(list);
    exit(1);
}

