#include <stdio.h>
#include <stdlib.h>

#include <string.h>
#include <iostream>

#include <mpi.h>
#include <math.h>

#include "vector3d.h"
#include "savebmp.h"
#include "properties.h"

#define epsilon 0.000000000000000222

int main(int argc, char* argv[]){
	
	if( argc != 10){
		printf("Usage: %s numParticlesLight numParticleMedium numParticleHeavy numSteps subSteps timeSubStep imageWidth imageHeight imageFilenamePrex\n", argv[0]);
	}

	MPI_Init(&argc,&argv);

	int p, my_rank;

	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &p);

	//variables
	int numParticlesLight = atof(argv[1]);
	int numParticleMedium = atof(argv[2]);
	int numParticleHeavy = atof(argv[3]);

	int numSteps = atof(argv[4]);
	int subSteps = atof(argv[5]);
	double timeSubStep =atof(argv[6]);

	int width =atof(argv[7]);
	int height =atof(argv[8]);

	unsigned char* image;

	double time, time2;
	int totalNumParticles=numParticlesLight+numParticleMedium+numParticleHeavy;

	//printf("testung");
	//Personal 
	double xdiff;
	double ydiff;
	double dist;
	double den;
	double size;

	//Constants
	double G=0.00000000006673;

	//arrays

	double *positionx = (double*)malloc(sizeof(double)*totalNumParticles);
	double *positiony = (double*)malloc(sizeof(double)*totalNumParticles);
	double *velocitiesx = (double*)malloc(sizeof(double)*totalNumParticles);
	double *velocitiesy = (double*)malloc(sizeof(double)*totalNumParticles);
	double *masses = (double*)malloc(sizeof(double)*totalNumParticles);
	double *ForcesX = (double*)malloc(sizeof(double)*totalNumParticles);
	double *ForcesY = (double*)malloc(sizeof(double)*totalNumParticles);
	double *ForcesIX = (double*)malloc(sizeof(double)*totalNumParticles);
	double *ForcesIY = (double*)malloc(sizeof(double)*totalNumParticles);
	double *ForcesKX = (double*)malloc(sizeof(double)*totalNumParticles);
	double *ForcesKY = (double*)malloc(sizeof(double)*totalNumParticles);

	for(int i=0; i<totalNumParticles; i++){
		positionx[i]=0;
		positiony[i]=0;
		velocitiesx[i]=0;
		velocitiesy[i]=0;
		masses[i]=0;
		ForcesX[i]=0;
		ForcesY[i]=0;
		ForcesIX[i]=0;
		ForcesIY[i]=0;
		ForcesKX[i]=0;
		ForcesKY[i]=0;
	}


	//printf("testung");
	for(int i=0; i<numParticlesLight; i++){
		masses[i]=drand48()+3;
		velocitiesx[i]=drand48()+12;
		velocitiesy[i]=drand48()+12;
		printf("full list %lf %lf %lf %lf %lf %lf %d \n", positionx[i],positiony[i],velocitiesx[i],velocitiesx[i],ForcesX[i],ForcesY[i],i);
	}
	for(int i=numParticlesLight; i<numParticleMedium+numParticlesLight; i++){
		masses[i]=drand48()+8;
		velocitiesx[i]=drand48()+8;
		velocitiesy[i]=drand48()+3;
		printf("full list %lf %lf %lf %lf %lf %lf %d \n", positionx[i],positiony[i],velocitiesx[i],velocitiesx[i],ForcesX[i],ForcesY[i],i);
	}
	for(int i=numParticleMedium+numParticlesLight; i<totalNumParticles; i++){
		masses[i]=drand48()+12;
		velocitiesx[i]=drand48()+3;
		velocitiesy[i]=drand48()+3;
		printf("full list %lf %lf %lf %lf %lf %lf %d \n", positionx[i],positiony[i],velocitiesx[i],velocitiesx[i],ForcesX[i],ForcesY[i],i);
	}
	for(int i=0; i<totalNumParticles; i++){
		positionx[i]=i;
		positiony[i]=i;
		printf("full list %lf %lf %lf %lf %lf %lf %d \n", positionx[i],positiony[i],velocitiesx[i],velocitiesx[i],ForcesX[i],ForcesY[i],i);
	}


	time=MPI_Wtime();

	double *loc_VeloX = (double*)malloc(sizeof(double)*size);
	double *loc_VeloY = (double*)malloc(sizeof(double)*size);
//  need two if else one for first one for gather 
//	
	//root node stuff goes here
	if(my_rank == 0){
		for(int i=0; i<timeSubStep; i++){
			//if(i%subSteps==0){
				//printf(" full list%lf %lf %lf %lf %lf %lf \n",positionx[i],positiony[i],velocitiesx[i],velocitiesx[i],ForcesX[i],ForcesY[i]);
			//}
			//change two latter to scalable
			size=totalNumParticles/p;
			//MPI_Bcast(masses, totalNumParticles, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			//MPI_Bcast(positionx, totalNumParticles, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			//MPI_Bcast(positiony, totalNumParticles, MPI_DOUBLE, 0, MPI_COMM_WORLD);

			//MPI_Scatter(velocitiesx,size, MPI_DOUBLE,loc_VeloX,size,MPI_DOUBLE,0,MPI_COMM_WORLD);
			//MPI_Scatter(velocitiesy,size, MPI_DOUBLE,loc_VeloY,size,MPI_DOUBLE,0,MPI_COMM_WORLD);
			for(int i=0; i<totalNumParticles; i++){
				for(int k=0; k<totalNumParticles; k++){
					if(i!=k){
						xdiff=velocitiesx[i]-velocitiesx[k];
						ydiff=velocitiesy[i]-velocitiesy[k];
						dist =sqrt(xdiff*xdiff+ydiff*ydiff);
						den=dist*dist*dist;
						ForcesX[i]+=masses[i]*masses[k]/den*xdiff;
						ForcesY[i]+=masses[i]*masses[k]/den*ydiff;
						printf("forces %lf %lf \n",ForcesX[i],ForcesY[i]);

					 	//ForcesIX[i]+=ForcesY[i];
					 	//ForcesIY[i]+=ForcesY[i];
		 				//ForcesKX[k]-=ForcesX[i];
						//ForcesKY[k]-=ForcesY[i];
					}
				}
				for(int i=0; i<totalNumParticles; i++){
					positionx[i]=timeSubStep*velocitiesx[i];
					positiony[i]=timeSubStep*velocitiesy[i];
					velocitiesx[i]=timeSubStep/(masses[i]*ForcesX[i]);
					velocitiesy[i]=timeSubStep/(masses[i]*ForcesY[i]);
					printf(" position %lf %lf %lf %lf \n",positionx[i],positiony[i],velocitiesx[i],velocitiesx[i]);
				}
				printf("full list %lf %lf %lf %lf %lf %lf \n", positionx[i],positiony[i],velocitiesx[i],velocitiesx[i],ForcesX[i],ForcesY[i]);
			}

			double *loc_posX = (double*)malloc(sizeof(double)*size);
			double *loc_posY = (double*)malloc(sizeof(double)*size);
		    //MPI_Allgather(loc_posX,size, MPI_DOUBLE, positionx, size, MPI_DOUBLE, MPI_COMM_WORLD);
			//MPI_Allgather(loc_posY,size, MPI_DOUBLE, positiony, size, MPI_DOUBLE, MPI_COMM_WORLD);
		} 
		//almost done, just save the image
		//saveBMP(argv[9], image, width, height);
		time2=MPI_Wtime();
	//}

	free(image);

	MPI_Finalize();
	return 0;
}