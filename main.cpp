#include <stdio.h>
#include <stdlib.h>

#include <string.h>

#include <mpi.h>
#include <math.h>

#include "vector3d.h"
#include "savebmp.h"
#include "properties.h"

#define epsilon 0.000000000000000222

void printDoubleArray(double* array, int n){
	for(int i = 0; i < n; ++i){
		printf("%f\n", array[i]);
	}
}

void drawImage(int i, unsigned char* image, double* x, double* y, int width, int numParticleLight, int numParticleMedium, int numParticleHeavy){
	if (i < numParticleLight){
		image[((int)x[i] + width*(int)y[i])*3] =  0;
		image[((int)x[i] + width*(int)y[i])*3+1] = 0;
		image[((int)x[i] + width*(int)y[i])*3+2] = 255;
	} else if (i >= numParticleLight && i < numParticleLight+numParticleMedium){
		image[((int)x[i] + width*(int)y[i])*3] =  0;
		image[((int)x[i] + width*(int)y[i])*3+1] = 255;
		image[((int)x[i] + width*(int)y[i])*3+2] = 0;
	} else{
		image[((int)x[i] + width*(int)y[i])*3] =  255;
		image[((int)x[i] + width*(int)y[i])*3+1] = 0;
		image[((int)x[i] + width*(int)y[i])*3+2] = 0;
	}
}

int main(int argc, char* argv[]){
	
	if( argc != 10){
		printf("Usage: %s numParticlesLight numParticleMedium numParticleHeavy numSteps subSteps timeSubStep imageWidth imageHeight imageFilenamePrefix\n", argv[0]);
		return -1;
	}

	MPI_Init(&argc,&argv);
	int p, my_rank;

	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &p);

	//variables
	int numParticleLight = atoi(argv[1]);
	int numParticleMedium = atoi(argv[2]);
	int numParticleHeavy = atoi(argv[3]);
	int n = numParticleLight + numParticleMedium + numParticleHeavy;

	double* forcex = (double*) calloc(sizeof(double) * n, sizeof(double));
	double* forcey = (double*) calloc(sizeof(double) * n, sizeof(double));
	double* mass = (double*) calloc(sizeof(double) * n, sizeof(double));
	double* velx = (double*) calloc(sizeof(double) * n, sizeof(double));
	double* vely = (double*) calloc(sizeof(double) * n, sizeof(double));
		
	double* x 	 = (double*) calloc(sizeof(double) * n, sizeof(double));
	double* y	 = (double*) calloc(sizeof(double) * n, sizeof(double));

	int steps = 0;
	int numSteps = atoi(argv[4]);
	int subSteps = atoi(argv[5]);
	double timeSubStep = atof(argv[6]);

	int width = atoi(argv[7]);
	int height = atoi(argv[8]);

	double G=0.00000000006673;

	unsigned char* image = NULL;

	//root node stuff goes here
	if(my_rank == 0){
		for(int i = 0; i <=numParticleLight; ++i){
			mass[i] = massLightMin + (drand48() * ((massLightMax - massLightMin)+1));
			velx[i] = velocityLightMin + (drand48() * ((velocityLightMax - velocityLightMin)+1));
			vely[i] = velocityLightMin + (drand48() * ((velocityLightMax - velocityLightMin)+1));
			x[i] = 1;//(drand48() * ((width) + 1));
			y[i] = 2;//(drand48() * ((height) + 1));
		}

		for(int i = numParticleLight; i <=numParticleLight + numParticleMedium; ++i){
			mass[i] = massMediumMin + (drand48() * ((massMediumMax - massMediumMin)+1));
			velx[i] = velocityMediumMin + (drand48() * ((velocityMediumMax - velocityMediumMin)+1));
			vely[i] = velocityMediumMin + (drand48() * ((velocityMediumMax - velocityMediumMin)+1));
			x[i] = 3;//(drand48() * ((width) + 1));
			y[i] = 4;//(drand48() * ((height) + 1));
		}

		for(int i = numParticleLight + numParticleMedium; i <=n; ++i){
			mass[i] = massHeavyMin + (drand48() * ((massHeavyMax - massHeavyMin)+1));
			velx[i] = velocityHeavyMin + (drand48() * ((velocityHeavyMax - velocityHeavyMin)+1));
			vely[i] = velocityHeavyMin + (drand48() * ((velocityHeavyMax - velocityHeavyMin)+1));
			x[i] = 5;//(drand48() * ((width) + 1));
			y[i] = 6;//(drand48() * ((height) + 1));
		}

		//almost done, just save the image
		//saveBMP(argv[9], image, width, height);
	}

	double* localvelx = (double*) calloc(sizeof(double) * n/p, sizeof(double));
	double* localvely = (double*) calloc(sizeof(double) * n/p, sizeof(double));
	double* locposx   = (double*) calloc(sizeof(double) * n/p, sizeof(double));
	double* locposy   = (double*) calloc(sizeof(double) * n/p, sizeof(double));

	int loc_n = p/n;

	MPI_Bcast(mass, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(x, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(y, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Scatter(velx, n/p, MPI_DOUBLE, localvelx, n/p, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Scatter(vely, n/p, MPI_DOUBLE, localvely, n/p, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	image = (unsigned char *) calloc(sizeof(unsigned char)* 3 * width * height, sizeof(unsigned char));		

	for(steps = 0; steps < numSteps; steps += timeSubStep){
		if(steps % subSteps == 0 && my_rank == 0){
			for(int i = 0; i <n; i++){
				printf("SubStep TIMESTEP: %d PARTICLE %d POSX: %f, POSY: %f VELX: %f VELY: %f\n", steps, i, x[i], y[i], velx[i], vely[i]);
				drawImage(i, image, x, y, width, numParticleLight, numParticleMedium, numParticleHeavy);
			}
			char integer_string[32];
			char filename[64] = "IMAGE_";
			sprintf(integer_string, "%d", steps);
			strcat(filename, integer_string);
			saveBMP(filename, image, width, height);
		}

		//Compute the forces on each particle
		for(int i = 0; i <loc_n; i++){
			for(int j = 0; j <loc_n; j++){
				if(i != j){
					double x_diff = x[i] - x[j];
					double y_diff = y[i] - y[j];
					double dist = sqrt(x_diff * x_diff + y_diff + y_diff);
					double dist_cubed = dist* dist * dist;
					forcex[i] -= mass[i]*mass[j] * G / dist_cubed * x_diff;
					forcey[i] -= mass[i]*mass[j] * G / dist_cubed * y_diff;
				}
			}

		}

		//Compute position and velocity
		for(int i = 0; i <=loc_n; ++i){
			locposx[i] += timeSubStep * localvelx[i];
			locposy[i] += timeSubStep * localvely[i];
			/// changed velx  to localvelx 
			localvelx[i] += timeSubStep/mass[i] * forcex[i];
			localvely[i] += timeSubStep/mass[i] * forcey[i];
			//printf("PosX: %f PosY: %f VelX: %f VelY: %f \n ", locposx[i], locposy[i], localvelx[i], localvely[i]);
		}

		//All gather
		MPI_Allgather(locposx, loc_n, MPI_DOUBLE, x, loc_n, MPI_DOUBLE, MPI_COMM_WORLD);
		MPI_Allgather(locposy, loc_n, MPI_DOUBLE, y, loc_n, MPI_DOUBLE, MPI_COMM_WORLD);
		// added these statements
		MPI_Allgather(localvelx, loc_n, MPI_DOUBLE, velx, n/p, MPI_DOUBLE, MPI_COMM_WORLD);
		MPI_Allgather(localvely, loc_n, MPI_DOUBLE, vely, n/p, MPI_DOUBLE, MPI_COMM_WORLD);
	}
	
	MPI_Barrier(MPI_COMM_WORLD);

	free(image);
	free(localvelx);
	free(localvely);
	free(locposx);
	free(locposy);
	free(forcex);
	free(forcey);
	free(mass);
	free(velx);
	free(vely);
	free(x);
	free(y);

	MPI_Finalize();
	return 0;
}
