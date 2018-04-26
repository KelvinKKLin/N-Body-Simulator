/**
 * main.cpp
 *
 * This file implements a parallel version of the n-body solver.
 * 
 * Authors: Kelvin Lin, Prabhbir Pooni, Yanting Zhang
 *
 **/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <math.h>
#include <omp.h>
#include <float.h>
#include "vector3d.h"
#include "savebmp.h"
#include "properties.h"

#define epsilon 0.000000000000000222

int main(int argc, char* argv[]){
	
	//Check desired number of arguments
	if( argc != 10){
		printf("Usage: %s numParticlesLight numParticleMedium numParticleHeavy numSteps subSteps timeSubStep imageWidth imageHeight imageFilenamePrefix\n", argv[0]);
		return -1;
	}

	//Set up OpenMP
	//omp_set_num_threads(4);

	//Initialize MPI
	MPI_Init(&argc,&argv);
	int p, my_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &p);

	//Extract and initialize variables for the solver
	int numParticleLight = atoi(argv[1]);
	int numParticleMedium = atoi(argv[2]);
	int numParticleHeavy = atoi(argv[3]);
	int n = numParticleLight + numParticleMedium + numParticleHeavy;
	double timeOne, timeTwo, timeFinal;
	
	double* forcex = (double*) calloc(sizeof(double) * n, sizeof(double));
	double* forcey = (double*) calloc(sizeof(double) * n, sizeof(double));
	double* mass = (double*) calloc(sizeof(double) * n, sizeof(double));
	double* velx = (double*) calloc(sizeof(double) * n, sizeof(double));
	double* vely = (double*) calloc(sizeof(double) * n, sizeof(double));
	double* x 	 = (double*) calloc(sizeof(double) * n, sizeof(double));
	double* y	 = (double*) calloc(sizeof(double) * n, sizeof(double));

	int steps = 0;
	int numSteps = atoi(argv[4]) + 1;
	int subSteps = atoi(argv[5]);
	double timeSubStep = atof(argv[6]);
	int width = atoi(argv[7]);
	int height = atoi(argv[8]);

	double G=0.673;

	unsigned char* image = NULL;

	//Randomly assign initial values with drand48()
	if(my_rank == 0){
		//#pragma omp parallel for
		for(int i = 0; i <=numParticleLight; ++i){
			mass[i] = massLightMin + (drand48() * ((massLightMax - massLightMin)+1));
			velx[i] = velocityLightMin + (drand48() * ((velocityLightMax - velocityLightMin)+1));
			vely[i] = velocityLightMin + (drand48() * ((velocityLightMax - velocityLightMin)+1));
			x[i] = drand48() * width;
			y[i] = drand48() * height;
		}
		//#pragma omp parallel for
		for(int i = numParticleLight; i <=numParticleLight + numParticleMedium; ++i){
			mass[i] = massMediumMin + (drand48() * ((massMediumMax - massMediumMin)+1));
			velx[i] = velocityMediumMin + (drand48() * ((velocityMediumMax - velocityMediumMin)+1));
			vely[i] = velocityMediumMin + (drand48() * ((velocityMediumMax - velocityMediumMin)+1));
			x[i] = drand48() * width;
			y[i] = drand48() * height;
		}
		//#pragma omp parallel for
		for(int i = numParticleLight + numParticleMedium; i <=n; ++i){
			mass[i] = massHeavyMin + (drand48() * ((massHeavyMax - massHeavyMin)+1));
			velx[i] = velocityHeavyMin + (drand48() * ((velocityHeavyMax - velocityHeavyMin)+1));
			vely[i] = velocityHeavyMin + (drand48() * ((velocityHeavyMax - velocityHeavyMin)+1));
			x[i] = drand48() * width;
			y[i] = drand48() * height;
		}

	}

	double* localvelx = (double*) calloc(sizeof(double) * n/p, sizeof(double));
	double* localvely = (double*) calloc(sizeof(double) * n/p, sizeof(double));
	double* locposx   = (double*) calloc(sizeof(double) * n/p, sizeof(double));
	double* locposy   = (double*) calloc(sizeof(double) * n/p, sizeof(double));

	int loc_n = n/p;

	image = (unsigned char *) calloc(sizeof(unsigned char)* 3 * width * height, sizeof(unsigned char));		

	//Broadcast variables to the individual nodes
	MPI_Bcast(mass, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(x, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(y, n, MPI_DOUBLE, 0,  MPI_COMM_WORLD);
	MPI_Bcast(velx, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(vely, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	int ultimate_counter = 0;

	//Compute the times
	double min_time = DBL_MAX;
	double max_time = DBL_MIN;
	double average_time = 0;

	//Core iterative process
	for(steps = 0; steps < numSteps*subSteps; steps++){
		//Save the image at each major timestep
		if(steps % subSteps == 0 && my_rank == 0){
			//Create the image
			//#pragma omp parallel for
			for(int i = 0; i <n; i++){
				if (i < numParticleLight){
	                image[((int)x[i] + width*(int)y[i])*3] =  0;
	                image[((int)x[i] + width*(int)y[i])*3+1] = 0;
	                image[((int)x[i] + width*(int)y[i])*3+2] = 255;
                } else if (i >= numParticleLight && i < numParticleLight+numParticleMedium){
                	int count = 0;
                	while(image[((int)x[i] + width*(int)y[i])*3+2] == 255){
                		++count;
                		x[i] = (int)(x[i] + 1)%width;
                		y[i] = (int)(y[i] + 1)%height;
                		if(count == width) break;
                	}
	                image[((int)x[i] + width*(int)y[i])*3] =  0;
	                image[((int)x[i] + width*(int)y[i])*3+1] = 255;
	                image[((int)x[i] + width*(int)y[i])*3+2] = 0;
                } else if(i > numParticleLight + numParticleMedium){
                	int count = 0;
                	while(image[((int)x[i] + width*(int)y[i])*3+2] == 255 || image[((int)x[i] + width*(int)y[i])*3+1] == 255){
                		++count;
                		x[i] = (int)(x[i] + 1)%width;
                		y[i] = (int)(y[i] + 1)%height;
                		if(count == width) break;
                	}
	                image[((int)x[i] + width*(int)y[i])*3] =  255;
	                image[((int)x[i] + width*(int)y[i])*3+1] = 0;
	                image[((int)x[i] + width*(int)y[i])*3+2] = 0;
                }
			}

			//Save the image
			char integer_string[32];
			char filename[64];
			sprintf(integer_string, "%d", steps);
			sprintf(filename, argv[9]);
			strcat(filename, integer_string);
			saveBMP(filename, image, width, height);
			printf("Saving picture #%d\n", ultimate_counter++);

			//Reset the image
			//#pragma omp parallel for
            for(int i = 0; i <n; i++){
				if (i < numParticleLight){
	                image[((int)x[i] + width*(int)y[i])*3] =  0;
	                image[((int)x[i] + width*(int)y[i])*3+1] = 0;
	                image[((int)x[i] + width*(int)y[i])*3+2] = 0;
                } else if (i >= numParticleLight && i < numParticleLight+numParticleMedium){
	                image[((int)x[i] + width*(int)y[i])*3] =  0;
	                image[((int)x[i] + width*(int)y[i])*3+1] = 0;
	                image[((int)x[i] + width*(int)y[i])*3+2] = 0;
                } else{
	                image[((int)x[i] + width*(int)y[i])*3] =  0;
	                image[((int)x[i] + width*(int)y[i])*3+1] = 0;
	                image[((int)x[i] + width*(int)y[i])*3+2] = 0;
                }

			}
		}

		//Start the timer
		if(my_rank==0){
			timeOne = MPI_Wtime();
		}

		//Compute the forces on each particle
		//#pragma omp parallel for
		for(int i = 0; i < n; i++){
			forcex[i] = 0;
			forcey[i] = 0;
			//#pragma omp parallel for
			for(int j = 0; j < n; j++){
				if(i != j){
					double x_diff = (x[j] - x[i]);
					double y_diff = (y[j] - y[i]);
					double dist = sqrt(x_diff*x_diff + y_diff*y_diff);
					double force = (G * mass[i] * mass[j]) / (dist*dist);
					forcex[i] += force * x_diff / dist;
					forcey[i] += force * y_diff / dist;

					//Check if forcex[i] is not a number, if so, reset it to 1
					if((forcex[i])!=forcex[i]){
						forcex[i] = 1;
					}
					if(forcey[i]!=forcey[i]){
						forcey[i] = 1;
					}

				}
			}
		}

		//Compute position and velocity
		//#pragma omp parallel for
		for(int i = 1; i <= loc_n; i++){
			localvelx[i] = velx[i*my_rank+i-1] + timeSubStep * forcex[i*my_rank+i-1] / mass[i];
			localvely[i] = vely[i*my_rank+i-1] + timeSubStep * forcey[i*my_rank+i-1] / mass[i];

			//Check if localvel[x] and localvely[i] are nan. If so, reset to 1.
			if((localvelx[i])!=localvelx[i]){
				localvelx[i] = 1;
			}
			if((localvely[i])!=localvely[i]){
				localvely[i] = 1;
			}

			//Warp the particle if it goes off the screen
			double distXToTravel = x[i*my_rank+i-1] + timeSubStep * localvelx[i];
			double distYToTravel = y[i*my_rank+i-1] + timeSubStep * localvely[i];
			int newx = (int) distXToTravel;
			int newy = (int) distYToTravel;

			if(newx < 0){
				newx *= -1;
			} 
			if(newy < 0){
				newy *= -1;
			}

			locposx[i] = newx % width;
			locposy[i] = newy % height;
		}

		//Collect and broadcast the time
		if(my_rank==0){
		   timeTwo = MPI_Wtime();
	       timeFinal = timeTwo-timeOne;
		   average_time += timeFinal;
		   if(timeFinal > max_time){
		   		max_time = timeFinal;
		   }
		   if(timeFinal < min_time){
		   		min_time = timeFinal;
		   }
		}

		//Wait until all nodes are done, then gather and send everything to every node
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Allgather(locposx, loc_n, MPI_DOUBLE, x, loc_n, MPI_DOUBLE, MPI_COMM_WORLD);
		MPI_Allgather(locposy, loc_n, MPI_DOUBLE, y, loc_n, MPI_DOUBLE, MPI_COMM_WORLD);
		MPI_Allgather(localvelx, loc_n, MPI_DOUBLE, velx, loc_n, MPI_DOUBLE, MPI_COMM_WORLD);
		MPI_Allgather(localvely, loc_n, MPI_DOUBLE, vely, loc_n, MPI_DOUBLE, MPI_COMM_WORLD);
	}
	
	//Wait until all nodes are done
	MPI_Barrier(MPI_COMM_WORLD);

	

	//Free the pointers
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

	//Finalize MPI
	MPI_Finalize();

	printf("%f %f %f\n", min_time, max_time, average_time/(numSteps*subSteps));

	//Return success
	return 0;
}
