#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <math.h>
#include "vector3d.h"
#include "savebmp.h"
#include "properties.h"
#define epsilon 0.000000000000000222

using namespace std;

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
	double timeOne, timeTwo, timeFinal;
	
	double* forcex = (double*) calloc(sizeof(double) * n, sizeof(double));
	double* forcey = (double*) calloc(sizeof(double) * n, sizeof(double));
	double* mass = (double*) calloc(sizeof(double) * n, sizeof(double));
	double* velx = (double*) calloc(sizeof(double) * n, sizeof(double));
	double* vely = (double*) calloc(sizeof(double) * n, sizeof(double));
	double* x 	 = (double*) calloc(sizeof(double) * n, sizeof(double));
	double* y	 = (double*) calloc(sizeof(double) * n, sizeof(double));

	int steps = 0;
	int numSteps = atoi(argv[4]);
        numSteps = numSteps+1;
	int subSteps = atoi(argv[5]);
	double timeSubStep = atof(argv[6]);
	int width = atoi(argv[7]);
	int height = atoi(argv[8]);

	double G=0.673;

	unsigned char* image = NULL;

	//root node stuff goes here
	if(my_rank == 0){
		for(int i = 0; i <=numParticleLight; ++i){
			mass[i] = massLightMin + (drand48() * ((massLightMax - massLightMin)+1));
			velx[i] = velocityLightMin + (drand48() * ((velocityLightMax - velocityLightMin)+1));
			vely[i] = velocityLightMin + (drand48() * ((velocityLightMax - velocityLightMin)+1));
			x[i] = drand48() * width;
			y[i] = drand48() * height;
		}

		for(int i = numParticleLight; i <=numParticleLight + numParticleMedium; ++i){
			mass[i] = massMediumMin + (drand48() * ((massMediumMax - massMediumMin)+1));
			velx[i] = velocityMediumMin + (drand48() * ((velocityMediumMax - velocityMediumMin)+1));
			vely[i] = velocityMediumMin + (drand48() * ((velocityMediumMax - velocityMediumMin)+1));
			x[i] = drand48() * width;
			y[i] = drand48() * height;
		}

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
	
	if(my_rank==0){
		timeOne = MPI_Wtime();
	}

	MPI_Bcast(mass, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(x, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(y, n, MPI_DOUBLE, 0,  MPI_COMM_WORLD);
	MPI_Bcast(velx, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(vely, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	int ultimate_counter = 0;

	for(steps = 0; steps < numSteps*subSteps; steps++){
		if(steps % subSteps == 0 && my_rank == 0){
			for(int i = 0; i <n; i++){
				//printf("SubStep TIMESTEP: %d PARTICLE %d POSX: %f, POSY: %f VELX: %f VELY: %f FORCEX: %f FORCEY: %f\n", steps, i, x[i], y[i], velx[i], vely[i], forcex[i], forcey[i]);
				//cout << "Substep Timestep "<< steps << " Particle " << i << " Pos X " << x[i] << " Pos Y " << y[i] << " VelX " << velx[i] << " VelY " << vely[i] << " ForceX " << forcex[i] << " ForceY " << forcey[i] << endl;
				
				//image array 
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

			//user input filename
			char integer_string[32];
			char filename[64];
			sprintf(integer_string, "%d", steps);
			sprintf(filename, argv[9]);
			strcat(filename, integer_string);
			saveBMP(filename, image, width, height);
			printf("Saving picture #%d\n", ultimate_counter++);

			//assigning pixels to black again
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

		//Compute the forces on each particle
		for(int i = 0; i < n; i++){
			forcex[i] = 0;
			forcey[i] = 0;
			for(int j = 0; j < n; j++){
				if(i != j){
					double x_diff = (x[j] - x[i]);
					double y_diff = (y[j] - y[i]);
					double dist = sqrt(x_diff*x_diff + y_diff*y_diff);
					double force = (G * mass[i] * mass[j]) / (dist*dist);
					forcex[i] += force * x_diff / dist;
					forcey[i] += force * y_diff / dist;

					//These guards may be removed when warping is complete.
					if(isnan(forcex[i])){
						forcex[i] = 0;
					}
					if(isnan(forcey[i])){
						forcey[i] = 0;
					}

				}
			}

		}

		//Compute position and velocity
		for(int i = 0; i < loc_n; i++){
			localvelx[i] = velx[i] + timeSubStep * forcex[i] / mass[i];
			localvely[i] = vely[i] + timeSubStep * forcey[i] / mass[i];

			//These guards may be removed when warping is complete
			if(isnan(localvelx[i])){
				localvelx[i] = 0;
			}
			if(isnan(localvely[i])){
				localvely[i] = 0;
			}

			//TODO: Implement intelligent warp
			double distXToTravel = x[i] + timeSubStep * localvelx[i];
			double distYToTravel = y[i] + timeSubStep * localvely[i];
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

		//All gather
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Allgather(locposx, loc_n, MPI_DOUBLE, x, loc_n, MPI_DOUBLE, MPI_COMM_WORLD);
		MPI_Allgather(locposy, loc_n, MPI_DOUBLE, y, loc_n, MPI_DOUBLE, MPI_COMM_WORLD);
		MPI_Allgather(localvelx, loc_n, MPI_DOUBLE, velx, loc_n, MPI_DOUBLE, MPI_COMM_WORLD);
		MPI_Allgather(localvely, loc_n, MPI_DOUBLE, vely, loc_n, MPI_DOUBLE, MPI_COMM_WORLD);
	}
	
	MPI_Barrier(MPI_COMM_WORLD);

	if(my_rank==0){
	   timeTwo = MPI_Wtime();
           timeFinal = timeTwo-timeOne;
	   printf("%lf \n", timeFinal);
	}

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
