#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <math.h>
#include "vector3d.h"
#include "savebmp.h"
#include "properties.h"
#define epsilon 0.000000000000000222

int main(int argc, char* argv[]){
	if( argc != 10){
		printf("Usage: %s numParticlesLight numParticleMedium numParticleHeavy numSteps subSteps timeSubStep imageWidth imageHeight imageFilenamePrefix\n", argv[0]);
		return -1;
	}

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

	int numSteps = atoi(argv[4]);
        numSteps = numSteps+1;
	int subSteps = atoi(argv[5]);
	double timeSubStep = atof(argv[6]);
	int width = atoi(argv[7]);
	int height = atoi(argv[8]);

	double G=6.673e-11;

	unsigned char* image = NULL;
	image = (unsigned char *) calloc(sizeof(unsigned char)* 3 * width * height, sizeof(unsigned char));		

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

	for(int steps = 1; steps < numSteps*subSteps; ++steps){
		if(steps % subSteps == 0){
			for(int i = 0; i <n; i++){
				printf("SubStep TIMESTEP: %d PARTICLE %d POSX: %f, POSY: %f VELX: %f VELY: %f FORCEX: %f FORCEY: %f\n", steps, i, x[i], y[i], velx[i], vely[i], forcex[i], forcey[i]);

				//image array 
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

			//user input filename
			char integer_string[32];
			char filename[64];
			sprintf(integer_string, "%d", steps);
			sprintf(filename, argv[9]);
			strcat(filename, integer_string);
			saveBMP(filename, image, width, height);

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

		for(int i = 0; i < n; i++){
			forcex[i] = 0;
			forcey[i] = 0;
			for(int j = 0; j < n; j++){
				if(i != j){
					double EPS = 3E4;
					double x_diff = x[j] - x[i];
					double y_diff = y[j] - y[i];
					double dist = sqrt(x_diff * x_diff + y_diff * y_diff);
					double force = (G * mass[i] * mass[j]) / (1e20 * dist*dist + EPS*EPS);
					forcex[i] += force * x_diff / dist;
					forcey[i] += force * y_diff / dist;
				}
			}
		}

		for(int i = 0; i < n; i++){
			velx[i] += timeSubStep * forcex[i] / mass[i];
			vely[i] += timeSubStep * forcey[i] / mass[i];
			x[i] += timeSubStep * velx[i];
			y[i] += timeSubStep * vely[i];
			x[i] = (int) x[i] % width;
			y[i] = (int) y[i] % height;
			if(x[i] < 0){
				x[i] *= -1;
			}
			if(y[i] < 0){
				y[i] *= -1;
			}
		}


	} //End of iteration loop

	free(image);
	free(forcex);
	free(forcey);
	free(mass);
	free(velx);
	free(vely);
	free(x);
	free(y);

	return 0;
}
