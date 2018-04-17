#ifndef properties_h
#define properties_h

#include "vector3d.h"

//light particles are the fastest
double velocityLightMin = -3000;
double velocityLightMax = 3000;

double velocityMediumMin = -10000;
double velocityMediumMax = 10000;

//heavy particles are the slowest
double velocityHeavyMin = -30000;
double velocityHeavyMax = 30000;

//mass
double massLightMin  = 0.5e30;
double massLightMax  = 1e30;

double massMediumMin = 1e30;
double massMediumMax = 1.5e30;

double massHeavyMin  = 1.5e30;
double massHeavyMax  = 2e30;


//colours
vec3 colourLight = vec3(0,0,1);
vec3 colourMedium = vec3(0,1,0);
vec3 colourHeavy = vec3(1,0,0);

#endif
