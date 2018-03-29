#ifndef properties_h
#define properties_h

#include "vector3d.h"

//light particles are the fastest
double velocityLightMin = 11;
double velocityLightMax = 15;

double velocityMediumMin = 6;
double velocityMediumMax = 10;

//heavy particles are the slowest
double velocityHeavyMin = 1;
double velocityHeavyMax = 5;

//mass
double massLightMin  = 3000000000;
double massLightMax  = 5999999999;

double massMediumMin = 6000000000;
double massMediumMax = 9999999999;

double massHeavyMin  = 10000000000;
double massHeavyMax  = 13999999999;


//colours
vec3 colourLight = vec3(0,0,1);
vec3 colourMedium = vec3(0,1,0);
vec3 colourHeavy = vec3(1,0,0);

#endif
