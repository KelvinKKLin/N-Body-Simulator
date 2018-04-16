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
double massLightMin  = 300000000000000000;
double massLightMax  = 599999999999999999;

double massMediumMin = 600000000000000000;
double massMediumMax = 999999999999999999;

double massHeavyMin  = 100000000000000000;
double massHeavyMax  = 139999999999999999;


//colours
vec3 colourLight = vec3(0,0,1);
vec3 colourMedium = vec3(0,1,0);
vec3 colourHeavy = vec3(1,0,0);

#endif
