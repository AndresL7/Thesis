#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#include "GalCos_variables.h"



char *data_prefix, *catalog_prefix;

int Npart_Total,MINIMUM_MEMBERS,NGB_MAX;
int FLAG_SUBFIND;
int MINIMUM_NSUBSTRUCT,NCLUSTERS, Npart_clustered, UNCLUSTERED;

double BoxSize,REDSHIFT,OMEGA_MATTER,OMEGALAMBDA,HUBBLEPARAM,OMEGABARYON;
int FlagSfr,FlagFeedback,FlagCooling;
int NumFiles;
double COSMIC_TIME;
float Llenght,b_Link;
double GRAV_SOFT,PARTMASS;

double G_INTERNAL_UNITS,LENGHT_INTERNAL_UNITS;
double VELOCITY_INTERNAL_UNITS,MASS_INTERNAL_UNITS,TIME_INTERNAL_UNITS;
double ENERGY_INTERNAL_UNITS,DENSITY_INTERNAL_UNITS,HUBBLE_INTERNAL_UNITS;

struct part *Particle=NULL;
struct halo *Halos=NULL;
