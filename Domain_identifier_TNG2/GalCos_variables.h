#ifndef INIT_GALCOS_H
#define INIT_GALCOS_H

#include<stdlib.h>
#include<stdio.h>
#include<math.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_sort_double.h>

#define EMPTY_FLAG -1
#define SUCCES 0
#define FAILURE 1
#define BH_OPENING 0.7

extern int Npart_Total,MINIMUM_MEMBERS,NGB_MAX;
extern int FLAG_SUBFIND;
extern int MINIMUM_NSUBSTRUCT,NCLUSTERS, UNCLUSTERED; 

extern double BoxSize,REDSHIFT,OMEGA_MATTER,OMEGALAMBDA,HUBBLEPARAM,OMEGABARYON,COSMIC_TIME,PARTMASS;
extern float Llenght,b_Link;
extern double GRAV_SOFT;

extern double G_INTERNAL_UNITS,LENGHT_INTERNAL_UNITS;
extern double VELOCITY_INTERNAL_UNITS,MASS_INTERNAL_UNITS,TIME_INTERNAL_UNITS;
extern double ENERGY_INTERNAL_UNITS,DENSITY_INTERNAL_UNITS,HUBBLE_INTERNAL_UNITS;

extern struct halo
{
  float mass;
  float pos[3];
  float Rvir,Mvir;
  int Nmembers;
  int NDomain_particles;
  int IDcluster;
  int *Halo_particles;
  int *Domain_particles;
}*Halos;

extern struct part
{
  float pos[3];
  int Oid;
  int Cluster_ID;
}*Particle;

float distance(float xi, float yi, float zi, float xj, float yj, float zj);
int gsl_fisort(int dimension,float *fvector,int *ivector);

// Funciones para grid espacial de halos
void build_halo_spatial_grid(void);
void free_halo_spatial_grid(void);

#endif /* INIT_GALCOS_H */
