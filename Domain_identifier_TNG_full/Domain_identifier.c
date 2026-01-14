/* 
   This program takes a snapshot of a simulations and the list of FOF
   halos. With that, the code identifies the domain PARTICLES of every
   FOF halo in that snapshot.
*/

#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<time.h>
#include<mpi.h>

#include "GalCos_variables.h"
#include "GalCos_variables.c"


#define MAX_FILENAME 80

int task, Number_Processes;

// structure to store the domain information of every processor
// Nparts_per_node: number of particles in the node
// istart: first particle in the node
// iend: last particle in the node
struct domain
{
  int Nparts_per_node, istart, iend;
} *domain_info;

void usage_main(void)
{
  
    if(task == 0)
    {
      printf("Domain_identifier = Identifies the domain of FOF halos.\n"); fflush(stdout);
      printf("Juan Carlos Munoz C.\n"); fflush(stdout);
      printf("Usage:./Domain_identifier <Gadget_snapshoth_file> <Gadget_catalog_file> <param_file>\n"); fflush(stdout);
    }
    
    MPI_Finalize();
    exit(0);
}

#include "GalCos_routines.c"
#include "GalCos_load_Gadget.c"

#include "GalCos_units.c"
// #include "GalCos_cola.c"
// #include "GalCos_clusterBound.c"

// #include "GalCos_Halo_prop.c"
#include "GalCos_domain.c"

#include "GalCos_PartitionParticles.c"
#include "GalCos_loadTNG_HaloCatalog.c"

// structure to store the information of every halo
struct data
{
  float mass;
  float pos[3];
  float Rvir;
  int Nmembers;
  int NDomain_particles;
  int IDcluster;
};

// Function to write a rescue file with the information of the halos
// and the particles of the halos in a binary file for backup purposes
// storing the halo mass, position, radius, number of members,
// number of domain particles, and the ID of the cluster
int write_rescue(char *infile)
{
  int i;
  FILE *pf;
  char buf[1000];
  struct data aux;
  
  sprintf(buf,"%s%s",infile,".rescue");
  pf = fopen(buf,"wb");
  
  fwrite(&NCLUSTERS,sizeof(int),1,pf);
  
  for(i=0; i<NCLUSTERS; i++)
    {
      aux.mass = Halos[i].mass;
      aux.pos[0] = Halos[i].pos[0];
      aux.pos[1] = Halos[i].pos[1];
      aux.pos[2] = Halos[i].pos[2];
      aux.Nmembers = Halos[i].Nmembers; // number of particles in the halo	
      aux.NDomain_particles = Halos[i].NDomain_particles; // number of particles in the domain of halo
      aux.IDcluster = Halos[i].IDcluster; 
      aux.Rvir = Halos[i].Rvir;
      
      fwrite(&aux,sizeof(struct data),1,pf);

      if(Halos[i].Nmembers > 0)
	fwrite(&(Halos[i].Halo_particles[0]),sizeof(int),Halos[i].Nmembers,pf);

      if(Halos[i].NDomain_particles > 0)
	fwrite(&(Halos[i].Domain_particles[0]),sizeof(int),Halos[i].NDomain_particles,pf);
    }
  
  fclose(pf);
  
  return 0;
}


int main(int argc, char *argv[])
{
  
  int i, Tini, counter, iadvance, NGROUPS=0, Number_Processes;
  int info[3], ihalo, istart, iend, Nparts_per_node, ipart, Nparts_in_node;
  char *param_file=NULL;
  char buf[1000], bufhalo[1000];
  FILE *aux_pf=NULL;
  FILE *pf=NULL;
  
  int COLLECT_TAG=1,Nparticles_in_Node,k,j,*Info_domains=NULL,sendbuff;
  int NDomain_parts_per_proc,Parts_in_Domain_per_node,l,m,Old_NDomain_particles;
  MPI_Status status;
  
  struct halo *auxHalos=NULL;

  // Reading input parameters (TNG)
  char infile[MAX_FILENAME], cat_file[MAX_FILENAME];
  
  
  // Initializing MPI
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD,&task);
  MPI_Comm_size(MPI_COMM_WORLD,&Number_Processes);
  
  if(task == 0)
    printf(" * Running on %d processors\n",Number_Processes); fflush(stdout);
  
  if(argc != 4)
    usage_main();

  ///// TNG
  data_prefix    = argv[1];    // Ej: "data_" (TNG)
  catalog_prefix = argv[2];    // Ej: "catalog_" (TNG)
  
  //¿¿¿¿ESTO SOLO LO DEBERIA HACER EL TASK 0?????????
  snprintf(infile,   MAX_FILENAME, "%s.%d.hdf5", data_prefix,    0); // Gadget snapshot file
  snprintf(cat_file, MAX_FILENAME, "%s.%d.hdf5", catalog_prefix, 0); // FOF catalog file
  

  param_file = argv[3];
  Tini = tiempoc();
  
  GalCos_units(param_file);          // Loads the units and parameters from the parameter file
  GalCos_preload_Gadget(infile);     // Preload the Gadget snapshot file to get header information
  GalCos_preload_Catalog(cat_file);  // Preload the FOF catalog file to get the halo information
  
  NGROUPS = NCLUSTERS; // Total number of clusters (halos)
  
  /////////////////////////////////////////////////////////////////////////////
  /*                          DOMAIN DECOMPOSITION                           */
  /////////////////////////////////////////////////////////////////////////////
  
  /* 
     The unclustered particles are distributed across the processes.
     We dont have to read the info of particles inFOF halos. We only
     read the unclustered particles, and compute the domains for those. 
  */
  
  domain_info = (struct domain *) malloc((size_t) Number_Processes*sizeof(struct domain));
  
  if(task == 0) // If task 0, initialize the domain information
    {
      
      /* Distributing unclustered particles */
      Nparts_in_node = (int) (UNCLUSTERED/Number_Processes);
      
      printf(" Le toca %d particulas a cada nodo\n",Nparts_in_node); fflush(stdout);
      printf(" sobran %d particles\n",(UNCLUSTERED % Number_Processes)); fflush(stdout);

      istart = Npart_clustered; // Start from the first unclustered particle
      for(i=0; i<Number_Processes; i++)
	{
	    
	  domain_info[i].istart = istart;
	  domain_info[i].iend   = istart + Nparts_in_node;
	  domain_info[i].Nparts_per_node = domain_info[i].iend - domain_info[i].istart;
	  istart = domain_info[i].iend  ;
	}
      
      i = Number_Processes-1;
      domain_info[i].iend = UNCLUSTERED+Npart_clustered; // Last process takes the remaining particles
      domain_info[i].Nparts_per_node = domain_info[i].iend - domain_info[i].istart;
      
      for(i=0; i<Number_Processes; i++)
	printf("task %d start at %d and ends at %d, nparts=%d\n",i,domain_info[i].istart,domain_info[i].iend,
	       domain_info[i].Nparts_per_node); fflush(stdout);
      
    }
  
  MPI_Bcast(&domain_info[0], Number_Processes*sizeof(struct domain), MPI_BYTE, 0, MPI_COMM_WORLD);
  
  // GalCos_load_Gadget(&UNCLUSTERED, &NGROUPS, infile);
  MPI_Bcast(&NGROUPS,     1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&UNCLUSTERED, 1, MPI_INT, 0, MPI_COMM_WORLD);

  
  MPI_Barrier(MPI_COMM_WORLD);
  
  /////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////
  
  if(task == 0) 
    printf(" >> LOADING PARTICLE DATA FROM GADGET FILE...\n"); fflush(stdout);
  
  
  if(task == 0)	
    {
      printf(" %g Memory in particles = %g %d\n", sizeof(struct part)/(1024*1024.0),
	     UNCLUSTERED*sizeof(struct part)/(1024*1024.0), Npart_Total); fflush(stdout);
      
      printf(" %g Memory in halos = %g\n", sizeof(struct halo)/(1024*1024.0),
	     NGROUPS*sizeof(struct halo)/(1024*1024.0)); fflush(stdout);
      
      printf(" THERE ARE %d UNCLUSTERED PARTICLES (%f percent)\n",UNCLUSTERED,
	     100.0*(1.0*UNCLUSTERED)/(Npart_Total*1.0)); fflush(stdout);
      
    }
  

  /////////////////////////////////////////////////////////////////////////////////////
  //                 FILLING HALO STRUCTURES WITH PROPERTIES                         //
  /////////////////////////////////////////////////////////////////////////////////////

    if (Particle != NULL)
    free(Particle);
    
  if (domain_info[task].Nparts_per_node > 0) {
    
    Particle = (struct part *) malloc((size_t) domain_info[task].Nparts_per_node * sizeof(struct part));
    if (Particle == NULL) {
      printf("No memory available to load particles\n");
      exit(0);
    }
    
    }

    if(Halos != NULL)
      free(Halos);
  
  Halos = (struct halo *) malloc((size_t) NGROUPS*sizeof(struct halo));
  if(Halos == NULL)
    {
      printf("There are no memory left to allocate Halos\n");
      MPI_Finalize();
      exit(0);
    }

    printf("Ngroups = %d\n", NGROUPS); fflush(stdout);
    Galcos_loadTNG_HaloCatalog();

  NCLUSTERS = NGROUPS;
  MPI_Bcast(&NCLUSTERS, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&Halos[0], NCLUSTERS*sizeof(struct halo), MPI_BYTE, 0, MPI_COMM_WORLD);
  //// Reading particle positions
  GalCos_PartitionParticles(domain_info[task].istart, domain_info[task].iend);
  
  MPI_Barrier(MPI_COMM_WORLD);
  
  /////////////////////////////////////////////////////////////////////////////////
  /*                       Evaluating the domain of every halo                   */
  /////////////////////////////////////////////////////////////////////////////////
  
  // PASS 1: Count particles per halo
  counter=0;
  printf(" Searching in %d halos (Pass 1: Counting)\n", NGROUPS);
  for(i=0; i<domain_info[task].Nparts_per_node; i++)
    {
      if((counter%500000) == 0)
	printf(" (%d) *Computed domains for %d/%d at %f min\n",task,counter,domain_info[task].Nparts_per_node,
	       (tiempoc()-Tini)/60.0); fflush(stdout);
      
      GalCos_domain(i); // Mode 0: Count only
      counter++;
    }

  // ALLOCATE MEMORY
  for(i=0; i<NCLUSTERS; i++)
    {
      if(Halos[i].NDomain_particles > 0)
        {
          Halos[i].Domain_particles = (int *) malloc((size_t) Halos[i].NDomain_particles * sizeof(int));
          if(Halos[i].Domain_particles == NULL) {
             printf("Error allocating memory for halo %d\n", i);
             exit(0);
          }
          Halos[i].NDomain_particles = 0; // Reset counter for Pass 2
        }
      else
        {
          Halos[i].Domain_particles = NULL;
        }
    }

  // PASS 2: Fill arrays
  counter=0;
  printf(" Searching in %d halos (Pass 2: Filling)\n", NGROUPS);
  for(i=0; i<domain_info[task].Nparts_per_node; i++)
    {
      if((counter%500000) == 0)
	printf(" (%d) *Filling domains for %d/%d at %f min\n",task,counter,domain_info[task].Nparts_per_node,
	       (tiempoc()-Tini)/60.0); fflush(stdout);
      
      // Direct assignment using pre-calculated ID (No re-calculation!)
      int cid = Particle[i].Cluster_ID;
      if(cid != EMPTY_FLAG) {
          Halos[cid].Domain_particles[Halos[cid].NDomain_particles] = Particle[i].Oid;
          Halos[cid].NDomain_particles++;
      }
      
      counter++;
    }
  
  MPI_Barrier(MPI_COMM_WORLD);

  /////////////////////////////////////////////////////////////////////////////////
  /*      Sending Particle information from the other machines to the 0 task     */
  /////////////////////////////////////////////////////////////////////////////////
  
  /*
    sprintf(buf,"%s%s%d",infile,".parts.rescue.",task);
    pf=fopen(buf,"w");
    for(i=0; i<Nparts_per_node; i++)
    fprintf(pf,"%d %g %g %g\n",Particle[i].Cluster_ID,Particle[i].pos[0],Particle[i].pos[1],
    Particle[i].pos[2]);
    fclose(pf);
  */
  
  for(l=0; l<Number_Processes; l++)
    {
      
      if(task == l)
	{
	  sprintf(buf,"%s%s",infile,".parts.rescue");
	  
	  if(task == 0)
	    pf=fopen(buf,"w");
	  else
	    pf=fopen(buf,"a");
	  
	  for(i=0; i<domain_info[task].Nparts_per_node; i++)
	      fprintf(pf,"%d\n",Particle[i].Cluster_ID);
	  
	  fclose(pf);
	}
      
      MPI_Barrier(MPI_COMM_WORLD);
      
    }
  
  MPI_Barrier(MPI_COMM_WORLD);
  
  
  for(i=1; i<Number_Processes; i++)
    {
      
      if(task == i)
	{
	  
	  /* Sending information of Halos */
	  
	  Parts_in_Domain_per_node=0;
	  for(j=0; j<NCLUSTERS; j++)
	    Parts_in_Domain_per_node += Halos[j].NDomain_particles; 
	  
	  Info_domains = (int *) malloc((size_t) (Parts_in_Domain_per_node + 2*NCLUSTERS)*sizeof(int));
	  if(Info_domains == NULL)
	    {
	      printf("No memory available for allocation\n");
	      exit(0);
	    }
	  
	  k=0;
	  for(j=0; j<NCLUSTERS; j++)
	    {
	      
	      Info_domains[k] = Halos[j].IDcluster;
	      k++;
	      Info_domains[k] = Halos[j].NDomain_particles;
	      k++;
	      
	      for(l=0; l<Halos[j].NDomain_particles; l++)
		{
		  Info_domains[k] = Halos[j].Domain_particles[l];
		  k++;
		}
	      
	    }
	  
	  sendbuff = Parts_in_Domain_per_node + 2*NCLUSTERS;
	  MPI_Ssend(&sendbuff,1,MPI_INT,0,COLLECT_TAG,MPI_COMM_WORLD);
	  MPI_Ssend(&Info_domains[0],sendbuff,MPI_INT,0,COLLECT_TAG,MPI_COMM_WORLD);
	  
	  free(Info_domains);
	  
	}
      
    }
  
  /* Now I have to collect the data of particles from the other tasks */
  
  if(task == 0)
    {
      
      for(i=1; i<Number_Processes; i++)
	{
	  
	  /* Receiving data of halos */
	  
	  MPI_Recv(&sendbuff,1,MPI_INT,i,COLLECT_TAG,MPI_COMM_WORLD,&status);
	  
	  printf("Receiving info of halos from %d %d\n",i,sendbuff); fflush(stdout);
	  
	  Info_domains = (int *) malloc((size_t) sendbuff*sizeof(int));
	  if(Info_domains == NULL)
	    {
	      printf("No memory available for allocation\n");
	      exit(0);
	    }
	  
	  MPI_Recv(&Info_domains[0],sendbuff,MPI_INT,i,COLLECT_TAG,MPI_COMM_WORLD,&status);
	  printf("recibidos\n"); fflush(stdout);
	  
	  k=0;
	  for(j=0; j<NCLUSTERS; j++)
	    {
	      ihalo = Info_domains[k];
	      k++;
	      
	      Old_NDomain_particles = Halos[ihalo].NDomain_particles;
	      NDomain_parts_per_proc = Info_domains[k];
	      
	      Halos[ihalo].NDomain_particles += Info_domains[k];
	      k++;
	      
	      Halos[ihalo].Domain_particles = realloc(Halos[ihalo].Domain_particles,(size_t) Halos[ihalo].NDomain_particles*sizeof(int));
	      
	      m = Old_NDomain_particles;
	      for(l=0; l<NDomain_parts_per_proc; l++)
		{
		  Halos[ihalo].Domain_particles[m] = Info_domains[k];
		  k++;
		  m++;
		}
	    }
	  
	  free(Info_domains);
	  
	}
      
    }
  
  MPI_Barrier(MPI_COMM_WORLD);

  if(task == 0) 
    {
      printf("Writing info files\n");
      write_rescue(infile);
      
      printf("all is done!\n");
      fflush(stdout);


      pf = fopen("Halo_Catalog.dat", "w");
      for(j=0; j<NCLUSTERS; j++)
	{
	  fprintf(pf,"%16.8lf %16.8lf %16.8lf %16.8lf %16.8lf %16.8lf %12d\n", Halos[j].pos[0], Halos[j].pos[1], Halos[j].pos[2],
		  Halos[j].Mvir, Halos[j].Rvir, Halos[j].NDomain_particles*PARTMASS, Halos[j].NDomain_particles);
	}
      fclose(pf);
      
    }
  
  free(Particle);
  free(Halos);
  
  MPI_Finalize();
  
  return 0;
  
}
