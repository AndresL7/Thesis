#include "GalCos_variables.h"
#include "hdf5.h"


#define MAXPATHLEN 256

struct gadget_head
{
  int      npart[6];
  double   mass[6];
  double   time;
  double   redshift;
  int      flag_sfr;
  int      flag_feedback;
  int      npartTotal[6];
  int      flag_cooling;
  int      num_files;
  double   BoxSize;
  double   Omega0;
  double   OmegaLambda;
  double   HubbleParam; 
  char     fill[256- 6*4- 6*8- 2*8- 2*4- 6*4- 2*4 - 4*8]; 
};

struct gadget_head header1;



//// For TNG
/////////////Preload Gadget new hdf5 ///////////////////////////////////////////////////////////

int GalCos_preload_Gadget(char *infile)
{

  hid_t file_id, header_group, attr_id;
  herr_t status;
  hsize_t dims[2];
  int npart[6];
  double mass_table[6];
  
  printf("%s",infile);
  if((file_id = H5Fopen(infile, H5F_ACC_RDONLY, H5P_DEFAULT)) < 0) {
    printf("Error opening file %s\n", infile);
    exit(FAILURE);
  }
  
  if((header_group = H5Gopen(file_id, "/Header", H5P_DEFAULT)) < 0) {
    printf("Error opening Header group in %s\n", infile);
    exit(FAILURE);
  }
    
  if((attr_id = H5Aopen(header_group, "NumPart_Total", H5P_DEFAULT)) < 0) {
    printf("Error opening NumPart_Total attribute in %s\n", infile);
    exit(FAILURE);
  }
  
  status = H5Aread(attr_id, H5T_STD_U32LE, &npart[0]); // number of particles in each category
  H5Aclose(attr_id);
  Npart_Total = npart[1]; //number of particles in the second category (DM particles)  
  if(task == 0) printf("Total number of particles: %d\n", Npart_Total); fflush(stdout);

  if((attr_id = H5Aopen(header_group, "MassTable", H5P_DEFAULT)) < 0) {
    printf("Error opening MassTable attribute in %s\n", infile);
    exit(FAILURE);
  }
  
  status = H5Aread(attr_id, H5T_IEEE_F64LE, &mass_table[0]); // mass of particles in each category
  H5Aclose(attr_id);

  PARTMASS = mass_table[1]; // mass of DM particles

  if((attr_id = H5Aopen(header_group, "Time", H5P_DEFAULT)) < 0) {
    printf("Error opening Time attribute in %s\n", infile);
    exit(FAILURE);
  }
  status = H5Aread(attr_id, H5T_IEEE_F64LE, &COSMIC_TIME); // cosmic time
  H5Aclose(attr_id);

  if((attr_id = H5Aopen(header_group, "Redshift", H5P_DEFAULT)) < 0) {
    printf("Error opening Redshift attribute in %s\n", infile);
    exit(FAILURE);
  }
  status = H5Aread(attr_id, H5T_IEEE_F64LE, &REDSHIFT); // redshift
  H5Aclose(attr_id);

  if((attr_id = H5Aopen(header_group, "Flag_Sfr", H5P_DEFAULT)) < 0) {
    printf("Error opening FlagSfr attribute in %s\n", infile);
    exit(FAILURE);
  }
  status = H5Aread(attr_id, H5T_STD_I32LE, &FlagSfr); // star formation flag
  H5Aclose(attr_id);

  if((attr_id = H5Aopen(header_group, "Flag_Feedback", H5P_DEFAULT)) < 0) {
    printf("Error opening FlagFeedback attribute in %s\n", infile);
    exit(FAILURE);
  }
  status = H5Aread(attr_id, H5T_STD_I32LE, &FlagFeedback); // feedback flag
  H5Aclose(attr_id);

  if((attr_id = H5Aopen(header_group, "Flag_Cooling", H5P_DEFAULT)) < 0) {
    printf("Error opening FlagCooling attribute in %s\n", infile);
    exit(FAILURE);
  }
  status = H5Aread(attr_id, H5T_STD_I32LE, &FlagCooling); // cooling flag
  H5Aclose(attr_id);

  if((attr_id = H5Aopen(header_group, "NumFilesPerSnapshot", H5P_DEFAULT)) < 0) {
    printf("Error opening NumFiles attribute in %s\n", infile);
    exit(FAILURE);
  }
  status = H5Aread(attr_id, H5T_STD_I32LE, &NumFiles); // number of files
  H5Aclose(attr_id);

  if((attr_id = H5Aopen(header_group, "BoxSize", H5P_DEFAULT)) < 0) {
    printf("Error opening BoxSize attribute in %s\n", infile);
    exit(FAILURE);
  }
  status = H5Aread(attr_id, H5T_IEEE_F64LE, &BoxSize); // box size
  H5Aclose(attr_id);

  if((attr_id = H5Aopen(header_group, "Omega0", H5P_DEFAULT)) < 0) {
    printf("Error opening Omega0 attribute in %s\n", infile);
    exit(FAILURE);
  }
  status = H5Aread(attr_id, H5T_IEEE_F64LE, &OMEGA_MATTER); // matter density parameter
  H5Aclose(attr_id);

  if((attr_id = H5Aopen(header_group, "OmegaLambda", H5P_DEFAULT)) < 0) {
    printf("Error opening OmegaLambda attribute in %s\n", infile);
    exit(FAILURE);
  }
  status = H5Aread(attr_id, H5T_IEEE_F64LE, &OMEGALAMBDA); // dark energy density parameter
  H5Aclose(attr_id);

  if((attr_id = H5Aopen(header_group, "HubbleParam", H5P_DEFAULT)) < 0) {
    printf("Error opening HubbleParam attribute in %s\n", infile);
    exit(FAILURE);
  }
  status = H5Aread(attr_id, H5T_IEEE_F64LE, &HUBBLEPARAM); // Hubble parameter
  H5Aclose(attr_id);

  if(task == 0)
    {
      printf("* Frame's Time... %g\n", COSMIC_TIME); fflush(stdout);
      printf("* Redshift... %g\n", REDSHIFT); fflush(stdout);
      printf("* FlagSfr... %d\n", FlagSfr); fflush(stdout);
      printf("* FlagFeedback... %d\n", FlagFeedback); fflush(stdout);
      printf("* FlagCooling... %d\n", FlagCooling); fflush(stdout);
      printf("* NumFiles... %d\n", NumFiles); fflush(stdout);
      printf("* BoxSize... %g\n", BoxSize); fflush(stdout);
      printf("* Omega0... %g\n", OMEGA_MATTER); fflush(stdout);
      printf("* OmegaLambda... %g\n", OMEGALAMBDA); fflush(stdout);
      printf("* HubbleParam... %g\n", HUBBLEPARAM); fflush(stdout);
      printf("* Particle mass... %g\n", PARTMASS); fflush(stdout);
      printf("====================================================\n");
    }
  H5Gclose(header_group);
  H5Fclose(file_id);

  return 0;
  
}

// Preload Catalog Gadget
int GalCos_preload_Catalog(char *cat_file)
{

  hid_t file_id, header_group, attr_id;
  herr_t status;
  hid_t group_id;

  char buff[100];
  int i, j, k, sorted_groups, auxint, counter, maxsend, maxrecv;
  struct part dp;
  FILE *aux_pf = NULL;
  FILE *fp_inp = NULL;

  // Open the catalog file
  if((file_id = H5Fopen(cat_file, H5F_ACC_RDONLY, H5P_DEFAULT)) < 0) {
    printf("Error opening catalog file %s\n", cat_file);
    exit(FAILURE);
  }

  // Open the group containing the halo data
  if((group_id = H5Gopen(file_id, "/Group", H5P_DEFAULT)) < 0) {
    printf("Error opening Group in %s\n", cat_file);
    exit(FAILURE);
  }

  //Header
  if((header_group = H5Gopen(file_id, "/Header", H5P_DEFAULT)) < 0) {
    printf("Error opening Header group in %s\n", cat_file);
    exit(FAILURE);
  }

  if((attr_id = H5Aopen(header_group, "Ngroups_Total", H5P_DEFAULT)) < 0) {
    printf("Error opening NgroupsTotal attribute in %s\n", cat_file);
    exit(FAILURE);
  }
  status = H5Aread(attr_id, H5T_NATIVE_INT, &NCLUSTERS); // total number of clusters
  H5Aclose(attr_id);

  //if((attr_id = H5Aopen(header_group, "Nids_ThisFile", H5P_DEFAULT)) < 0) {
  //printf("Error opening NidsThisFile attribute in %s\n", cat_file);
  //exit(FAILURE);
  //}

  if((attr_id = H5Aopen(header_group, "Nids_Total", H5P_DEFAULT)) < 0) {
    printf("Error opening Nids_Total attribute in %s\n", cat_file);
    exit(FAILURE);
  }
  
  status = H5Aread(attr_id, H5T_NATIVE_INT, &Npart_clustered); // number of clustered particles
  H5Aclose(attr_id);
  
  UNCLUSTERED = Npart_Total - Npart_clustered; // number of TOTAL unclustered particles
  
  if(task == 0) {
    // printf("Total number of clusters: %d\n", NCLUSTERS); fflush(stdout);
    // printf("Number of clustered particles: %d\n", Npart_clustered); fflush(stdout);
    printf("Number of TOTAL unclustered particles: %d\n", UNCLUSTERED); fflush(stdout);
    printf("====================================================\n");
  }

  H5Gclose(group_id);
  H5Fclose(file_id);

  return 0;
  
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Load Gadget New
/*
  int GalCos_load_Gadget(int *unclust, int *sacum, char *infile)
  {
  // New loading mechanism for Gadget files
  hid_t file_id, group_id, dataset_id, dataspace_id;
  herr_t status;
  hsize_t dims[2];
  int i, j, k, sorted_groups, auxint, counter, maxsend, maxrecv;
  struct part dp;
  char buff[100];
  long int bloksize;
  FILE *aux_pf = NULL;
  FILE *fp_inp = NULL;
  
  // Open the input HDF5 file
  // if((file_id = H5Fopen(infile, H5F_ACC_RDONLY, H5P_DEFAULT)) < 0) {
  //   printf("read_gadget cannot open %s\n", infile);
  //   exit(FAILURE);
  // }
  
  // Unclustered particles and clusters
  (*unclust) = Npart_unclustered; // set unclustered particles
  (*sacum) = NCLUSTERS; // set the number of clusters
  
  if (task == 0) {
  domain_info[task].istart = 0;
  domain_info[task].iend = Npart_Total;
  
  domain_info[task].Nparts_per_node = Npart_Total - (*unclust);
  printf("There are %d clustered particles\n", domain_info[task].Nparts_per_node); fflush(stdout);
  printf("There are %d FOF groups\n", (*sacum)); fflush(stdout);
  }
  
  
  // creo que ya no necesito esto
  if (Particle != NULL) 
  free(Particle);
  
  if (domain_info[task].Nparts_per_node > 0) {
  Particle = (struct part *) malloc((size_t) domain_info[task].Nparts_per_node * sizeof(struct part));
  if (Particle == NULL) {
  printf("No memory available to load particles\n");
  exit(0);
  }
  }
  
  MPI_Barrier(MPI_COMM_WORLD);
  
  return 0;  
  
  }
*/
