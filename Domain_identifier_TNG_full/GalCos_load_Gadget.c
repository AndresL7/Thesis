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
  
  // status = H5Aread(attr_id, H5T_NATIVE_INT, &Npart_clustered); // number of clustered particles
  H5Aclose(attr_id);

  // Calculate Npart_clustered by summing DM particles (Type 1) from GroupLenType across all catalog files
  // This is necessary for Full Physics runs where Nids_Total includes baryons
  
  int NumFilesCatalog;
  if((attr_id = H5Aopen(header_group, "NumFiles", H5P_DEFAULT)) < 0) {
      printf("Error opening NumFiles attribute in %s\n", cat_file);
      exit(FAILURE);
  }
  status = H5Aread(attr_id, H5T_NATIVE_INT, &NumFilesCatalog);
  H5Aclose(attr_id);

  Npart_clustered = 0;
  char current_cat_file[MAXPATHLEN];
  int ngroups_this;
  int f;
  
  // We close the initial file/group handles to avoid conflicts or just for cleanliness before the loop
  H5Gclose(header_group);
  H5Gclose(group_id);
  H5Fclose(file_id);

  for(f=0; f<NumFilesCatalog; f++) {
      sprintf(current_cat_file, "%s.%d.hdf5", catalog_prefix, f);
      hid_t fid = H5Fopen(current_cat_file, H5F_ACC_RDONLY, H5P_DEFAULT);
      if(fid < 0) {
          printf("Error opening catalog file %s\n", current_cat_file);
          continue;
      }
      
      hid_t hgrp = H5Gopen(fid, "/Header", H5P_DEFAULT);
      hid_t aid = H5Aopen(hgrp, "Ngroups_ThisFile", H5P_DEFAULT);
      status = H5Aread(aid, H5T_NATIVE_INT, &ngroups_this);
      H5Aclose(aid);
      H5Gclose(hgrp);
      
      if(ngroups_this > 0) {
          hid_t ggrp = H5Gopen(fid, "/Group", H5P_DEFAULT);
          hid_t did = H5Dopen(ggrp, "GroupLenType", H5P_DEFAULT);
          
          if(did < 0) {
             printf("Error opening GroupLenType in %s\n", current_cat_file);
             exit(FAILURE);
          }

          int *glt = (int*)malloc(ngroups_this * 6 * sizeof(int));
          status = H5Dread(did, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, glt);
          
          for(int k=0; k<ngroups_this; k++) {
              Npart_clustered += glt[k*6 + 1]; // Column 1 is DM
          }
          
          free(glt);
          H5Dclose(did);
          H5Gclose(ggrp);
      }
      H5Fclose(fid);
  }

  // Re-open file_id and header_group/group_id because the function expects to close them at the end?
  // The original code closes them at the end.
  // Let's just return 0 here since we are done.
  
  UNCLUSTERED = Npart_Total - Npart_clustered; // number of TOTAL unclustered particles
  
  if(task == 0) {
    // printf("Total number of clusters: %d\n", NCLUSTERS); fflush(stdout);
    // printf("Number of clustered particles: %d\n", Npart_clustered); fflush(stdout);
    printf("Number of TOTAL unclustered particles: %d\n", UNCLUSTERED); fflush(stdout);
    printf("====================================================\n");
  }

  return 0;

  /* Original cleanup code is unreachable now, which is fine as we closed everything */

  
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
