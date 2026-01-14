#include "GalCos_variables.h"

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
  char     fill[256- 6*4- 6*8- 2*8- 2*4- 6*4- 2*4 - 4*8];  /* fills to 
							      256 Bytes */
};

struct gadget_head header1;

int GalCos_preload_Gadget(char *infile)
{
  //Leer hdf5
  hid_t file_id, dataset_id, dataspace_id, header_group, attr_id;
  herr_t status;
  hsize_t dims[2];

  // Abre el archivo HDF5
  file_id = H5Fopen(infile, H5F_ACC_RDONLY, H5P_DEFAULT);
  header_group = H5Gopen(file_id, "Header", H5P_DEFAULT);

  //Redshift
  attr_id = H5Aopen(header_group, "Redshift", H5P_DEFAULT);
  status = H5Aread(attr_id, H5T_NATIVE_FLOAT, &REDSHIFT);
  H5Aclose(attr_id);
  //Omega0
  attr_id = H5Aopen(header_group, "Omega0", H5P_DEFAULT);
  status = H5Aread(attr_id, H5T_NATIVE_FLOAT, &OMEGA_MATTER);
  H5Aclose(attr_id);
  //OmegaLambda
  attr_id = H5Aopen(header_group, "OmegaLambda", H5P_DEFAULT);
  status = H5Aread(attr_id, H5T_NATIVE_FLOAT, &OMEGALAMBDA);
  H5Aclose(attr_id);
  //HubbleParam
  attr_id = H5Aopen(header_group, "HubbleParam", H5P_DEFAULT);
  status = H5Aread(attr_id, H5T_NATIVE_FLOAT, &HUBBLEPARAM);
  H5Aclose(attr_id);
  //Time???? (Cosmic Time)
  attr_id = H5Aopen(header_group, "Time", H5P_DEFAULT);
  status = H5Aread(attr_id, H5T_NATIVE_FLOAT, &COSMIC_TIME);
  H5Aclose(attr_id);
  //FlagSfr
  attr_id = H5Aopen(header_group, "Flag_Sfr", H5P_DEFAULT);
  status = H5Aread(attr_id, H5T_NATIVE_INT, &FlagSfr);
  H5Aclose(attr_id);
  //FlagFeedback
  attr_id = H5Aopen(header_group, "Flag_Feedback", H5P_DEFAULT);
  status = H5Aread(attr_id, H5T_NATIVE_INT, &FlagFeedback);
  H5Aclose(attr_id);
  //FlagCooling
  attr_id = H5Aopen(header_group, "Flag_Cooling", H5P_DEFAULT);
  status = H5Aread(attr_id, H5T_NATIVE_INT, &FlagCooling);
  H5Aclose(attr_id);
  //NumFiles
  attr_id = H5Aopen(header_group, "NumFilesPerSnapshot", H5P_DEFAULT);
  status = H5Aread(attr_id, H5T_NATIVE_INT, &NumFiles);
  H5Aclose(attr_id);
  //BoxSize
  attr_id = H5Aopen(header_group, "BoxSize", H5P_DEFAULT);
  status = H5Aread(attr_id, H5T_NATIVE_FLOAT, &BoxSize);
  H5Aclose(attr_id);
  //NpartTotal
  int particles[6];
  // attr_id = H5Aopen(header_group, "NumPart_ThisFile", H5P_DEFAULT);
  attr_id = H5Aopen(header_group, "NumPart_Total", H5P_DEFAULT);
  status = H5Aread(attr_id, H5T_NATIVE_INT, &particles[0]);
  H5Aclose(attr_id);
  Npart_Total=particles[1];
  printf("Npart_Total: %d\n", Npart_Total);
  
  printf("\n");
  printf(" * Frame's Time... %g\n",COSMIC_TIME); 
  printf(" * Redshift... %g\n",REDSHIFT);
  printf(" * Flagsfr... %d\n",FlagSfr);
  printf(" * Flagfed... %d\n",FlagFeedback);
  printf(" * Flagcool... %d\n",FlagCooling);
  printf(" * Numfiles... %d\n",NumFiles);
  printf(" * Boxsize... %g\n",BoxSize);
  printf(" * Omega0... %g\n",OMEGA_MATTER);
  printf(" * OmageLambda... %g\n",OMEGALAMBDA);
  printf(" * Hubbleparam... %g\n",HUBBLEPARAM);
      
  H5Gclose(header_group);
  H5Fclose(file_id);

  return 0;
  
}


/// @brief 
/// @param infile 
/// @param overlap 
/// @param l 
/// @return 
int GalCos_load_Gadget(char *infile,float overlap,int l)
{
  
  int i,NPH,j,oNPH,task;
  float Plenght,Npos[3],posa[3];
  struct part dp;
  hid_t file_id, dataset_id, dataspace_id;
  herr_t status;
  hsize_t dims[2];

  task = l;
  ///////////////////////////////////Leer HDF5///////////////////

  printf("Number of particles: %d\n", Npart_Total);
  
  NPH = (int) (2*(Npart_Total/Number_Processes));
  
  if(task == 0)
      printf("Initially %d particles per node\n",NPH); fflush(stdout);
  
  if(NPH != 0)
    {
      Particle = (struct Vpart *) malloc((size_t) NPH*sizeof(struct Vpart)); 
      if(Particle == NULL){
	printf("No memory available to load particles \n");
	exit(0);
      }
    }
    
  //////////////////// positions
    
  domain_info[task].Nparticles = 0;
  Plenght = pow(0.5*BoxSize,2);
  printf("\nPlenght: %f\n", Plenght);
   
  float coord[3];

  posa[0] = domain_info[task].pos[0];
  posa[1] = domain_info[task].pos[1];
  posa[2] = domain_info[task].pos[2];

for(i=0; i<Npart_Total; i++) 
    {
      dp.pos[0] = Rpart[i].pos[0];
      dp.pos[1] = Rpart[i].pos[1];
      dp.pos[2] = Rpart[i].pos[2];

      Npos[0] = dp.pos[0];
      Npos[1] = dp.pos[1];
      Npos[2] = dp.pos[2];
      

      if(pow(dp.pos[0]-posa[0],2) > Plenght)
	{
	  if(posa[0] > dp.pos[0])
	    Npos[0] = dp.pos[0] + BoxSize;
	  else
	    Npos[0] = dp.pos[0] - BoxSize;
	}
      
      if(pow(dp.pos[1]-posa[1],2) > Plenght)
	{
	  if(posa[1] > dp.pos[1])
		Npos[1] = dp.pos[1] + BoxSize;
	  else
	    Npos[1] = dp.pos[1] - BoxSize;
	}
      
      if(pow(dp.pos[2]-posa[2],2) > Plenght)
	{
	  if(posa[2] > dp.pos[2])
	    Npos[2] = dp.pos[2] + BoxSize;
	  else
	    Npos[2] = dp.pos[2] - BoxSize;
	}
    if((Npos[0] >= domain_info[task].xmin) && (Npos[0] <= domain_info[task].xmax)){
     
	    if((Npos[1] >= domain_info[task].ymin) && (Npos[1] <= domain_info[task].ymax)){
       
	      if((Npos[2] >= domain_info[task].zmin) && (Npos[2] <= domain_info[task].zmax))
	    {
      
	      for(j=0; j<3; j++)
		    Particle[domain_info[task].Nparticles].pos[j] = Npos[j];
	      Particle[domain_info[task].Nparticles].Oid = i;
	      Particle[domain_info[task].Nparticles].tag = 0;
	      	      
	      /////////******************///////
	      
	      if((dp.pos[0] >= domain_info[task].xmind) && (dp.pos[0] <= domain_info[task].xmaxd)){
		  
		  if((dp.pos[1] >= domain_info[task].ymind) && (dp.pos[1] <= domain_info[task].ymaxd)){
		      
		      if((dp.pos[2] >= domain_info[task].zmind) && (dp.pos[2] <= domain_info[task].zmaxd))
		      {
			  Particle[domain_info[task].Nparticles].tag = 1;
		      }
		      
		  }
		  
	      }
	      
	      if(domain_info[task].Nparticles >= (NPH-1))
	      {
		  oNPH=NPH;
		  NPH = NPH + (int) (0.1*NPH);
		  printf("%d Increasing the value for particle allocation from %d to %d\n",task,oNPH,NPH);
		  Particle = (struct Vpart *) realloc(Particle,(size_t) NPH*sizeof(struct Vpart));
		  
		  if(Particle == NULL)
		  {
		      printf("%d problem of reallocation\n",task);
		      exit(0);
		  }
	      }
	      
	      domain_info[task].Nparticles++;
	      
	    }
	}
      }
    }
  
  printf("Task %d has %d particles\n",task,domain_info[task].Nparticles);
  
  Particle = (struct Vpart *) realloc(Particle,(size_t) domain_info[task].Nparticles*sizeof(struct Vpart));
  
  return domain_info[task].Nparticles;
 
  
}
