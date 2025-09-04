#include "hdf5.h"
#include "GalCos_variables.h"

int Galcos_loadTNG_HaloCatalog(void)
{  
  
  int Halo_counter, Nmembers, fileid, NGROUPS_THISFILE, i;
  char cat_file[1000];
  hid_t file_id, group_id, dataset_id, dataspace_id, header_group, attr_id;
  herr_t status;
  hsize_t dims[1];
  int *group_members = NULL;
  float *group_mass = NULL, *group_mass_top_hat = NULL;
  float (*group_pos)[3] = NULL;
  float *group_radius = NULL;
  
  
  if(task == 0)
    {
      
      Halo_counter = 0; //Contador de halos
      
      for(fileid=0; fileid<NumFiles; fileid++)
        {
          
          printf("\nReading catalog %d.....\n", fileid); fflush(stdout);

          sprintf(cat_file, "%s.%d.hdf5", catalog_prefix, fileid);
          file_id = H5Fopen(cat_file, H5F_ACC_RDONLY, H5P_DEFAULT);
          
          if(file_id < 0) {
            printf("Cannot open %s\n", cat_file);
            MPI_Finalize();
            exit(0);
          }
	  
          //Read Header to get the number of groups in this file
	  
          if((header_group = H5Gopen(file_id, "/Header", H5P_DEFAULT)) < 0) {
            printf("Error opening Header group in %s\n", cat_file);
            exit(FAILURE);
          }
	  
          if((attr_id = H5Aopen(header_group, "Ngroups_ThisFile", H5P_DEFAULT)) < 0) {
            printf("Error opening NgroupsThisFile attribute in %s\n", cat_file);
            exit(FAILURE);
          }
	  
          status = H5Aread(attr_id, H5T_NATIVE_INT, &NGROUPS_THISFILE);
	  
          H5Aclose(attr_id);
          H5Gclose(header_group);
	  
          printf("Number of groups in this file: %d\n", NGROUPS_THISFILE); fflush(stdout);

          group_members = (int *) malloc(NGROUPS_THISFILE * sizeof(int));
          group_mass = (float *) malloc(NGROUPS_THISFILE * sizeof(float));
          group_pos = malloc(NGROUPS_THISFILE * sizeof(*group_pos));
          group_radius = (float *) malloc(NGROUPS_THISFILE * sizeof(float));
          group_mass_top_hat = (float *) malloc(NGROUPS_THISFILE * sizeof(float));
	  
	  
          /////////////////////////////////////////////////////////////////////////////////////
	  
          group_id = H5Gopen(file_id, "Group", H5P_DEFAULT);
          if(group_id < 0) {
            printf("Cannot open Group in %s\n", cat_file);
            MPI_Finalize();
            exit(0);
          }
	  
          //Open the GroupLen dataset to get the number of members in each group
          dataset_id = H5Dopen(group_id, "GroupLen", H5P_DEFAULT);
          if(dataset_id < 0) {
            printf("Cannot open GroupLen in %s\n", cat_file);
            MPI_Finalize();
            exit(0);
          }

          if(group_members == NULL) {
            printf("No memory available for group members\n");
            MPI_Finalize();
            exit(0);
          }
          status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, group_members);
          if(status < 0) {
            printf("Cannot read GroupLen in %s\n", cat_file);
            MPI_Finalize();
            exit(0);
          }
          H5Dclose(dataset_id);

          //Open the GroupMass dataset to get the mass of each group                                                                                                       
          dataset_id = H5Dopen(group_id, "GroupMass", H5P_DEFAULT);
          if(dataset_id < 0) {
            printf("Cannot open GroupMass in %s\n", cat_file);
            MPI_Finalize();
            exit(0);
          }


          if(group_mass == NULL) {
            printf("No memory available for group mass\n");
            MPI_Finalize();
            exit(0);
          }
          status = H5Dread(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, group_mass);
          if(status < 0) {
            printf("Cannot read GroupMass in %s\n", cat_file);
            MPI_Finalize();
            exit(0);
          }
          H5Dclose(dataset_id);

          //Open the GroupPos dataset to get the positions of each group                                                                                                   
          dataset_id = H5Dopen(group_id, "GroupPos", H5P_DEFAULT);
          if(dataset_id < 0) {
            printf("Cannot open GroupPos in %s\n", cat_file);
            MPI_Finalize();
            exit(0);
          }

          if(group_pos == NULL) {
            printf("No memory available for group positions\n");
            MPI_Finalize();
            exit(0);
          }
          status = H5Dread(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, group_pos);
          if(status < 0) {
            printf("Cannot read GroupPos in %s\n", cat_file);
            MPI_Finalize();
            exit(0);
          }
          H5Dclose(dataset_id);
	  
          //Open the Group_R_TopHat200 dataset to get the radius of each group                                                                                             
          dataset_id = H5Dopen(group_id, "Group_R_TopHat200", H5P_DEFAULT);
          if(dataset_id < 0) {
            printf("Cannot open Group_R_TopHat200 in %s\n", cat_file);
            MPI_Finalize();
            exit(0);
          }

          if(group_radius == NULL) {
            printf("No memory available for group radius\n");
            MPI_Finalize();
            exit(0);
          }
          status = H5Dread(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, group_radius);
          if(status < 0) {
            printf("Cannot read Group_R_TopHat200 in %s\n", cat_file);
            MPI_Finalize();
            exit(0);
          }
          H5Dclose(dataset_id);

          //Open the Group_M_TopHat200 dataset to get the mass of each group                                                                                               
          dataset_id = H5Dopen(group_id, "Group_M_TopHat200", H5P_DEFAULT);
          if(dataset_id < 0) {
            printf("Cannot open Group_M_TopHat200 in %s\n", cat_file);
            MPI_Finalize();
            exit(0);
          }

          if(group_mass_top_hat == NULL) {
            printf("No memory available for group mass top hat\n");
            MPI_Finalize();
            exit(0);
          }
          status = H5Dread(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, group_mass_top_hat);
          if(status < 0) {
            printf("Cannot read Group_M_TopHat200 in %s\n", cat_file);
            MPI_Finalize();
            exit(0);
          }
          H5Dclose(dataset_id);

          H5Gclose(group_id);
          H5Fclose(file_id);

	  
          for (i=0; i<NGROUPS_THISFILE; i++)
            {
              

	      
              Halos[Halo_counter].Nmembers = 0; //group_members[i];
              Halos[Halo_counter].mass     = group_mass[i];
              Halos[Halo_counter].NDomain_particles   = 0;
              Halos[Halo_counter].pos[0]   = group_pos[i][0];
              Halos[Halo_counter].pos[1]   = group_pos[i][1];
              Halos[Halo_counter].pos[2]   = group_pos[i][2];
              Halos[Halo_counter].Rvir     = group_radius[i]; // Assign the radius of the group                                                                            
              Halos[Halo_counter].Mvir     = group_mass_top_hat[i]; // Assign the mass of the group                                                                        
              Halos[Halo_counter].IDcluster = Halo_counter; // Assign the cluster ID                                                                                       
              Halos[Halo_counter].Domain_particles = NULL;
	      
              Halo_counter++;
	                                                                                                                                                        

              if (i==0)
                {
                  printf("Group %d has %d members\n", Halo_counter-1, group_members[i]); fflush(stdout);
                  printf("Group %d has mass %g\n", Halo_counter-1, group_mass[i]); fflush(stdout);
                  printf("Group %d position: (%g, %g, %g)\n", Halo_counter-1, group_pos[i][0], group_pos[i][1], group_pos[i][2]); fflush(stdout);
                  printf("Group %d radius (TopHat200): %g kpc\n", Halo_counter-1, group_radius[i]); fflush(stdout);
                  printf("Group %d mass (TopHat200): %ge10 Msun\n", Halo_counter-1, group_mass_top_hat[i]); fflush(stdout);
                }

            }

          
            free(group_members);
            group_members = NULL;
          
            free(group_mass); 
            group_mass = NULL;
          
            free(group_pos);
            group_pos = NULL;
          
            free(group_radius);
            group_radius = NULL;
          
            free(group_mass_top_hat);
            group_mass_top_hat = NULL;
          

        }// for field
      
      
    }// if task

  return 0;
}
