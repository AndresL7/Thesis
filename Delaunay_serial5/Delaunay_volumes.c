#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include <gsl/gsl_linalg.h>
#include <libqhull_r/qhull_ra.h>

int mybsearch(int *array,int N,int elemento)
{
  
  int desde,hasta,medio,posicion=-1;
  
  for(desde=0,hasta=N-1;desde<=hasta;)
    {
      
      if(desde==hasta)
        {
          if(array[desde]==elemento)
            posicion=desde;
          else
            posicion=-1;
	  
          break;
        }
      
      medio=(desde+hasta)/2;
      
      if(array[medio]==elemento)
        {
          posicion=medio;
          break;
	}
      else if(array[medio]>elemento)
	hasta=medio-1;
      else
        desde=medio+1;
    }
  
  return posicion;
}


// int fill_matrix(double *distances,int **facets,int ifacet)
// {
//   ////////////////////////////
  
//   distances[0] = 0.0;
//   distances[1] = 1.0;    
//   distances[2] = 1.0;    
//   distances[3] = 1.0;    
//   distances[4] = 1.0;    
  
//   ///////////////////////////
  
//   distances[5] = 1.0;    
//   distances[6] = 0.0;
//   //1,2//
//   //Particle[i].pos[0]
//   //coords[i][0]
//   distances[7] = distance(Particle[facets[ifacet][0]].pos[0],Particle[facets[ifacet][0]].pos[1],Particle[facets[ifacet][0]].pos[2],
//   Particle[facets[ifacet][1]].pos[0],Particle[facets[ifacet][1]].pos[1],Particle[facets[ifacet][1]].pos[2]);
//   distances[7] = distances[7]*distances[7];
  
//   //1,3//
//   distances[8] = distance(Particle[facets[ifacet][0]].pos[0],Particle[facets[ifacet][0]].pos[1],Particle[facets[ifacet][0]].pos[2],
//   Particle[facets[ifacet][2]].pos[0],Particle[facets[ifacet][2]].pos[1],Particle[facets[ifacet][2]].pos[2]);
//   distances[8] = distances[8]*distances[8];
  
//   //1,4//
//   distances[9] = distance(Particle[facets[ifacet][0]].pos[0],Particle[facets[ifacet][0]].pos[1],Particle[facets[ifacet][0]].pos[2],
//   Particle[facets[ifacet][3]].pos[0],Particle[facets[ifacet][3]].pos[1],Particle[facets[ifacet][3]].pos[2]);
//   distances[9] = distances[9]*distances[9];
  
//   //////////////////////////
  
//   distances[10] = 1.0;
//   //2,1//
//   distances[11] = distances[7];
    
//   distances[12] = 0.0;
  
//   //2,3//
//   distances[13] = distance(Particle[facets[ifacet][1]].pos[0],Particle[facets[ifacet][1]].pos[1],Particle[facets[ifacet][1]].pos[2],
// 			   Particle[facets[ifacet][2]].pos[0],Particle[facets[ifacet][2]].pos[1],Particle[facets[ifacet][2]].pos[2]);
//   distances[13] = distances[13]*distances[13];
  
//   //2,4//
//   distances[14] = distance(Particle[facets[ifacet][1]].pos[0],Particle[facets[ifacet][1]].pos[1],Particle[facets[ifacet][1]].pos[2],
// 			   Particle[facets[ifacet][3]].pos[0],Particle[facets[ifacet][3]].pos[1],Particle[facets[ifacet][3]].pos[2]);
//   distances[14] = distances[14]*distances[14];
  
//   ///////////////////////////
    
//   distances[15] = 1.0;
//   distances[16] = distances[8];
//   distances[17] = distances[13]; 
//   distances[18] = 0.0;
  
//   //3,4//
//   distances[19] = distance(Particle[facets[ifacet][2]].pos[0],Particle[facets[ifacet][2]].pos[1],Particle[facets[ifacet][2]].pos[2],
// 			   Particle[facets[ifacet][3]].pos[0],Particle[facets[ifacet][3]].pos[1],Particle[facets[ifacet][3]].pos[2]);
//   distances[19] = distances[19]*distances[19];
  
//   ///////////////////////////    
  
//   distances[20] = 1.0;
//   distances[21] = distances[9]; 
//   distances[22] = distances[14]; 
//   distances[23] = distances[19]; 
//   distances[24] = 0.0;
  
//   return 0;
  
// }
int fill_matrix(double *distances,int **facets,double **coords,int ifacet)
{
  ////////////////////////////
  
  distances[0] = 0.0;
  distances[1] = 1.0;    
  distances[2] = 1.0;    
  distances[3] = 1.0;    
  distances[4] = 1.0;    
  
  ///////////////////////////
  
  distances[5] = 1.0;    
  distances[6] = 0.0;
  //1,2//
  distances[7] = distance(coords[facets[ifacet][0]][0],coords[facets[ifacet][0]][1],coords[facets[ifacet][0]][2],
  coords[facets[ifacet][1]][0],coords[facets[ifacet][1]][1],coords[facets[ifacet][1]][2]);
  distances[7] = distances[7]*distances[7];
  
  //1,3//
  distances[8] = distance(coords[facets[ifacet][0]][0],coords[facets[ifacet][0]][1],coords[facets[ifacet][0]][2],
  coords[facets[ifacet][2]][0],coords[facets[ifacet][2]][1],coords[facets[ifacet][2]][2]);
  distances[8] = distances[8]*distances[8];
  
  //1,4//
  distances[9] = distance(coords[facets[ifacet][0]][0],coords[facets[ifacet][0]][1],coords[facets[ifacet][0]][2],
  coords[facets[ifacet][3]][0],coords[facets[ifacet][3]][1],coords[facets[ifacet][3]][2]);
  distances[9] = distances[9]*distances[9];
  
  //////////////////////////
  
  distances[10] = 1.0;
  //2,1//
  distances[11] = distances[7];
    
  distances[12] = 0.0;
  
  //2,3//
  distances[13] = distance(coords[facets[ifacet][1]][0],coords[facets[ifacet][1]][1],coords[facets[ifacet][1]][2],
			   coords[facets[ifacet][2]][0],coords[facets[ifacet][2]][1],coords[facets[ifacet][2]][2]);
  distances[13] = distances[13]*distances[13];
  
  //2,4//
  distances[14] = distance(coords[facets[ifacet][1]][0],coords[facets[ifacet][1]][1],coords[facets[ifacet][1]][2],
			   coords[facets[ifacet][3]][0],coords[facets[ifacet][3]][1],coords[facets[ifacet][3]][2]);
  distances[14] = distances[14]*distances[14];
  
  ///////////////////////////
    
  distances[15] = 1.0;
  distances[16] = distances[8];
  distances[17] = distances[13]; 
  distances[18] = 0.0;
  
  //3,4//
  distances[19] = distance(coords[facets[ifacet][2]][0],coords[facets[ifacet][2]][1],coords[facets[ifacet][2]][2],
			   coords[facets[ifacet][3]][0],coords[facets[ifacet][3]][1],coords[facets[ifacet][3]][2]);
  distances[19] = distances[19]*distances[19];
  
  ///////////////////////////    
  
  distances[20] = 1.0;
  distances[21] = distances[9]; 
  distances[22] = distances[14]; 
  distances[23] = distances[19]; 
  distances[24] = 0.0;
  
  return 0;
  
}

double compute_det(double *data)
{
  int s;
  double Det;
  
  gsl_matrix_view m = gsl_matrix_view_array(data,5,5);
  gsl_permutation *p = gsl_permutation_alloc(5);
  gsl_linalg_LU_decomp(&m.matrix, p, &s);
  
  Det=gsl_linalg_LU_det(&m.matrix,s);
  gsl_permutation_free(p);
    
  return Det;
  
}


int print_delaunuay(int **facets, int i, int j)
{
  
  FILE *pf;
  
  pf=fopen("plotador.gpl","w");
  
  fprintf(pf,"set arrow from %f,%f,%f to %f,%f,%f nohead\n",Particle[facets[i][0]].pos[0],Particle[facets[i][0]].pos[1],
	  Particle[facets[i][0]].pos[2],Particle[facets[i][1]].pos[0],Particle[facets[i][1]].pos[1],Particle[facets[i][1]].pos[2]);
  
  fprintf(pf,"set arrow from %f,%f,%f to %f,%f,%f nohead\n",Particle[facets[i][0]].pos[0],Particle[facets[i][0]].pos[1],
	  Particle[facets[i][0]].pos[2],Particle[facets[i][2]].pos[0],Particle[facets[i][2]].pos[1],Particle[facets[i][2]].pos[2]);
  
  fprintf(pf,"set arrow from %f,%f,%f to %f,%f,%f nohead\n",Particle[facets[i][0]].pos[0],Particle[facets[i][0]].pos[1],
	  Particle[facets[i][0]].pos[2],Particle[facets[i][3]].pos[0],Particle[facets[i][3]].pos[1],Particle[facets[i][3]].pos[2]);
    
  fprintf(pf,"set arrow from %f,%f,%f to %f,%f,%f nohead\n",Particle[facets[i][1]].pos[0],Particle[facets[i][1]].pos[1],
	  Particle[facets[i][1]].pos[2],Particle[facets[i][3]].pos[0],Particle[facets[i][3]].pos[1],Particle[facets[i][3]].pos[2]);
  
  fprintf(pf,"set arrow from %f,%f,%f to %f,%f,%f nohead\n",Particle[facets[i][1]].pos[0],Particle[facets[i][1]].pos[1],
	  Particle[facets[i][1]].pos[2],Particle[facets[i][2]].pos[0],Particle[facets[i][2]].pos[1],Particle[facets[i][2]].pos[2]);
  
  fprintf(pf,"set arrow from %f,%f,%f to %f,%f,%f nohead\n",Particle[facets[i][2]].pos[0],Particle[facets[i][2]].pos[1],
	  Particle[facets[i][2]].pos[2],Particle[facets[i][3]].pos[0],Particle[facets[i][3]].pos[1],Particle[facets[i][3]].pos[2]);
  
  fprintf(pf,"set label 'Here' at %f,%f,%f\n",Particle[facets[i][j]].pos[0],Particle[facets[i][j]].pos[1],Particle[facets[i][j]].pos[2]);
  
  fclose(pf);
  
  return 0;
}


int Delaunay_volumes(int Htask)
{
    
  int i,j,**facets=NULL,k,SIZEMATRIX,NFACETS,counter,dim=3;
  double **coords=NULL,*distances=NULL,*Volumenes=NULL,Vi;
  
  task = Htask;
  printf("V2");
  printf("entering to delaunay volumes in task %d\n",task); fflush(stdout);

  //printf("domain_info[task].Nparticles = %d\n",domain_info[task].Nparticles);
  
  /*Points*/
  int n=27;
  int num_points=n*n*n;
  printf("num_points = %d\n",num_points);
  
  //Particle = (struct part *) malloc((size_t) *sizeof(struct part)); 

  // Allocating space for coords
  coords=(double **) malloc((size_t) num_points*sizeof(double *));
  if(coords == NULL){
      printf("No memory available to allocate coords\n");
      exit(0);
  }
  
  for(i=0; i<num_points; i++)
  {
      coords[i]=(double *) malloc((size_t) 3*sizeof(double));
      if(coords[i] == NULL){
	  printf("No memory available to allocate coords2\n");
	  exit(0);
      }
  }

  //Delaunay magic

    coordT *points = malloc(num_points * dim * sizeof(coordT));
    
    if (points == NULL) {
        fprintf(stderr, "Error allocating memory for points.\n");
        exit(0);
    }

    int K = 0;
    int ID = 0;
    double dx = 100.0/(n-1);
    double x,y,z;
    for(int i = 0;i < n; i++){
      x = i * dx;
      for(int j = 0;j < n; j++){
        y = j * dx;
        for(int k = 0;k < n; k++){
          z = k * dx;
          // Particle[ID].pos[0] = x;
          // Particle[ID].pos[1] = y;
          // Particle[ID].pos[2] = z;
          coords[ID][0] = x;
          coords[ID][1] = y;
          coords[ID][2] = z;
          ID++;
          points[K++] = x;
          points[K++] = y;
          points[K++] = z;
        }
      }
    }
    printf("Last particle: %d %f %f %f\n",ID-1,coords[ID-1][0],coords[ID-1][1],coords[ID-1][2]);

    qhT qh_qh;
    qhT *qh = &qh_qh;
    QHULL_LIB_CHECK
    qh_zero(qh, stderr);
    const char *options = "qhull d QJ"; // 'd' for Delaunay, 'QJ' for exactitude 
    int exitcode = qh_new_qhull(qh, dim, num_points, points, False, (char*)options, NULL, stderr);

    if (!exitcode) {
        facetT *facet;
        vertexT *vertex, **vertexp;
        
        //qh_getarea(qh, qh->facet_list);
        NFACETS = qh->num_good;

    //Total Volume
        facets = (int **) malloc((size_t) NFACETS*sizeof(int *));
        int facetid=0;
        FORALLfacets {
          if(!facet->upperdelaunay) {
            facets[facetid] = (int *) malloc((size_t) 4*sizeof(int));
            
            int i=0;
            FOREACHvertex_(facet->vertices) {
                pointT *point = vertex->point;
    //indices each vertex
                int idv = qh_pointid(qh, point);
    //Here are stored the indices of the vertices of each facet
                facets[facetid][i] = idv;
                i++;
            }
            facetid++;
        }
        }   
        
    } else {
        fprintf(stderr, "Error computing convex hull.\n");
    }

    qh_freeqhull(qh, !qh_ALL);
    int curlong, totlong;
    qh_memfreeshort(qh, &curlong, &totlong);
    if (curlong || totlong) {
        fprintf(stderr, "Qhull internal warning: did not free %d bytes of long memory (%d pieces)\n", totlong, curlong);
    }
  free(points);
  //printf("Delaunay Computed\n");

  // Allocating space for volumes
  Volumenes = (double *) malloc((size_t) NFACETS*sizeof(double));
  if(Volumenes == NULL){
      printf("No memory available to allocate Volumenes\n");
      exit(0);
  }
  //printf("\nVolumes Alocated\n");
  
  // Allocating space for matrix distances
  SIZEMATRIX = 5;
  distances = (double *) malloc((size_t) SIZEMATRIX*SIZEMATRIX*sizeof(double));
  if(distances == NULL){
      printf("No memory available to allocate distances\n");
      exit(0);
  }
  //printf("Matrix Distances Alocated\n");
  
  for(i=0; i<num_points; i++)
  {      
      Particle[i].Volume = 0;
  }
  printf("Hasta aqui todo bien\n");
  for(i=0; i<NFACETS; i++)
  {        
      fill_matrix(distances, facets, coords, i);
      
      Vi = compute_det(distances);
      if(Vi < 0)
      {
	  Vi = 0.0;
      }
      Volumenes[i] = sqrt(Vi/288.0);
      if (i == 0 || i==NFACETS-1) {
          printf("Volumen %d = %f\n", i, Volumenes[i]);
      }
      
  }
  
  /*to speed up*/
  
  int *linear_facets = (int *) malloc((size_t) (4*NFACETS)*sizeof(int));
  int *index = (int *) malloc((size_t) (4*NFACETS)*sizeof(int));
  
  k=0;
  for(i=0; i<NFACETS; i++)
  {
      for(j=0; j<4; j++)
      {
	  linear_facets[k] = facets[i][j];
	  index[k]=i;
	  k++;
      }
  }
  
  gsl_int_int_sort(4*NFACETS,linear_facets,index);
  
  /*end to speed up*/
  
  int pos,start,jfacet;
  for(k=0; k<num_points; k++)
  {
      
      Particle[k].Volume = 0.0;
      counter = 0;
      
      pos = mybsearch(linear_facets,4*NFACETS,k);
      
      //////////////////////////////////////////////////////
      //looking for the starting point
      
      if(pos != EMPTY_FLAG)
      {
	  
	  start=pos;
	  for(i=pos; i>=0; i--) 
	  {
	      if(linear_facets[i] != k)
	      {
		  start=i+1;
		  break;
	      }
	      else
		  start=0;
	  }
	  
	  //now adding the volumes of pieces with 
	  i=start;
	  while((linear_facets[i] == k) && (i < 4*NFACETS))
	  {
	      
	      jfacet = index[i];
	      Particle[k].Volume += Volumenes[jfacet];
        counter++;
	      i++;
	  }
	  
      }
      //////////////////////////////////////////////////////
      
      if(counter > 0)
	      {Particle[k].Volume = Particle[k].Volume/4.0;
        }
      else
      {
	  Particle[k].Volume = 0;
      }
      
  }
  
  free(linear_facets);
  free(index);
  
  for(i=0; i<NFACETS; i++)
      free(facets[i]);
  
  free(facets);
  
  for(i=0; i<num_points; i++)
      free(coords[i]);
  
  free(coords);
  free(distances);
  free(Volumenes);
  return 0;
  
}
