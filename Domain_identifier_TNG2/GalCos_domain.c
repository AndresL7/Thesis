#include "GalCos_variables.h"

extern int task;  // Variable MPI definida en Domain_identifier.c

static float BoxHalf;

// Grid espacial para optimizar búsqueda de halos
#define GRID_SIZE 32  // Número de celdas por dimensión (ajustable)

struct grid_cell {
    int *halo_indices;  // Índices de halos en esta celda
    int num_halos;      // Número de halos en esta celda
    int capacity;       // Capacidad del array
};

static struct grid_cell ***spatial_grid = NULL;  // Grid 3D
static float grid_cell_size;
static int *top_massive_halo_indices = NULL;  // Índices de los halos más masivos
static int num_top_massive = 0;

#define TOP_MASSIVE_FRACTION 0.30  // Fracción de halos más masivos a revisar siempre (30%)

static float ngb_periodic(float x)
{
  while(x > BoxHalf)
    x-=BoxSize;
  while(x < -BoxHalf)
    x+=BoxSize;
  return x;
}


static float periodic_distance(float *posa, float *posb)
{
  
  float dist,Npos[3],Plenght;
  
  dist = distance(posa[0],posa[1],posa[2],posb[0],posb[1],posb[2]);
  
  if(dist <= 0.5*BoxSize)
    return dist;
  
  // Including corrections for periodic boundary conditions /
  
  Plenght = pow(0.5*BoxSize,2);
  
  Npos[0] = posb[0];
  Npos[1] = posb[1];
  Npos[2] = posb[2];
  
  if(pow(posb[0]-posa[0],2) > Plenght)
    {
      if(posa[0] > posb[0])
	Npos[0] = posb[0] + BoxSize;
      else
	Npos[0] = posb[0] - BoxSize;
    }
  
  if(pow(posb[1]-posa[1],2) > Plenght)
    {
      if(posa[1] > posb[1])
	Npos[1] = posb[1] + BoxSize;
      else
	Npos[1] = posb[1] - BoxSize;
    }
  
  if(pow(posb[2]-posa[2],2) > Plenght)
    {
      if(posa[2] > posb[2])
	Npos[2] = posb[2] + BoxSize;
      else
	Npos[2] = posb[2] - BoxSize;
    }
  
  dist = distance(posa[0],posa[1],posa[2],Npos[0],Npos[1],Npos[2]);
  
  return dist;
}

// Comparador para ordenar halos por masa (de mayor a menor)
static int compare_halo_mass(const void *a, const void *b)
{
  int idx_a = *(int*)a;
  int idx_b = *(int*)b;
  
  if(Halos[idx_b].Mvir > Halos[idx_a].Mvir) return 1;
  if(Halos[idx_b].Mvir < Halos[idx_a].Mvir) return -1;
  return 0;
}

// Función para construir el grid espacial con los halos
void build_halo_spatial_grid(void)
{
  int i, ix, iy, iz;
  int ix_min, ix_max, iy_min, iy_max, iz_min, iz_max;
  int cx, cy, cz;
  
  grid_cell_size = BoxSize / GRID_SIZE;
  
  // Crear lista de índices de halos ordenados por masa
  int *sorted_indices = (int*) malloc(NCLUSTERS * sizeof(int));
  for(i=0; i<NCLUSTERS; i++) {
    sorted_indices[i] = i;
  }
  qsort(sorted_indices, NCLUSTERS, sizeof(int), compare_halo_mass);
  
  // Guardar el TOP_MASSIVE_FRACTION % más masivos
  num_top_massive = (int)(NCLUSTERS * TOP_MASSIVE_FRACTION);
  if(num_top_massive < 1) num_top_massive = 1;
  if(num_top_massive > NCLUSTERS) num_top_massive = NCLUSTERS;
  
  top_massive_halo_indices = (int*) malloc(num_top_massive * sizeof(int));
  for(i=0; i<num_top_massive; i++) {
    top_massive_halo_indices[i] = sorted_indices[i];
  }
  free(sorted_indices);
  
  if(task == 0) {
    printf(" >> Building spatial grid (optimized with top massive halos)...\n");
    printf(" >> Grid: %dx%dx%d cells, cell_size=%.2f\n", GRID_SIZE, GRID_SIZE, GRID_SIZE, grid_cell_size);
    printf(" >> Will always check top %d most massive halos (%.0f%% of total)\n", 
           num_top_massive, TOP_MASSIVE_FRACTION*100.0);
    fflush(stdout);
  }
  
  // Alocar memoria para el grid 3D
  spatial_grid = (struct grid_cell ***) malloc(GRID_SIZE * sizeof(struct grid_cell **));
  for(ix=0; ix<GRID_SIZE; ix++) {
    spatial_grid[ix] = (struct grid_cell **) malloc(GRID_SIZE * sizeof(struct grid_cell *));
    for(iy=0; iy<GRID_SIZE; iy++) {
      spatial_grid[ix][iy] = (struct grid_cell *) malloc(GRID_SIZE * sizeof(struct grid_cell));
      for(iz=0; iz<GRID_SIZE; iz++) {
        spatial_grid[ix][iy][iz].halo_indices = NULL;
        spatial_grid[ix][iy][iz].num_halos = 0;
        spatial_grid[ix][iy][iz].capacity = 0;
      }
    }
  }
  
  // Insertar cada halo en todas las celdas que toca su esfera de influencia (Rvir)
  for(i=0; i<NCLUSTERS; i++) {
    
    if(Halos[i].Rvir <= 0.0) continue;
    
    ix_min = (int) floorf((Halos[i].pos[0] - Halos[i].Rvir) / grid_cell_size);
    ix_max = (int) ceilf((Halos[i].pos[0] + Halos[i].Rvir) / grid_cell_size);
    iy_min = (int) floorf((Halos[i].pos[1] - Halos[i].Rvir) / grid_cell_size);
    iy_max = (int) ceilf((Halos[i].pos[1] + Halos[i].Rvir) / grid_cell_size);
    iz_min = (int) floorf((Halos[i].pos[2] - Halos[i].Rvir) / grid_cell_size);
    iz_max = (int) ceilf((Halos[i].pos[2] + Halos[i].Rvir) / grid_cell_size);
    
    // Agregar este halo a todas las celdas que toca
    for(ix=ix_min; ix<=ix_max; ix++) {
      for(iy=iy_min; iy<=iy_max; iy++) {
        for(iz=iz_min; iz<=iz_max; iz++) {
          
          // Manejar condiciones de frontera periódicas
          cx = ((ix % GRID_SIZE) + GRID_SIZE) % GRID_SIZE;
          cy = ((iy % GRID_SIZE) + GRID_SIZE) % GRID_SIZE;
          cz = ((iz % GRID_SIZE) + GRID_SIZE) % GRID_SIZE;
          
          // Expandir array si es necesario
          if(spatial_grid[cx][cy][cz].num_halos >= spatial_grid[cx][cy][cz].capacity) {
            spatial_grid[cx][cy][cz].capacity = (spatial_grid[cx][cy][cz].capacity == 0) ? 4 : spatial_grid[cx][cy][cz].capacity * 2;
            spatial_grid[cx][cy][cz].halo_indices = (int *) realloc(spatial_grid[cx][cy][cz].halo_indices, 
                                                                      spatial_grid[cx][cy][cz].capacity * sizeof(int));
          }
          
          // Agregar índice del halo
          spatial_grid[cx][cy][cz].halo_indices[spatial_grid[cx][cy][cz].num_halos] = i;
          spatial_grid[cx][cy][cz].num_halos++;
        }
      }
    }
  }
  
  if(task == 0) {
    printf(" >> Spatial grid built successfully\n");
    fflush(stdout);
  }
}


// Función para liberar memoria del grid
void free_halo_spatial_grid(void)
{
  int ix, iy, iz;
  
  if(spatial_grid == NULL) return;
  
  for(ix=0; ix<GRID_SIZE; ix++) {
    for(iy=0; iy<GRID_SIZE; iy++) {
      for(iz=0; iz<GRID_SIZE; iz++) {
        if(spatial_grid[ix][iy][iz].halo_indices != NULL) {
          free(spatial_grid[ix][iy][iz].halo_indices);
        }
      }
      free(spatial_grid[ix][iy]);
    }
    free(spatial_grid[ix]);
  }
  free(spatial_grid);
  spatial_grid = NULL;
  
  if(top_massive_halo_indices != NULL) {
    free(top_massive_halo_indices);
    top_massive_halo_indices = NULL;
  }
}


int GalCos_domain(int ipart)
{
  int i, j, counter, halo_idx;
  int ix, iy, iz, dx, dy, dz, cx, cy, cz;
  float dist, mindist;
  
  counter = EMPTY_FLAG;
  mindist = 1e30;
  BoxHalf = 0.5*BoxSize;
  
  // Encontrar la celda de la partícula (usar floor para manejar correctamente negativos)
  ix = (int) floorf(Particle[ipart].pos[0] / grid_cell_size);
  iy = (int) floorf(Particle[ipart].pos[1] / grid_cell_size);
  iz = (int) floorf(Particle[ipart].pos[2] / grid_cell_size);
  
  // Asegurar que estén en rango [0, GRID_SIZE)
  ix = ((ix % GRID_SIZE) + GRID_SIZE) % GRID_SIZE;
  iy = ((iy % GRID_SIZE) + GRID_SIZE) % GRID_SIZE;
  iz = ((iz % GRID_SIZE) + GRID_SIZE) % GRID_SIZE;
  
  // FASE 1: Búsqueda local en 15x15x15 celdas vecinas (xx celdas)
  for(dx=-7; dx<=7; dx++) {
    for(dy=-7; dy<=7; dy++) {
      for(dz=-7; dz<=7; dz++) {
        
        cx = ((ix + dx) % GRID_SIZE + GRID_SIZE) % GRID_SIZE;
        cy = ((iy + dy) % GRID_SIZE + GRID_SIZE) % GRID_SIZE;
        cz = ((iz + dz) % GRID_SIZE + GRID_SIZE) % GRID_SIZE;
        
        for(j=0; j<spatial_grid[cx][cy][cz].num_halos; j++) {
          
          halo_idx = spatial_grid[cx][cy][cz].halo_indices[j];
          
          if(Halos[halo_idx].Rvir <= 0.0) continue;
          
          dist = periodic_distance(Halos[halo_idx].pos, Particle[ipart].pos);
          dist = dist / Halos[halo_idx].Rvir;
          
          if(dist < mindist) {
            mindist = dist;
            counter = halo_idx;
          }
        }
      }
    }
  }
  
  // FASE 2: Revisar siempre los halos más masivos (pueden ganar aunque estén lejos)
  // COMENTADO para probar solo búsqueda local 4x4x4
  /*
  for(i=0; i<num_top_massive; i++) {
    halo_idx = top_massive_halo_indices[i];
    
    if(Halos[halo_idx].Rvir <= 0.0) continue;
    
    dist = periodic_distance(Halos[halo_idx].pos, Particle[ipart].pos);
    dist = dist / Halos[halo_idx].Rvir;
    
    if(dist < mindist) {
      mindist = dist;
      counter = halo_idx;
    }
  }
  */
  
  // Fallback de seguridad: si no encontró ningún halo, buscar en todos
  // Esto solo pasaría si el grid está mal configurado o no hay halos
  if(counter == EMPTY_FLAG) {
    for(i=0; i<NCLUSTERS; i++) {
      // Skip halos with invalid Rvir
      if(Halos[i].Rvir <= 0.0) continue;
      
      dist = periodic_distance(Halos[i].pos, Particle[ipart].pos);
      dist = dist / Halos[i].Rvir;
      
      if(dist < mindist) {
        mindist = dist;
        counter = i;
      }
    }
  }
  
  // Asignar partícula al dominio del halo
  Particle[ipart].Cluster_ID = counter;
  
  Halos[counter].NDomain_particles++;
  Halos[counter].Domain_particles = realloc(Halos[counter].Domain_particles,
                                            (size_t) Halos[counter].NDomain_particles*sizeof(int));
  Halos[counter].Domain_particles[Halos[counter].NDomain_particles-1] = Particle[ipart].Oid;
  
  return 0;
}
