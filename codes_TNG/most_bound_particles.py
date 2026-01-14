#!/usr/bin/env python3
import numpy as np
import illustris_python as il
import h5py
import os
import time

# ============ CONFIGURACIÓN ============
SIMULATION = "TNG50-1" 
BASE_PATH = "../"
BASE_PATH_OUTPUT = "../TNG50-1/output/"
SNAP_Z0 = 99    # z=0
# =======================================

time_start = time.time()


def get_most_bound_particles(basePath, subhalo_id, fraction=0.3):
    """
    Obtiene los IDs de las partículas más ligadas de un subhalo.
    
    """
    
    print(f"Cargando partículas DM del subhalo {subhalo_id}...")
    cutout = il.snapshot.loadSubhalo(basePath, SNAP_Z0, subhalo_id, 'dm', 
                                     fields=['ParticleIDs', 'Potential'])
    
    particle_ids = cutout['ParticleIDs']
    potential = cutout['Potential']
    
    n_particles = len(particle_ids)
    print(f"  Total: {n_particles:,} partículas")
    
    # Ordenar por potencial (más negativo = más ligado)
    sorted_indices = np.argsort(potential)
    
    # Tomar el top X% más ligado
    n_most_bound = int(n_particles * fraction)
    most_bound_indices = sorted_indices[:n_most_bound]
    most_bound_ids = particle_ids[most_bound_indices]
    
    print(f"  Seleccionadas: {n_most_bound:,} más ligadas ({fraction:.1%})")
    
    return most_bound_ids


def save_results(subhalo_id, most_bound_ids, fraction):
    """
    Guarda los IDs de partículas más ligadas.
    """
    output_file = f"{subhalo_id}_most_bound_{int(fraction*100)}.txt"
    np.savetxt(output_file, most_bound_ids, fmt='%d')
    print(f"\n Guardado: {output_file}")


def main():
    subhalo_id = int(input("ID del subhalo: "))
    
    fraction = 0.3  # 30% más ligadas
    
    sim_path = f"{BASE_PATH}{SIMULATION}/output"
    
    most_bound_ids = get_most_bound_particles(sim_path, subhalo_id, fraction)
    
    save_results(subhalo_id, most_bound_ids, fraction)
    
    print(f"Tiempo: {(time.time() - time_start) / 60:.2f} minutos")


if __name__ == "__main__":
    main()
