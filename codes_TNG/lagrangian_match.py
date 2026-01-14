import numpy as np
import illustris_python as il
import h5py
import os
import time
# ============ CONFIGURACIÓN ============
SIMULATION = "TNG50-1"  # o "TNG50-3"
BASE_PATH = "../"
BASE_PATH_OUTPUT = "../TNG50-1/output/"
SNAP_Z0 = 99    # z=0
SNAP_Z20 = 0    # z=20.05 (snapdir_000)
# =======================================
time_start = time.time()

def get_particle_positions(basePath, subhalo_id):
    """
    Obtiene las posiciones en z=20.05 (snap 0) de las partículas DM
    que pertenecen a un subhalo en z=0 (snap=99).
    
    Solo usa el 30% de partículas más ligadas (menor potencial).
    """

    print(f"Cargando partículas DM del subhalo {subhalo_id} en z=0 (snap {SNAP_Z0})...")
    print(f"  Ruta base: {basePath}")

    
    cutout = il.snapshot.loadSubhalo(basePath, SNAP_Z0, subhalo_id, 'dm', 
                                        fields=['ParticleIDs', 'Coordinates', 'Potential'])
    
    
    all_ids = cutout['ParticleIDs']
    all_coords = cutout['Coordinates']
    potential_z0_all = cutout['Potential']
    
    print(f"  Encontradas {len(all_ids)} partículas DM en z=0")
    
    # Seleccionar solo el 30% más ligado (menor potencial)
    sorted_indices = np.argsort(potential_z0_all)
    n_most_bound = int(len(all_ids) * 0.3)
    most_bound_indices = sorted_indices[:n_most_bound]
    
    dm_ids = all_ids[most_bound_indices]
    coords_z0 = all_coords[most_bound_indices]
    potential_z0 = potential_z0_all[most_bound_indices]
    
    print(f"  Usando solo el 30% más ligado: {len(dm_ids)} partículas")
    print(f"  Rango de potenciales en z=0: [{potential_z0.min():.6e}, {potential_z0.max():.6e}]")

    # Cargar posiciones del snapshot 0 (z=20.05) directamente del HDF5
    print(f"\nCargando posiciones en z=20.05 (snap {SNAP_Z20})...")
    
    # Acceso directo al HDF5
    snap_path = os.path.join(basePath, f"snapdir_{SNAP_Z20:03d}")
    
    # Usar set de ParticleIDs para búsqueda O(1)
    ids_needed = set(dm_ids)
    
    # Diccionario para mapear ParticleID -> coordenadas
    coords_dict = {}
    
    # Leer archivos HDF5 del snapshot
    chunk_files = sorted([f for f in os.listdir(snap_path) if f.startswith('snap_') and f.endswith('.hdf5')])
    
    print(f"  Snapshot {SNAP_Z20} dividido en {len(chunk_files)} archivos")
    print(f"  Buscando {len(ids_needed)} partículas específicas...")
    
    for i, chunk_file in enumerate(chunk_files):
        chunk_path = os.path.join(snap_path, chunk_file)
        
        with h5py.File(chunk_path, 'r') as f:
            # Cargar IDs y coordenadas de este chunk
            ids_chunk = f['PartType1/ParticleIDs'][:]
            coords_chunk = f['PartType1/Coordinates'][:]
            
            # Buscar qué IDs de este chunk están en los que necesitamos
            mask = np.isin(ids_chunk, list(ids_needed))
            
            if np.any(mask):
                # Guardar las coordenadas con su ParticleID como clave
                matching_ids = ids_chunk[mask]
                matching_coords = coords_chunk[mask]
                
                for pid, coord in zip(matching_ids, matching_coords):
                    coords_dict[pid] = coord
            
            # Mostrar progreso
            if (i + 1) % 20 == 0 or i == len(chunk_files) - 1:
                print(f"    Progreso: {i+1}/{len(chunk_files)} archivos, {len(coords_dict)} partículas encontradas, {(time.time() - time_start)/60:.2f} minutos")
            
            # Early exit si ya encontramos todas
            if len(coords_dict) >= len(ids_needed):
                print(f"  Todas las partículas encontradas (archivo {i+1}/{len(chunk_files)})")
                break
    
    # Reconstruir array de coordenadas en el mismo orden que dm_ids
    coords_z20 = np.array([coords_dict[pid] for pid in dm_ids], dtype=np.float32)
    
    print(f"  {len(coords_dict)} posiciones obtenidas del snapshot {SNAP_Z20}")

    return dm_ids, coords_z20, coords_z0, potential_z0


def main():
    # Pedir el subhalo ID al usuario
    subhalo_id = int(input("Ingresa el ID del subhalo: "))
    
    # Path de la simulación
    sim_path = f"{BASE_PATH}/{SIMULATION}/output"
    
    # Obtener IDs, posiciones y potencial lagrangianos
    particle_ids, positions_z20, positions_z0, potential_z0 = get_particle_positions(sim_path, subhalo_id)
    
    print(f"\n{'='*60}")
    print(f"RESULTADOS")
    print(f"Simulación: {SIMULATION}")
    print(f"Subhalo ID en z=0: {subhalo_id}")
    print(f"Partículas encontradas: {len(particle_ids)}")
    print(f"Tiempo total de ejecución: {(time.time() - time_start) / 60:.2f} minutos")
    print(f"{'-'*60}")
    
    # Guardar en archivo CSV
    output_file = f"subhalo_{subhalo_id}_z20.csv"
    print(f"\nGuardando resultados completos en: {output_file}")
    
    with open(output_file, 'w') as f:
        f.write("ParticleID,X_z20,Y_z20,Z_z20,X_z0,Y_z0,Z_z0,Potential_z0\n")
        for pid, pos_z20, pos_z0, pot_z0 in zip(particle_ids, positions_z20, positions_z0, potential_z0):
            f.write(f"{pid},{pos_z20[0]:.6f},{pos_z20[1]:.6f},{pos_z20[2]:.6f},{pos_z0[0]:.6f},{pos_z0[1]:.6f},{pos_z0[2]:.6f},{pot_z0:.6e}\n")
    
    print(f"Archivo guardado exitosamente")

if __name__ == "__main__":
    main()
