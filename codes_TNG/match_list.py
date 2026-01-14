import numpy as np
import illustris_python as il
import h5py
import os
from scipy.spatial import cKDTree

# ============ CONFIGURACIÓN ============
SIMULATION = "TNG50-1" 
SIMULATION3 = "TNG50-3"
BASE_PATH = "../"
BASE_PATH_OUTPUT = "../TNG50-1/output/"
BASE_PATH_OUTPUT_3 = "../TNG50-3/output/"
SNAP_Z0 = 99    # z=0
BOXSIZE = 35000.0  # ckpc/h 
subhalosdisk_file = "./halos_nodisk.txt"

subhalosdisk = np.loadtxt(subhalosdisk_file, dtype=int)
# =======================================

#para el tiempo de ejecución
import time

time_start = time.time()

def periodic_distance(pos1, pos2, boxsize):
    """
    Calcula la distancia mínima entre dos posiciones considerando 
    condiciones periódicas de frontera.
    """
    delta = np.abs(pos1 - pos2)
    # Si la diferencia es mayor que la mitad de la caja, usar el camino corto
    delta = np.where(delta > 0.5 * boxsize, boxsize - delta, delta)
    return np.linalg.norm(delta, axis=-1)

def get_particle_positions_at_ic(basePath, subhalo_list):
    """
    Obtiene las posiciones en z=20.05 (snap 0) de las partículas DM
    que pertenecen a un subhalo en z=0 (snap=99).
    
    Solo usa el 10% de partículas más ligadas (menor potencial).
    Usa acceso directo por índice: los ParticleIDs de DM son secuenciales,
    entonces el ID corresponde directamente al índice (ID - 1).
    """

    print(f"Cargando partículas DM de {len(subhalo_list)} subhalos en z=0 (snap {SNAP_Z0})...")

    parts_all_sub = []
    offsets = []

    for k, subhalo_id in enumerate(subhalo_list):
        
        print(f"  Cargando particulas DM de subhalo {k+1}/{len(subhalo_list)} (ID: {subhalo_id})")
        print(f"    Tiempo transcurrido: {(time.time() - time_start) / 60:.2f} minutos")

        cutout = il.snapshot.loadSubhalo(basePath, SNAP_Z0, subhalo_id, 'dm',
                                            fields=['ParticleIDs', 'Potential'])
        all_ids = cutout['ParticleIDs']
        potential_z0_all = cutout['Potential']
        
        # Seleccionar solo el 10% más ligado (menor potencial)
        sorted_indices = np.argsort(potential_z0_all)
        n_most_bound = int(len(all_ids) * 0.1)
        most_bound_indices = sorted_indices[:n_most_bound]
        
        dm_ids = all_ids[most_bound_indices]

        parts_all_sub.extend(dm_ids)
        offsets.append(len(dm_ids))

    
    # Cargar posiciones del snapshot 0 (z=20.05) buscando por IDs
    snap_z20 = 0
    
    # Acceso directo al HDF5
    snap_path = os.path.join(basePath, f"snapdir_{snap_z20:03d}")

    
    # Usar set de ParticleIDs para búsqueda O(1)
    ids_needed = set(parts_all_sub)
        
    # Diccionario para mapear ParticleID -> coordenadas
    coords_dict = {}
    
    # Leer archivos HDF5 del snapshot
    chunk_files = sorted([f for f in os.listdir(snap_path) if f.startswith('snap_') and f.endswith('.hdf5')])
    
    print(f"  Snapshot {snap_z20} dividido en {len(chunk_files)} archivos")
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
            if (i + 1) % 10 == 0 or i == len(chunk_files) - 1:
                print(f"    Progreso: {i+1}/{len(chunk_files)} archivos, {len(coords_dict)} partículas encontradas, {(time.time() - time_start)/60:.2f} minutos")
            
            # Early exit si ya encontramos todas
            if len(coords_dict) >= len(ids_needed):
                print(f"  Todas las partículas encontradas (archivo {i+1}/{len(chunk_files)})")
                break
    
    # Reconstruir array de coordenadas en el mismo orden que dm_ids
    positions_z20 = np.array([coords_dict[pid] for pid in parts_all_sub], dtype=np.float32)
    
    print(f"  {len(coords_dict)} posiciones obtenidas del snapshot {snap_z20}")
    print(f"  Tiempo transcurrido: {(time.time() - time_start) / 60:.2f} minutos")

    return parts_all_sub, positions_z20, offsets


def sphere_around_set(points, boxsize):
    """
    Calcula el centro y radio de la esfera mínima que contiene todos los puntos dados,
    considerando condiciones periódicas de frontera.
    """
    if len(points) == 0:
        raise ValueError("No se pueden calcular esfera para un conjunto vacío de puntos")
    
    center = np.mean(points, axis=0)
    # Calcular distancias periódicas al centro
    radii = periodic_distance(points, center, boxsize)
    radius = np.max(radii)
    return center, radius

def find_particles_in_sphere_TNG503(center, radius, tree, ids_all):
    """
    Encuentra ids de partículas DM dentro de una esfera dada en z=20.05 (snap 0)
    de la simulación TNG50-3, usando un KDTree pre-calculado.
    """
    # Usar KDTree para buscar índices de partículas dentro del radio
    # query_ball_point devuelve una lista de índices
    indices = tree.query_ball_point(center, radius)
    
    # Obtener IDs de partículas
    particle_ids_in_sphere = ids_all[indices]
    
    # print(f"  {len(particle_ids_in_sphere):,} partículas encontradas dentro de la esfera")

    return particle_ids_in_sphere

def cost_function(distance, mass_orig, mass_cand, radius, spin_orig, spin_cand, veldisp_orig, veldisp_cand):
    """
    Calcula una función de coste basada en distancia, diferencia de masas, spin y dispersión de velocidades.
    Coste = sqrt((distance/Rhalfmass)^2 + (diferencia de masas/Msubhalo)^2 + (diferencia de spin/Spin_orig)^2 + (diferencia de veldisp/VelDisp_orig)^2)
    """
    if radius <= 0 or mass_orig <= 0:
        return float('inf')
        
    mass_diff = mass_orig - mass_cand
    
    # Spin magnitude
    spin_mag_orig = np.linalg.norm(spin_orig)
    spin_mag_cand = np.linalg.norm(spin_cand)
    spin_diff = spin_mag_orig - spin_mag_cand
    
    # VelDisp
    veldisp_diff = veldisp_orig - veldisp_cand
    
    # Avoid division by zero if original values are 0
    term_mass = (mass_diff / mass_orig)**2
    
    term_spin = 0
    if spin_mag_orig > 0:
        term_spin = (spin_diff / spin_mag_orig)**2
        
    term_veldisp = 0
    if veldisp_orig > 0:
        term_veldisp = (veldisp_diff / veldisp_orig)**2
        
    return np.sqrt((distance/radius)**2 + term_mass + term_spin + term_veldisp)


def main():
    
    # Path de la simulación
    sim_path = f"{BASE_PATH}{SIMULATION}/output"
    
    # Obtener IDs y posiciones iniciales
    particle_ids, positions_ic, offsets = get_particle_positions_at_ic(sim_path, subhalosdisk)
    print(f"Obtenidas posiciones iniciales de partículas DM de los subhalos de interes en z=20.05")
    print(f"Tiempo transcurrido: {(time.time() - time_start) / 60:.2f} minutos")

    # --- PRE-CARGA DE DATOS TNG50-3 ---
    print("\nCargando datos globales de TNG50-3 para optimización...")
    sim_path3 = f"{BASE_PATH}{SIMULATION3}/output"
    
    # 1. Cargar partículas DM de TNG50-3 en z=20.05 para búsqueda espacial
    snap_z20 = 0
    print(f"  Cargando snapshot {snap_z20} (z=20.05)...")
    snap_data_z20 = il.snapshot.loadSubset(sim_path3, snap_z20, 'dm', fields=['ParticleIDs', 'Coordinates'])
    tng50_3_ids_z20 = snap_data_z20['ParticleIDs']
    tng50_3_coords_z20 = snap_data_z20['Coordinates']
    
    print("  Construyendo KDTree para búsqueda espacial...")
    # Construir KDTree con condiciones periódicas
    tree = cKDTree(tng50_3_coords_z20, boxsize=BOXSIZE*1.001)
    
    # 2. Cargar datos de subhalos de TNG50-3 en z=0 para matching
    snap_z0 = 99
    print(f"  Cargando catálogo de subhalos (z=0)...")
    subhalos_tng50_3 = il.groupcat.loadSubhalos(sim_path3, snap_z0, fields=['SubhaloPos', 'SubhaloMass', 'SubhaloSpin', 'SubhaloVelDisp'])
    subhalo_pos_z0 = subhalos_tng50_3['SubhaloPos']
    subhalo_mass_z0 = subhalos_tng50_3['SubhaloMass']
    subhalo_spin_z0 = subhalos_tng50_3['SubhaloSpin']
    subhalo_veldisp_z0 = subhalos_tng50_3['SubhaloVelDisp']
    
    print(f"Datos globales cargados. Tiempo: {(time.time() - time_start) / 60:.2f} minutos")
    # ----------------------------------

    # Convertir offsets de longitudes a índices acumulativos
    cumulative_offsets = np.cumsum([0] + offsets)
    
    # Abrir archivo de salida para escribir en tiempo real
    output_file = "sphere_match_list_results.txt"
    with open(output_file, 'w') as f:
        f.write("id id_match mass_original mass_match distance_ckpc/h shared_fraction\n")
    
    for i, subhalo_id in enumerate(subhalosdisk):
        
        print(f"\nProcesando subhalo {i+1}/{len(subhalosdisk)} (ID: {subhalo_id})")
        print(f"  Tiempo transcurrido: {(time.time() - time_start) / 60:.2f} minutos")

        # Cargar info del subhalo original para diagnóstico
        subhalo_info = il.groupcat.loadSingle(sim_path, SNAP_Z0, subhaloID=subhalo_id)
        subhalo_mass = subhalo_info['SubhaloMass']
        subhalo_pos = subhalo_info['SubhaloPos']
        subhalo_rhalf = subhalo_info['SubhaloHalfmassRad']
        subhalo_spin = subhalo_info['SubhaloSpin']
        subhalo_veldisp = subhalo_info['SubhaloVelDisp']
    
        # print(f"\nCalculando esfera mínima que contiene las partículas...")
        start_idx = cumulative_offsets[i]
        end_idx = cumulative_offsets[i + 1]
        center, radius = sphere_around_set(positions_ic[start_idx:end_idx], BOXSIZE)

        # Usar versión optimizada con KDTree
        particle_ids_in_sphere = find_particles_in_sphere_TNG503(center, radius, tree, tng50_3_ids_z20)

        # Usar versión optimizada de matching - obtener top 5
        # Filtrar subhalos cercanos
        search_radius = 2000.0
        dists = periodic_distance(subhalo_pos_z0, subhalo_pos, BOXSIZE)
        candidate_indices = np.where(dists < search_radius)[0]
        
        if len(candidate_indices) == 0:
            # Sin candidatos cercanos, escribir línea con -1
            with open(output_file, 'a') as f:
                f.write(f"{subhalo_id} -1 {subhalo_mass} 0.0 -1.0 0.0\n")
            print(f"  Sin matches cercanos")
            continue
        
        # Contar partículas compartidas con cada candidato
        ids_target = set(particle_ids_in_sphere)
        subhalo_counts = {}
        
        for sub_id in candidate_indices:
            if subhalo_mass_z0[sub_id] == 0:
                continue
            try:
                dm_ids = il.snapshot.loadSubhalo(sim_path3, SNAP_Z0, sub_id, 'dm', 
                                                fields=['ParticleIDs'])
                matches_in_subhalo = ids_target.intersection(set(dm_ids))
                if len(matches_in_subhalo) > 0:
                    subhalo_counts[sub_id] = len(matches_in_subhalo)
            except Exception:
                continue
        
        if not subhalo_counts:
            # Sin matches con partículas compartidas
            with open(output_file, 'a') as f:
                f.write(f"{subhalo_id} -1 {subhalo_mass} 0.0 -1.0 0.0\n")
            print(f"  Sin matches con partículas compartidas")
            continue
        
        # Ordenar por número de partículas compartidas y tomar top 5
        sorted_by_particles = sorted(subhalo_counts.items(), key=lambda x: x[1], reverse=True)
        top_5_candidates = sorted_by_particles[:5]
        
        # Construir info de candidatos con distancias y coste
        candidate_info = []
        for sub_id, count in top_5_candidates:
            distance = dists[sub_id]
            candidate_mass = subhalo_mass_z0[sub_id]
            candidate_spin = subhalo_spin_z0[sub_id]
            candidate_veldisp = subhalo_veldisp_z0[sub_id]
            fraction = count / len(ids_target) if len(ids_target) > 0 else 0.0
            
            cost = cost_function(distance, subhalo_mass, candidate_mass, subhalo_rhalf, subhalo_spin, candidate_spin, subhalo_veldisp, candidate_veldisp)
            
            candidate_info.append({
                'sub_id': sub_id,
                'count': count,
                'distance': distance,
                'mass': candidate_mass,
                'fraction': fraction,
                'cost': cost
            })
        
        # Ordenar por coste (el menor coste será el mejor match)
        candidate_info.sort(key=lambda x: x['cost'])
        
        # Escribir top 5 al archivo (primero el mejor match)
        with open(output_file, 'a') as f:
            for cand in candidate_info:
                f.write(f"{subhalo_id} {cand['sub_id']} {subhalo_mass} {cand['mass']} {cand['distance']} {cand['fraction']:.4f}\n")
        
        print(f"  Best match: ID={candidate_info[0]['sub_id']}, mass={candidate_info[0]['mass']:.2e}, dist={candidate_info[0]['distance']:.2f}, cost={candidate_info[0]['cost']:.2f}, shared={candidate_info[0]['fraction']:.2%}")
        print(f"  Top 5 matches guardados")

    print(f"\nTiempo total de ejecución: {(time.time() - time_start) / 60:.2f} minutos")
    print(f"Resultados guardados en: {output_file}")
    

if __name__ == "__main__":
    main()
