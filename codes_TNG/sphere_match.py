import numpy as np
import illustris_python as il
import h5py
import os

# ============ CONFIGURACIÓN ============
SIMULATION = "TNG50-1" 
SIMULATION3 = "TNG50-3"
BASE_PATH = "../"
BASE_PATH_OUTPUT = "../TNG50-1/output/"
BASE_PATH_OUTPUT_3 = "../TNG50-3/output/"
SNAP_Z0 = 99    # z=0
BOXSIZE = 35000.0  # ckpc/h 
# =======================================

#para el tiempo de ejecución
import time

time_start = time.time()

def periodic_distance(pos1, pos2, boxsize):
    """
    Calcula la distancia mínima entre dos posiciones considerando 
    condiciones periódicas de frontera.
    
    Si una partícula está en x=1000 y otra en x=50000 en una caja de 51700,
    la distancia real es 1700 (cruzando el borde), no 49000.
    """
    delta = np.abs(pos1 - pos2)
    # Si la diferencia es mayor que la mitad de la caja, usar el camino corto
    delta = np.where(delta > 0.5 * boxsize, boxsize - delta, delta)
    return np.linalg.norm(delta, axis=-1)

def get_particle_positions_at_ic(basePath, subhalo_id):
    """
    Obtiene las posiciones en z=20.05 (snap 0) de las partículas DM
    que pertenecen a un subhalo en z=0 (snap=99).
    
    Solo usa el 10% de partículas más ligadas (menor potencial).
    Usa acceso directo por índice: los ParticleIDs de DM son secuenciales,
    entonces el ID corresponde directamente al índice (ID - 1).
    """

    print(f"Cargando partículas DM del subhalo {subhalo_id} en z=0 (snap {SNAP_Z0})...")

    # Obtener IDs, coordenadas y potencial de partículas DM del subhalo en z=0
    cutout = il.snapshot.loadSubhalo(basePath, SNAP_Z0, subhalo_id, 'dm', 
                                     fields=['ParticleIDs', 'Potential'])
    all_ids = cutout['ParticleIDs']
    potential_z0_all = cutout['Potential']
    
    print(f"  Encontradas {len(all_ids)} partículas DM en z=0")
    
    # Seleccionar solo el 10% más ligado (menor potencial)
    sorted_indices = np.argsort(potential_z0_all)
    n_most_bound = int(len(all_ids) * 0.1)
    most_bound_indices = sorted_indices[:n_most_bound]
    
    dm_ids = all_ids[most_bound_indices]
    
    print(f"  Usando solo el 10% más ligado: {len(dm_ids)} partículas")

    # Cargar posiciones del snapshot 0 (z=20.05) buscando por IDs
    snap_z20 = 0
    print(f"\nCargando posiciones en z=20.05 (snap {snap_z20})...")
    
    # Acceso directo al HDF5
    snap_path = os.path.join(basePath, f"snapdir_{snap_z20:03d}")
    
    # Usar set de ParticleIDs para búsqueda O(1)
    ids_needed = set(dm_ids)
    
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
            if (i + 1) % 20 == 0 or i == len(chunk_files) - 1:
                print(f"    Progreso: {i+1}/{len(chunk_files)} archivos, {len(coords_dict)} partículas encontradas, {(time.time() - time_start)/60:.2f} minutos")
            
            # Early exit si ya encontramos todas
            if len(coords_dict) >= len(ids_needed):
                print(f"  ✓ Todas las partículas encontradas (archivo {i+1}/{len(chunk_files)})")
                break
    
    # Reconstruir array de coordenadas en el mismo orden que dm_ids
    positions_z20 = np.array([coords_dict[pid] for pid in dm_ids], dtype=np.float32)
    
    print(f"  ✓ {len(coords_dict)} posiciones obtenidas del snapshot {snap_z20}")
    print(f"  Tiempo transcurrido: {(time.time() - time_start) / 60:.2f} minutos")

    return dm_ids, positions_z20


def sphere_around_set(points, boxsize):
    """
    Calcula el centro y radio de la esfera mínima que contiene todos los puntos dados,
    considerando condiciones periódicas de frontera.
    """
    center = np.mean(points, axis=0)
    # Calcular distancias periódicas al centro
    radii = periodic_distance(points, center, boxsize)
    radius = np.max(radii)
    return center, radius

def find_particles_in_sphere_TNG503(center, radius, boxsize):
    """
    Encuentra ids de partículas DM dentro de una esfera dada en z=20.05 (snap 0)
    de la simulación TNG50-3, considerando condiciones periódicas de frontera.
    """
    snap_z20 = 0
    print(f"\nBuscando partículas en TNG50-3 dentro de la esfera (snap {snap_z20}, z=20.05)...")

    # Path de la simulación TNG50-3
    sim_path3 = f"{BASE_PATH}{SIMULATION3}/output"
    
    # Cargar coordenadas e IDs del snapshot 0 usando illustris_python
    snap_data = il.snapshot.loadSubset(sim_path3, snap_z20, 'dm', 
                                       fields=['ParticleIDs', 'Coordinates'])
    ids_all = snap_data['ParticleIDs']
    coords_all = snap_data['Coordinates']

    print(f"  Snapshot {snap_z20} TNG50-3 contiene {len(ids_all):,} partículas DM totales")
    print(f"  Buscando partículas dentro de la esfera...")

    # Calcular distancias PERIÓDICAS de todas las partículas al centro de la esfera
    distances = periodic_distance(coords_all, center, boxsize)
    mask = distances <= radius

    # Obtener IDs de partículas dentro de la esfera
    particle_ids_in_sphere = ids_all[mask]
    
    print(f"  ✓ {len(particle_ids_in_sphere):,} partículas encontradas dentro de la esfera")

    return particle_ids_in_sphere


def find_matching_subhalo_TNG503(basepath, particle_ids_in_sphere, target_position, original_mass):
    """
    Encuentra el subhalo de TNG50-3 (z=0) que contiene la mayor cantidad
    de partículas cuyos IDs están en 'particle_ids_in_sphere'.
    
    Además considera la proximidad espacial al subhalo original para
    desambiguar casos donde múltiples subhalos contienen partículas de la región.
    """
    print("\nBuscando subhalo correspondiente en TNG50-3 (z=0)…")

    snapnum = 99

    # Cargar longitudes y posiciones de subhalos (DM partículas)
    print("  Cargando longitudes y posiciones de subhalos…")
    subhalos = il.groupcat.loadSubhalos(basepath, snapnum,
                                        fields=['SubhaloLenType', 'SubhaloPos'])
    subhalo_len_dm = subhalos['SubhaloLenType'][:, 1]  # columna 1 = DM
    subhalo_positions = subhalos['SubhaloPos']
    n_subhalos = len(subhalo_len_dm)
    print(f"  Hay {n_subhalos} subhalos en TNG50-3 (z=0)")

    # Preparar conteo
    ids_target = set(particle_ids_in_sphere)
    subhalo_counts = {}  # usar diccionario para solo guardar subhalos con matches
    
    # Crear diccionario: particle_id -> subhalo_id
    print("  Construyendo diccionario de partículas -> subhalos...")
    
    particles_checked = 0
    subhalos_with_matches = 0
    
    for sub_id in range(n_subhalos):
        # Solo procesar subhalos con partículas DM
        if subhalo_len_dm[sub_id] == 0:
            continue
        
        # Mostrar progreso cada 20000 subhalos
        if sub_id % 20_000 == 0:
            print(f"    Procesando subhalo {sub_id}/{n_subhalos}, "
                  f"subhalos con coincidencias: {subhalos_with_matches}")
        
        try:
            # Cargar IDs de partículas DM de este subhalo
            dm_ids = il.snapshot.loadSubhalo(basepath, snapnum, sub_id, 'dm', 
                                            fields=['ParticleIDs'])
            
            particles_checked += len(dm_ids)
            
            # Contar cuántas de estas partículas están en particle_ids_in_sphere
            # Usar intersección de conjuntos para eficiencia
            matches_in_subhalo = ids_target.intersection(set(dm_ids))
            
            if len(matches_in_subhalo) > 0:
                subhalo_counts[sub_id] = len(matches_in_subhalo)
                subhalos_with_matches += 1
                
        except Exception as e:
            # Algunos subhalos pueden no tener partículas DM
            continue
    
    total_found = sum(subhalo_counts.values())
    print(f"\n  Total de partículas emparejadas: {total_found} de {len(ids_target)}")
    print(f"  Partículas DM revisadas: {particles_checked:,}")
    print(f"  Subhalos con coincidencias: {subhalos_with_matches}")

    if total_found == 0:
        print("  No se encontraron coincidencias en TNG50-3.")
        return -1, 0, 0.0

    # Ordenar por número de partículas compartidas
    sorted_by_particles = sorted(subhalo_counts.items(), key=lambda x: x[1], reverse=True)
    
    # Tomar los 5 mejores por número de partículas
    top_5_candidates = sorted_by_particles[:5]
    
    print(f"\n  Top 5 subhalos por partículas compartidas:")
    print(f"  {'Rank':<6} {'SubID':<8} {'Count':<7} {'Frac':<7} {'Dist[ckpc/h]':<14} {'MassRatio':<10}")
    print(f"  {'-'*75}")
    
    candidate_info = []
    for i, (sub_id, count) in enumerate(top_5_candidates, 1):
        # Calcular distancia PERIÓDICA entre posiciones en z=0
        distance = periodic_distance(subhalo_positions[sub_id], target_position, BOXSIZE)
        particle_frac = count / len(particle_ids_in_sphere)
        
        # Obtener masa del candidato
        candidate_mass = il.groupcat.loadSingle(basepath, 99, subhaloID=sub_id)['SubhaloMass']
        mass_ratio = original_mass / candidate_mass if candidate_mass > 0 else np.inf
        
        candidate_info.append({
            'sub_id': sub_id,
            'count': count,
            'distance': distance,
            'particle_frac': particle_frac,
            'mass': candidate_mass,
            'mass_ratio': mass_ratio
        })
        print(f"  #{i:<5} {sub_id:<8} {count:<7} {particle_frac:<7.3f} {distance:<14.1f} {mass_ratio:<10.2f}")
    
    # De los top 5, elegir el más cercano
    best_candidate = min(candidate_info, key=lambda x: x['distance'])
    
    print(f"\n  → Mejor match (más cercano de los top 5): Subhalo {best_candidate['sub_id']}")
    
    return best_candidate['sub_id'], best_candidate['count'], best_candidate['particle_frac'], best_candidate['mass']




def main():
    # Pedir el subhalo ID al usuario
    subhalo_id = int(input("Ingresa el ID del subhalo: "))
    
    # Path de la simulación
    sim_path = f"{BASE_PATH}{SIMULATION}/output"
    
    # Obtener IDs y posiciones iniciales
    particle_ids, positions_ic = get_particle_positions_at_ic(sim_path, subhalo_id)

    # Cargar info del subhalo original para diagnóstico
    subhalo_mass = il.groupcat.loadSingle(sim_path, SNAP_Z0, subhaloID=subhalo_id)['SubhaloMass']
    subhalo_pos = il.groupcat.loadSingle(sim_path, SNAP_Z0, subhaloID=subhalo_id)['SubhaloPos']
    
    # Sphere
    print(f"\nCalculando esfera mínima que contiene las partículas...")
    center, radius = sphere_around_set(positions_ic, BOXSIZE)

    print(f"\n{'='*60}")
    print(f"RESULTADOS - DIAGNÓSTICO")
    print(f"Simulación: {SIMULATION}")
    print(f"Subhalo ID en z=0: {subhalo_id}")
    print(f"Subhalo Mass (z=0): {subhalo_mass:.2e} (1e10 Msun/h)")
    print(f"Subhalo Pos (z=0): {subhalo_pos}")
    print(f"Partículas DM encontradas: {len(particle_ids)}")
    print(f"\nESFERA EN z=20.05 (con condiciones periódicas):")
    print(f"  Centro: [{center[0]:.2f}, {center[1]:.2f}, {center[2]:.2f}] ckpc/h")
    print(f"  Radio: {radius:.2f} ckpc/h")
    print(f"  BoxSize: {BOXSIZE:.2f} ckpc/h")
    print(f"  Tiempo transcurrido: {(time.time() - time_start) / 60:.2f} minutos")
    print(f"{'-'*60}")

    particle_ids_in_sphere = find_particles_in_sphere_TNG503(center, radius, BOXSIZE)
    print(f"Partículas DM de TNG50-3 encontradas en la esfera: {len(particle_ids_in_sphere)}")
    print(f" Tiempo transcurrido: {(time.time() - time_start) / 60:.2f} minutos")
    print(f"{'-'*60}")

    # Path de la simulación TNG50-3
    sim_path3 = f"{BASE_PATH}{SIMULATION3}/output"

    best_subhalo, max_count, frac, best_mass = find_matching_subhalo_TNG503(sim_path3, particle_ids_in_sphere, subhalo_pos, subhalo_mass)

    # Cargar info del mejor match para comparación
    best_pos = il.groupcat.loadSingle(sim_path3, SNAP_Z0, subhaloID=best_subhalo)['SubhaloPos']
    
    print(f"{'='*60}")
    print(f"MEJOR MATCH EN TNG50-3:")
    print(f"   Subhalo ID: {best_subhalo}")
    print(f"   Partículas compartidas: {max_count}")
    print(f"   Fracción compartida: {frac:.3f}")
    print(f"   Masa (z=0): {best_mass:.2e} (1e10 Msun/h)")
    print(f"   Posición (z=0): {best_pos}")
    print(f"\nCOMPARACIÓN CON SUBHALO ORIGINAL (TNG50-1):")
    print(f"   Masa original: {subhalo_mass:.2e} (1e10 Msun/h)")
    if best_mass > 0:
        print(f"   Ratio de masas (TNG50-1/TNG50-3): {subhalo_mass/best_mass:.2f}")
    else:
        print(f"   Ratio de masas: N/A (masa TNG50-3 = 0)")
    
    # Calcular distancia periódica entre posiciones en z=0
    dist_periodic = periodic_distance(subhalo_pos, best_pos, BOXSIZE)
    print(f"   Distancia entre posiciones (periódica): {dist_periodic:.2f} ckpc/h")

    print(f"\nTiempo total de ejecución: {(time.time() - time_start) / 60:.2f} minutos")
    
    # # ========== GUARDAR TODAS LAS POSICIONES DE PARTÍCULAS EN ARCHIVOS ==========
    # print(f"\n{'='*60}")
    # print(f"GUARDANDO POSICIONES DE PARTÍCULAS EN ARCHIVOS")
    # print(f"{'='*60}")
    
    # # Posiciones de partículas del subhalo original en TNG50-1 (30% más ligado en z=0)
    # print(f"\nCargando posiciones de partículas del subhalo original en TNG50-1 (z=0)...")
    # cutout_tng1 = il.snapshot.loadSubhalo(sim_path, SNAP_Z0, subhalo_id, 'dm', 
    #                                       fields=['ParticleIDs', 'Coordinates', 'Potential'])
    # all_ids_tng1 = cutout_tng1['ParticleIDs']
    # all_coords_tng1 = cutout_tng1['Coordinates']
    # potential_tng1 = cutout_tng1['Potential']
    
    # # Seleccionar el 30% más ligado
    # sorted_indices_tng1 = np.argsort(potential_tng1)
    # n_most_bound_tng1 = int(len(all_ids_tng1) * 0.3)
    # most_bound_indices_tng1 = sorted_indices_tng1[:n_most_bound_tng1]
    
    # coords_tng1_bound = all_coords_tng1[most_bound_indices_tng1]
    # ids_tng1_bound = all_ids_tng1[most_bound_indices_tng1]
    
    # # Guardar en archivo
    # output_file_tng1 = f'particles_TNG50-1_subhalo_{subhalo_id}.dat'
    # data_tng1 = np.column_stack([ids_tng1_bound, coords_tng1_bound])
    # np.savetxt(output_file_tng1, data_tng1,
    #            header=f'Subhalo {subhalo_id} en TNG50-1 (z=0) - 30% más ligado\nParticleID   X[ckpc/h]   Y[ckpc/h]   Z[ckpc/h]',
    #            fmt='%d %.6f %.6f %.6f')
    # print(f"  Posiciones TNG50-1 guardadas en '{output_file_tng1}'")
    # print(f"  Total de partículas DM (30% más ligado): {len(coords_tng1_bound)}")
    
    # # Posiciones de partículas del mejor match en TNG50-3 (z=0)
    # print(f"\nCargando posiciones de partículas del mejor match en TNG50-3 (z=0)...")
    # cutout_tng3 = il.snapshot.loadSubhalo(sim_path3, SNAP_Z0, best_subhalo, 'dm', 
    #                                       fields=['ParticleIDs', 'Coordinates'])
    # ids_tng3 = cutout_tng3['ParticleIDs']
    # coords_tng3 = cutout_tng3['Coordinates']
    
    # # Guardar en archivo
    # output_file_tng3 = f'particles_TNG50-3_subhalo_{best_subhalo}.dat'
    # data_tng3 = np.column_stack([ids_tng3, coords_tng3])
    # np.savetxt(output_file_tng3, data_tng3,
    #            header=f'Subhalo {best_subhalo} en TNG50-3 (z=0) - Match del subhalo {subhalo_id} de TNG50-1\nParticleID   X[ckpc/h]   Y[ckpc/h]   Z[ckpc/h]',
    #            fmt='%d %.6f %.6f %.6f')
    # print(f"  Posiciones TNG50-3 guardadas en '{output_file_tng3}'")
    # print(f"  Total de partículas DM: {len(coords_tng3)}")

if __name__ == "__main__":
    main()
