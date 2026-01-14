import illustris_python as il
import numpy as np
import os

# CONFIGURACIÓN 
basePath = '../TNG50-1/output' 
snap_z0 = 99  # snapshot de z=0 para TNG50-1

print(f"Cargando datos de TNG50-1 desde: {basePath}")

# CARGAR DATOS DE SUBHALOS EN Z=0
print("\nCargando subhalos en z=0...")
subhalos_z0 = il.groupcat.loadSubhalos(basePath, snap_z0, fields=['SubhaloGrNr', 'SubhaloMass'])

# CARGAR MASAS DE HALOS FOF EN Z=0
print("\nCargando masas de halos FOF en z=0...")
# Al pedir solo 1 campo, devuelve un array directamente
halos_masses = il.groupcat.loadHalos(basePath, snap_z0, fields=['GroupMass'])
print(f"{len(halos_masses)} halos FOF cargados")

# CARGAR LISTAS DE HALOS
print("\nCargando listas de halos (halos.txt y halos_nodisk.txt)...")
ids_disk = np.loadtxt('halos.txt', dtype=int)
if ids_disk.ndim == 0: ids_disk = np.array([ids_disk])

ids_nodisk = np.loadtxt('halos_nodisk.txt', dtype=int)
if ids_nodisk.ndim == 0: ids_nodisk = np.array([ids_nodisk])

# Combinar listas
target_ids = np.concatenate([ids_disk, ids_nodisk])
has_disk = np.concatenate([np.ones(len(ids_disk), dtype=int), np.zeros(len(ids_nodisk), dtype=int)])
n_targets = len(target_ids)

print(f"Total de subhalos a procesar: {n_targets} ({len(ids_disk)} con disco, {len(ids_nodisk)} sin disco)")

# ARRAYS DE RESULTADOS
snap_indep = np.full(n_targets, -1)
res_fof_indep = np.full(n_targets, -1)
res_fof_z0 = np.full(n_targets, -1)
res_mass_z0 = np.full(n_targets, np.nan)
res_mass_fof_z0 = np.full(n_targets, np.nan)

# LOOP PRINCIPAL
print("\nProcesando árboles de merger (Rama Principal - MPB)...")
processed = 0

for i, sid in enumerate(target_ids):

    # Obtener datos en z=0
    try:
        fof_id = subhalos_z0['SubhaloGrNr'][sid]
        mass = subhalos_z0['SubhaloMass'][sid]
        
        res_fof_z0[i] = fof_id
        res_mass_z0[i] = mass
        
        if fof_id >= 0 and fof_id < len(halos_masses):
            res_mass_fof_z0[i] = halos_masses[fof_id]
    except IndexError:
        print(f"Warning: SubhaloID {sid} fuera de rango.")
        continue

    # Mostrar progreso
    if i % 100 == 0:
        print(f"Procesando {i+1}/{n_targets} (SubhaloID: {sid})...")
    
    # Cargar árbol de merger (Rama Principal) del subhalo con id sid (onlyMPB=True)
    tree = il.sublink.loadTree(basePath, snap_z0, sid,
                                fields=['SubhaloGrNr', 'SnapNum', 'SubhaloID', 
                                        'FirstSubhaloInFOFGroupID'],
                                onlyMPB=True)

    if tree is None or len(tree['SnapNum']) == 0:
        snap_indep[i] = -2  # no tiene historia de fusión
        continue

    grps = tree['SubhaloGrNr']
    snaps = tree['SnapNum']
    subhalo_ids = tree['SubhaloID']
    first_in_fof = tree['FirstSubhaloInFOFGroupID']
    
    # Determinar si es central en cada snapshot: True si SubhaloID == FirstSubhaloInFOFGroupID
    is_central = (subhalo_ids == first_in_fof)
    
    # Buscar momentos donde el subhalo (o sus progenitores) fue central
    if np.any(is_central):
        # Encontrar el snapshot más reciente donde fue central
        central_snaps = snaps[is_central]
        max_snap = np.max(central_snaps)
        snap_indep[i] = max_snap
        
        # Encontrar el ID del halo FOF en ese momento
        idx_indep = np.where(snaps == max_snap)[0]
        if len(idx_indep) > 0:
            res_fof_indep[i] = grps[idx_indep[0]]
    else:
        # Nunca fue central en ninguna rama del árbol
        snap_indep[i] = -2  # Marca cuando nunca fue independiente
    
    processed += 1

print(f"\n Procesamiento completado:")
print(f"  - Subhalos procesados correctamente: {processed}")

# GUARDAR RESULTADOS COMPLETOS
output_file = 'subhalo_full_info_TNG50-1_disks.dat'
print(f"\nGuardando resultados completos en '{output_file}'...")

# Crear array con todos los datos
data = np.column_stack([
    target_ids,            # SubhaloID en z=0
    res_fof_z0,            # FOF group en z=0
    snap_indep,            # Snapshot de independencia (árbol completo)
    res_fof_indep,         # FOF group en el momento de independencia
    res_mass_z0,           # Masa del subhalo en z=0
    res_mass_fof_z0,       # Masa del FOF en z=0
    has_disk               # 1=Disco, 0=No Disco
])

np.savetxt(output_file,
           data,
           header='SubhaloID_z0   FOF_z0   Snap_indep   FOF_indep   Mass_sub_z0   Mass_FOF_z0   Has_Disk',
           fmt='%d %d %d %d %.6e %.6e %d')

print(f" Archivo '{output_file}' guardado correctamente.")

print("\n=== ESTADÍSTICAS ===")

# Clasificar subhalos según su historia
nunca_central = np.sum(snap_indep == -2)
fueron_centrales = np.sum(snap_indep >= 0)

print(f"Subhalos totales procesados: {n_targets}")
print(f"  - Fueron centrales en algún momento: {fueron_centrales} ({100*fueron_centrales/n_targets if n_targets>0 else 0:.1f}%)")
print(f"  - Nunca fueron centrales (siempre satélites): {nunca_central} ({100*nunca_central/n_targets if n_targets>0 else 0:.1f}%)")

print("\n Análisis completo finalizado.")