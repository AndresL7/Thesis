import numpy as np
import illustris_python as il

# Configuración
basePath = "../TNG50-1/output"
snapNum = 99
halos_file = "halos.txt"
output_file = "fof_disk_stats.txt"

print(f"Cargando IDs de subhalos con disco desde {halos_file}...")
subhalos_disc = np.loadtxt(halos_file, dtype=int)

print(f"Cargando catálogo de subhalos (SubhaloGrNr) del snapshot {snapNum}...")
# Solo necesitamos SubhaloGrNr para saber quién es el padre
sub_fields = ['SubhaloGrNr']
subs = il.groupcat.loadSubhalos(basePath, snapNum, fields=sub_fields)
sub_gr_nr = subs['SubhaloGrNr']

# Identificar a qué halo pertenece cada subhalo con disco
fof_ids_of_discs = sub_gr_nr[subhalos_disc]

# Contar cuántos discos hay en cada halo
# unique_fof_ids: IDs de los halos que tienen al menos un disco
# counts: Cuántos subhalos con disco tiene cada uno
unique_fof_ids, counts = np.unique(fof_ids_of_discs, return_counts=True)

print(f"Se encontraron {len(unique_fof_ids)} halos únicos que contienen los {len(subhalos_disc)} discos.")

print("Cargando masas de los halos...")
halo_fields = ['Group_M_TopHat200']
halos = il.groupcat.loadHalos(basePath, snapNum, fields=halo_fields)
all_masses = halos['Group_M_TopHat200']

# Obtener masas de los halos relevantes
relevant_masses = all_masses[unique_fof_ids]

print(f"Guardando estadísticas en {output_file}...")
# Columnas: HaloID, N_subhalos_con_disco, Masa_Halo
data = np.column_stack((unique_fof_ids, counts, relevant_masses))

header = "HaloID N_subhalos_with_disk Mass_1e10_Msun_h"
np.savetxt(output_file, data, fmt="%d %d %.6f", header=header)

print("Listo.")
