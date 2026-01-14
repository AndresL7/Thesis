import numpy as np
import illustris_python as il
from scipy.spatial import cKDTree

# -------------------------
# Config
# -------------------------
basePath_1 = "../work/TNG50-1"
basePath_2 = "../work/TNG50-3/output"
Lbox = 35000.0   # ckpc/h
snapNum = 99

out_npy_candidates = "disco2.npy"
out_npy_rest = "no_disco2.npy"

# -------------------------
# Cargar catálogos
# -------------------------
print("Loading halos...")
fof1 = il.groupcat.loadHalos(basePath_1, snapNum)
fof2 = il.groupcat.loadHalos(basePath_2, snapNum)

pos1  = fof1['GroupPos'] % Lbox       # (N_halos1, 3)
mass1 = fof1['Group_M_TopHat200']       # unidades Illustris (1e10 Msun/h)
rvir1 = fof1['Group_R_TopHat200']       # ckpc/h
first_sub1 = fof1['GroupFirstSub']      # ID del subhalo central

pos2  = fof2['GroupPos'] % Lbox
mass2 = fof2['Group_M_TopHat200']
rvir2 = fof2['Group_R_TopHat200']

print("N halos:", len(mass1), len(mass2))

# -------------------------
# Filtrar halos válidos en TNG50-3 (masa > 0)
# -------------------------
valid = mass2 > 0
pos2_valid  = pos2[valid]
mass2_valid = mass2[valid]
rvir2_valid = rvir2[valid]
orig_idx2   = np.nonzero(valid)[0]   # índices originales en catálogo fof3

print(f"Halos válidos en TNG50-2: {len(pos2_valid)} / {len(pos2)}")

# KDTree con periodicidad (boxsize)
tree = cKDTree(pos2_valid, boxsize=Lbox * 1.0001)

# -------------------------
# Leer subhalos de TNG50-1 y construir lista de halos objetivo
# -------------------------
print("Loading subhalos...")
subs1 = il.groupcat.loadSubhalos(basePath_1, snapNum)
sub_to_fof = subs1['SubhaloGrNr']   # map subhalo -> halo FoF index

all_subhalos = np.arange(len(sub_to_fof))

subhalos_disc = np.loadtxt("../work/halos.txt", dtype=int)   # lista de subhalo IDs
# Halos (uno por subhalo en la lista).
halos_disc = sub_to_fof[subhalos_disc]

# Identificar halos únicos que contienen discos
unique_halos_with_disk = np.unique(halos_disc)

# Definir el "resto" como halos que NO contienen ningún disco
all_halos_ids = np.arange(len(pos1))
halos_rest = np.setdiff1d(all_halos_ids, unique_halos_with_disk)

# Para el grupo "rest", usamos el subhalo central como referencia
subhalos_rest = first_sub1[halos_rest]

print(f"Halos en halos.txt (len lista_subhalos): {len(subhalos_disc)}")
print(f"Halos FOF únicos con discos: {len(unique_halos_with_disk)}")
print(f"Halos restantes (sin discos): {len(halos_rest)}")

# -------------------------
# Función helper: distancia mínima periódica (vector)
# -------------------------
def min_image_delta(vec_a, vec_b, boxsize):
    # devuelve delta = vec_b - vec_a en mínima imagen
    delta = vec_b - vec_a
    # shift a la mínima imagen
    delta = (delta + 0.5 * boxsize) % boxsize - 0.5 * boxsize
    return delta

# -------------------------
# Función de matching para una lista de halos (hid_list)
# -------------------------
def match_halos(subhalos_list, hid_list):
    rows = []

    for i, hid in enumerate(hid_list):
        center = pos1[hid]
        Rvir   = rvir1[hid]
        Mvir   = mass1[hid]

        # Buscar candidatos dentro de 2*Rvir usando la KDTree periódica
        idxs_local = tree.query_ball_point(center, r=2.0 * Rvir)

        best = None
        best_cost = np.inf

        # si no hay candidatos, guardaremos campos None
        if len(idxs_local) == 0:
            best = None
        else:
            for idx_local in idxs_local:
                # idx_local indexa pos2_valid; mapear a índice original con orig_idx2
                orig_idx = int(orig_idx2[idx_local])
                pos2_cand = pos2_valid[idx_local]
                M2 = mass2_valid[idx_local]
                R2 = rvir2_valid[idx_local]

                # calcular delta periódico y distancia
                delta = min_image_delta(center, pos2_cand, Lbox)
                dist = np.linalg.norm(delta)

                # diferencias normalizadas
                dM = (M2 - Mvir) / Mvir if Mvir != 0 else np.inf

                # función de coste: raíz de suma de cuadrados (aquí usamos dist/Rvir y dM)
                cost = np.sqrt((dist / Rvir) ** 2 + dM ** 2)

                if cost < best_cost:
                    best_cost = cost
                    best = {
                        "halo_TNG50-2": orig_idx,
                        "Mvir_TNG50-2": float(M2),
                        "dist": float(dist),
                        "cost": float(cost),
                    }

        # Construir fila (llenar None si no hay match)
        if best is None:
            row = {
                "subhalo_TNG50-1": int(subhalos_list[i]),
                "halo_TNG50-1": int(hid),
                "Mvir_TNG50-1": float(Mvir),
                "halo_TNG50-2": None,
                "Mvir_TNG50-2": None,
                "dist": None,
                "cost": None,
            }
        else:
            row = {
                "subhalo_TNG50-1": int(subhalos_list[i]),
                "halo_TNG50-1": int(hid),
                "Mvir_TNG50-1": float(Mvir),
                "halo_TNG50-2": best["halo_TNG50-2"],

                **best
            }

        rows.append(row)

    return rows

# -------------------------
# Ejecutar matching para los dos conjuntos y guardar
# -------------------------
print("Matching halos in halos.txt (candidates)...")
matches_candidates = match_halos(subhalos_disc, halos_disc)
np.save(out_npy_candidates, matches_candidates, allow_pickle=True)

print("Matching halos_rest (rest)...")
matches_rest = match_halos(subhalos_rest, halos_rest)
np.save(out_npy_rest, matches_rest, allow_pickle=True)

# -------------------------
# Estadísticas finales
# -------------------------
masses_disc = mass1[halos_disc]
if len(masses_disc) > 0:
    min_mass_disc = np.min(masses_disc)
    max_mass_disc = np.max(masses_disc)
    
    masses_rest = mass1[halos_rest]
    count_rest_in_range = np.sum((masses_rest >= min_mass_disc) & (masses_rest <= max_mass_disc))
    
    print(f"Masa mínima de halos con disco: {min_mass_disc} (1e10 Msun/h)")
    print(f"Masa máxima de halos con disco: {max_mass_disc} (1e10 Msun/h)")
    print(f"Cantidad de halos sin disco en el rango [{min_mass_disc}, {max_mass_disc}]: {count_rest_in_range}")
else:
    print("No hay halos con disco para calcular rango de masas.")
