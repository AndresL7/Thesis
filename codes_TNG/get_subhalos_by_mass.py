import numpy as np
import illustris_python as il

# ============ CONFIGURACIÓN ============
SIMULATION = "TNG50-1" 
BASE_PATH = "../"
SNAP_Z0 = 99    # z=0
# =======================================

# Rango de masas subhalos en unidades de 10^10 Msun/h
MASS_MIN_sub = 1.466765
MASS_MAX_sub = 1013.933

# Rango de masas halos FOF en unidades de 10^10 Msun/h
MASS_MIN_fof = 107.67
MASS_MAX_fof = 15167.14

def get_subhalos_in_mass_range_with_fof_filter(basepath, snapnum, mass_min_sub, mass_max_sub, 
                                                mass_min_fof, mass_max_fof):
    """
    Obtiene los IDs de todos los subhalos cuya masa está en el rango especificado
    Y que pertenecen a halos FOF cuya masa también está en el rango especificado.
    
    Parameters:
    -----------
    basepath : str
        Ruta base de la simulación
    snapnum : int
        Número del snapshot (99 para z=0)
    mass_min_sub : float
        Masa mínima de subhalos en unidades de 10^10 Msun/h
    mass_max_sub : float
        Masa máxima de subhalos en unidades de 10^10 Msun/h
    mass_min_fof : float
        Masa mínima de halos FOF en unidades de 10^10 Msun/h
    mass_max_fof : float
        Masa máxima de halos FOF en unidades de 10^10 Msun/h
    
    Returns:
    --------
    subhalo_ids : array
        Array con los IDs de los subhalos que cumplen ambos criterios
    subhalo_masses : array
        Array con las masas de los subhalos
    fof_masses : array
        Array con las masas de los halos FOF correspondientes
    """
    print(f"Cargando datos de subhalos y halos FOF del snapshot {snapnum}...")
    
    # Cargar masas de subhalos y su grupo FOF padre
    subhalos = il.groupcat.loadSubhalos(basepath, snapnum, 
                                        fields=['SubhaloMass', 'SubhaloGrNr'])
    subhalo_masses = subhalos['SubhaloMass']
    subhalo_grnr = subhalos['SubhaloGrNr']  # ID del halo FOF al que pertenece
    
    # Cargar masas de halos FOF
    halos = il.groupcat.loadHalos(basepath, snapnum, fields=['GroupMass'])
    fof_masses = halos['GroupMass']
    
    print(f"Total de subhalos: {len(subhalo_masses)}")
    print(f"Total de halos FOF: {len(fof_masses)}")
    
    print(f"\nFiltros aplicados:")
    print(f"  Subhalos - Masa mínima: {mass_min_sub:.4e} (10^10 Msun/h)")
    print(f"  Subhalos - Masa máxima: {mass_max_sub:.4e} (10^10 Msun/h)")
    print(f"  Halos FOF - Masa mínima: {mass_min_fof:.4e} (10^10 Msun/h)")
    print(f"  Halos FOF - Masa máxima: {mass_max_fof:.4e} (10^10 Msun/h)")
    
    # Filtrar subhalos por masa
    mask_subhalo = (subhalo_masses >= mass_min_sub) & (subhalo_masses <= mass_max_sub)
    
    # Filtrar halos FOF por masa
    mask_fof = (fof_masses >= mass_min_fof) & (fof_masses <= mass_max_fof)
    fof_ids_valid = np.where(mask_fof)[0]
    
    # Combinar filtros: subhalo debe estar en rango de masas Y pertenecer a un halo FOF válido
    mask_combined = mask_subhalo & np.isin(subhalo_grnr, fof_ids_valid)
    
    subhalo_ids = np.where(mask_combined)[0]
    subhalo_masses_filtered = subhalo_masses[mask_combined]
    fof_masses_filtered = fof_masses[subhalo_grnr[mask_combined]]
    
    print(f"\n Subhalos en rango de masas (antes de filtro FOF): {np.sum(mask_subhalo)}")
    print(f" Halos FOF en rango de masas: {len(fof_ids_valid)}")
    print(f" Subhalos que cumplen AMBOS criterios: {len(subhalo_ids)}")
    
    return subhalo_ids, subhalo_masses_filtered, fof_masses_filtered


def main():
    # Path de la simulación
    sim_path = f"{BASE_PATH}{SIMULATION}/output"
    
    # Obtener subhalos que cumplen ambos criterios
    subhalo_ids, subhalo_masses, fof_masses = get_subhalos_in_mass_range_with_fof_filter(
        sim_path, SNAP_Z0, 
        MASS_MIN_sub, MASS_MAX_sub,
        MASS_MIN_fof, MASS_MAX_fof
    )
    
    # Mostrar estadísticas
    print(f"\n{'='*60}")
    print(f"RESULTADOS")
    print(f"{'='*60}")
    print(f"Simulación: {SIMULATION}")
    print(f"Snapshot: {SNAP_Z0} (z=0)")
    print(f"\nSubhalos encontrados: {len(subhalo_ids)}")
    
    if len(subhalo_ids) > 0:
               
        # Mostrar primeros 10 IDs como muestra
        print(f"\nPrimeros 10 IDs de subhalos (muestra):")
        print(f"  {'#':<4} {'SubID':<8} {'Masa Sub':<15} {'Masa FOF':<15}")
        for i, (sid, smass, fmass) in enumerate(zip(subhalo_ids[:10], subhalo_masses[:10], fof_masses[:10])):
            print(f"  {i+1:<4} {sid:<8} {smass:.4e}     {fmass:.4e}")
        
        if len(subhalo_ids) > 10:
            print(f"  ... ({len(subhalo_ids) - 10} subhalos más)")
        
        # Guardar solo los IDs en archivo
        output_file = 'subhalo_ids.txt'
        np.savetxt(output_file, subhalo_ids, fmt='%d')
        print(f"\n✓ IDs guardados en '{output_file}'")
    else:
        print("\nNo se encontraron subhalos que cumplan ambos criterios.")
    
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
