import numpy as np
import h5py

def get_halfmass_radius(halo_ids, snapshot_path):
    """
    Obtiene el SubhaloHalfmassRad para una lista de subhalos.
    
    Parameters:
    -----------
    halo_ids : list or array
        Lista de IDs de subhalos
    snapshot_path : str
        Ruta al archivo HDF5 del snapshot de Illustris
    
    Returns:
    --------
    dict : Diccionario con {subhalo_id: halfmass_radius}
    """
    results = {}
    
    with h5py.File(snapshot_path, 'r') as f:
        # Leer SubhaloHalfmassRad del snapshot
        halfmass_rad = f['Subhalos']['SubhaloHalfmassRad'][:]
        
        for halo_id in halo_ids:
            if halo_id < len(halfmass_rad):
                results[halo_id] = halfmass_rad[halo_id]
            else:
                results[halo_id] = None
                print(f"Warning: Subhalo ID {halo_id} out of range")
    
    return results


def main():
    # Leer los IDs de subhalos desde halos.txt
    halo_ids = np.loadtxt('halos.txt', dtype=int)
    
    # Ruta al snapshot (ajustar según tu ubicación)
    snapshot_path = 'snapdir_099/snap_099.0.hdf5'  # Modificar según necesites
    
    # Obtener los radios de semimasa
    halfmass_radii = get_halfmass_radius(halo_ids, snapshot_path)
    
    # Guardar resultados en un archivo
    with open('halfmass_radius_results.txt', 'w') as f:
        f.write("SubhaloID SubhaloHalfmassRad\n")
        for halo_id, radius in halfmass_radii.items():
            f.write(f"{halo_id} {radius}\n")
    
    # Imprimir resultados
    print("SubhaloID | SubhaloHalfmassRad")
    print("-" * 35)
    for halo_id, radius in halfmass_radii.items():
        print(f"{halo_id:8d} | {radius}")
    
    print(f"\nResultados guardados en 'halfmass_radius_results.txt'")


if __name__ == "__main__":
    main()
