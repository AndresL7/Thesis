import numpy as np
import illustris_python as il
import sys

# ============ CONFIGURACIÓN ============
SIMULATION = "TNG50-1" 
BASE_PATH = "../"
SNAP_Z0 = 99    # z=0
# =======================================

def get_gas_particles_from_subhalo(basePath, subhalo_id, snap=99):
    """
    Obtiene los IDs y posiciones de las partículas de gas de un subhalo específico.
    
    Parameters:
    -----------
    basePath : str
        Ruta base de la simulación
    subhalo_id : int
        ID del subhalo
    snap : int
        Número del snapshot (default: 99 para z=0)
        
    Returns:
    --------
    gas_ids : array
        IDs de las partículas de gas
    gas_positions : array
        Posiciones [x, y, z] de las partículas de gas en ckpc/h
    """
    
    print(f"Cargando partículas de gas del subhalo {subhalo_id} en snapshot {snap}...")
    
    try:
        # Cargar IDs y coordenadas de partículas de gas (PartType0)
        cutout = il.snapshot.loadSubhalo(basePath, snap, subhalo_id, 'gas', 
                                         fields=['ParticleIDs', 'Coordinates'])
        
        gas_ids = cutout['ParticleIDs']
        gas_positions = cutout['Coordinates']
        
        print(f"  ✓ {len(gas_ids)} partículas de gas encontradas")
        
        return gas_ids, gas_positions
        
    except Exception as e:
        print(f"  ✗ Error al cargar partículas de gas: {e}")
        return None, None


def save_gas_particles(subhalo_id, gas_ids, gas_positions, filename=None):
    """
    Guarda los IDs y posiciones de las partículas de gas en un archivo.
    
    Parameters:
    -----------
    subhalo_id : int
        ID del subhalo
    gas_ids : array
        IDs de las partículas de gas
    gas_positions : array
        Posiciones de las partículas de gas
    filename : str, optional
        Nombre del archivo de salida. Si no se proporciona, se genera automáticamente.
    """
    
    if filename is None:
        filename = f'gas_particles_subhalo_{subhalo_id}.dat'
    
    # Combinar IDs y posiciones en una sola matriz
    data = np.column_stack([gas_ids, gas_positions])
    
    # Guardar en archivo
    np.savetxt(filename, data,
               header=f'Subhalo {subhalo_id} - Partículas de gas (PartType0)\nParticleID   X[ckpc/h]   Y[ckpc/h]   Z[ckpc/h]',
               fmt='%d %.6f %.6f %.6f')
    
    print(f"\n  Archivo guardado: '{filename}'")
    print(f"  Total de partículas: {len(gas_ids)}")


def main():
    # Verificar si se proporciona el ID del subhalo como argumento
    if len(sys.argv) > 1:
        subhalo_id = int(sys.argv[1])
    else:
        # Pedir el subhalo ID al usuario
        subhalo_id = int(input("Ingresa el ID del subhalo: "))
    
    # Path de la simulación
    sim_path = f"{BASE_PATH}{SIMULATION}/output"
    
    print(f"\n{'='*60}")
    print(f"EXTRAYENDO PARTÍCULAS DE GAS")
    print(f"Simulación: {SIMULATION}")
    print(f"Snapshot: {SNAP_Z0} (z=0)")
    print(f"Subhalo ID: {subhalo_id}")
    print(f"{'='*60}\n")
    
    # Obtener partículas de gas
    gas_ids, gas_positions = get_gas_particles_from_subhalo(sim_path, subhalo_id, SNAP_Z0)
    
    if gas_ids is not None and len(gas_ids) > 0:
        # Mostrar estadísticas básicas
        print(f"\nEstadísticas de posición:")
        print(f"  X: [{gas_positions[:,0].min():.2f}, {gas_positions[:,0].max():.2f}] ckpc/h")
        print(f"  Y: [{gas_positions[:,1].min():.2f}, {gas_positions[:,1].max():.2f}] ckpc/h")
        print(f"  Z: [{gas_positions[:,2].min():.2f}, {gas_positions[:,2].max():.2f}] ckpc/h")
        print(f"  Centro de masa: [{gas_positions[:,0].mean():.2f}, {gas_positions[:,1].mean():.2f}, {gas_positions[:,2].mean():.2f}] ckpc/h")
        
        # Guardar en archivo
        save_gas_particles(subhalo_id, gas_ids, gas_positions)
        
        print(f"\n{'='*60}")
        print("COMPLETADO")
        print(f"{'='*60}")
    else:
        print("\n  ✗ No se encontraron partículas de gas o el subhalo no existe")


if __name__ == "__main__":
    main()
