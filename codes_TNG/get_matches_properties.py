import illustris_python as il
import numpy as np
import pandas as pd
import sys

# ================= CONFIGURACIÓN =================
# Ruta a la simulación (AJUSTAR ESTA RUTA EN LA MÁQUINA DE DATOS)
BASE_PATH = "output"
SNAP_NUM = 99  # z=0

# Archivos de entrada y salida
INPUT_FILE = "subhalo_ids.txt"  # Archivo con una columna de IDs
OUTPUT_FILE = "subhalo_properties.csv"
# =================================================

def main():
    print(f"Leyendo IDs de: {INPUT_FILE}")
    print(f"Ruta TNG50-1: {BASE_PATH}")
    
    # Leer archivo asumiendo que no tiene encabezado y es una sola columna
    try:
        df_ids = pd.read_csv(INPUT_FILE, sep=r'\s+', header=None)
        ids = df_ids.iloc[:, 0].astype(int).values
        print(f"IDs a procesar: {len(ids)}")
    except Exception as e:
        print(f"Error leyendo archivo de entrada {INPUT_FILE}: {e}")
        sys.exit(1)

    print("Cargando catálogo completo de subhalos (esto puede tardar un poco)...")
    
    fields = ['SubhaloMassType', 'SubhaloSpin', 'SubhaloVelDisp', 'SubhaloVmax', 'SubhaloVmaxRad']
    
    try:
        # Cargar todos los subhalos de una vez
        subhalos = il.groupcat.loadSubhalos(BASE_PATH, SNAP_NUM, fields=fields)
        
        if not subhalos or ('count' in subhalos and subhalos['count'] == 0):
             print("Error: No se pudieron cargar los subhalos o el catálogo está vacío.")
             sys.exit(1)
             
        n_subhalos = subhalos['count']
        print(f"Catálogo cargado. Total subhalos en simulación: {n_subhalos}")
        
    except Exception as e:
        print(f"Error cargando catálogo de subhalos: {e}")
        sys.exit(1)

    results = []
    
    print("Extrayendo propiedades...")
    for i, subhalo_id in enumerate(ids):
        if subhalo_id >= n_subhalos or subhalo_id < 0:
            print(f"  Warning: ID {subhalo_id} fuera de rango (0-{n_subhalos-1})")
            continue
            
        try:
            # Extraer datos usando el ID como índice
            # SubhaloMassType es (N, 6)
            gas_mass = subhalos['SubhaloMassType'][subhalo_id][0]
            stellar_mass = subhalos['SubhaloMassType'][subhalo_id][4]
            
            # SubhaloSpin es (N, 3)
            spin = subhalos['SubhaloSpin'][subhalo_id]
            spin_mag = np.linalg.norm(spin)
            
            entry = {
                'SubhaloID': subhalo_id,
                'GasMass': gas_mass,
                'StellarMass': stellar_mass,
                'SubhaloSpin_x': spin[0],
                'SubhaloSpin_y': spin[1],
                'SubhaloSpin_z': spin[2],
                'SubhaloSpin_mag': spin_mag,
                'SubhaloVelDisp': subhalos['SubhaloVelDisp'][subhalo_id],
                'SubhaloVmax': subhalos['SubhaloVmax'][subhalo_id],
                'SubhaloVmaxRad': subhalos['SubhaloVmaxRad'][subhalo_id]
            }
            results.append(entry)
            
        except Exception as e:
            print(f"Error procesando ID {subhalo_id}: {e}")

    # Crear DataFrame final
    df_results = pd.DataFrame(results)
    
    # Guardar
    df_results.to_csv(OUTPUT_FILE, index=False)
    print(f"\nResultados guardados en {OUTPUT_FILE}")
    print(f"Total procesados exitosamente: {len(df_results)}")

if __name__ == "__main__":
    main()
