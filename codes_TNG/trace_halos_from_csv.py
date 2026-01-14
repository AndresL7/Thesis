import illustris_python as il
import numpy as np
import os

# --- CONFIGURACIÓN ---
basepath = "output/"
input_file = "halos.txt"  # Archivo de entrada: IDs de subhalos
output_file = "subhalos_host_history.txt"
snap_z0 = 99 # Snapshot de referencia para los IDs de entrada (z=0)

# Snapshots de interés para rastrear
relevant_snaps = [99, 91, 81, 67, 59, 50, 44, 40, 33, 29, 25, 21]
# ---------------------

def run_tracking():
    # Verificar archivo de entrada
    if not os.path.exists(input_file):
        print(f"Error: No se encuentra el archivo '{input_file}'.")
        return

    print(f"Leyendo lista de subhalos desde {input_file}...")
    
    # Cargar datos
    # Asumimos archivo de texto con IDs, sin cabecera, comentarios con #
    try:
        data = np.loadtxt(input_file, dtype=int)
    except:
        # Intento alternativo con coma si falla el espacio
        data = np.loadtxt(input_file, delimiter=',', dtype=int)

    # Asegurar que sea array 1D
    data = np.atleast_1d(data)
    
    # Si tiene múltiples columnas (shape N, M), tomamos la primera columna
    if data.ndim > 1:
        data = data[:, 0]

    print(f"Se procesarán {len(data)} subhalos (IDs correspondientes al Snap {snap_z0}).")

    with open(output_file, "w") as f:
        # Escribir cabecera solicitada
        f.write("subhaloID_z0,snap,HaloFoF_host\n")

        for i, sub_id in enumerate(data):
            if i % 100 == 0:
                print(f"Procesando {i}/{len(data)}... (Subhalo {sub_id})")

            # Cargar el árbol (Main Progenitor Branch) para este subhalo desde z=0
            try:
                tree = il.sublink.loadTree(basepath, snap_z0, sub_id,
                                           fields=['SnapNum', 'SubhaloGrNr'],
                                           onlyMPB=True)
            except Exception as e:
                print(f"Error cargando árbol para subhalo {sub_id}: {e}")
                continue

            if tree is not None:
                # Recorrer la historia del subhalo
                for snap, host_halo in zip(tree['SnapNum'], tree['SubhaloGrNr']):
                    # Guardar solo si es un snapshot relevante
                    if snap in relevant_snaps:
                        f.write(f"{sub_id},{snap},{host_halo}\n")

    print(f"Listo. Archivo guardado en: {output_file}")

if __name__ == "__main__":
    run_tracking()
