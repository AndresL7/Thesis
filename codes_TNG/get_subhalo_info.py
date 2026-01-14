import illustris_python as il
import numpy as np
import sys

# CONFIGURACIÓN
basePath = '../TNG50-1/output'
snap = 99  # snapshot de z=0 para TNG50-1

def get_subhalo_info(subhalo_id):
    """
    Obtiene el número de partículas y la masa de un subhalo dado su ID.
    
    Parámetros:
    -----------
    subhalo_id : int
        ID del subhalo en el snapshot especificado
        
    Retorna:
    --------
    dict con las siguientes claves:
        - 'SubhaloID': ID del subhalo
        - 'SubhaloMass': Masa total del subhalo (10^10 Msun/h)
        - 'SubhaloLen': Número total de partículas
        - 'SubhaloLenType': Número de partículas por tipo [gas, DM, ..., stars, BH]
    """
    
    print(f"\nCargando información del subhalo {subhalo_id}...")
    print(f"Ruta: {basePath}, Snapshot: {snap}")
    
    try:
        # Cargar datos del subhalo específico
        subhalo = il.groupcat.loadSingle(basePath, snap, subhaloID=subhalo_id)
        
        # Extraer información relevante
        info = {
            'SubhaloID': subhalo_id,
            'SubhaloMass': subhalo['SubhaloMass'],
            'SubhaloLen': subhalo['SubhaloLen'],
            'SubhaloLenType': subhalo['SubhaloLenType']
        }
        
        return info
        
    except Exception as e:
        print(f"Error al cargar el subhalo {subhalo_id}: {e}")
        return None


def print_subhalo_info(info):
    """Imprime la información del subhalo de forma legible."""
    
    if info is None:
        print("No se pudo obtener información del subhalo.")
        return
    
    print("\n" + "="*60)
    print(f"INFORMACIÓN DEL SUBHALO {info['SubhaloID']}")
    print("="*60)
    print(f"\nMasa total: {info['SubhaloMass']:.6e} (10^10 Msun/h)")
    print(f"             = {info['SubhaloMass'] * 1e10:.6e} Msun/h")
    print(f"\nNúmero total de partículas: {info['SubhaloLen']}")
    print(f"\nNúmero de partículas por tipo:")
    print(f"  - Gas (tipo 0):           {info['SubhaloLenType'][0]}")
    print(f"  - Materia oscura (tipo 1): {info['SubhaloLenType'][1]}")
    print(f"  - (tipo 2):               {info['SubhaloLenType'][2]}")
    print(f"  - (tipo 3):               {info['SubhaloLenType'][3]}")
    print(f"  - Estrellas (tipo 4):     {info['SubhaloLenType'][4]}")
    print(f"  - Agujeros negros (tipo 5): {info['SubhaloLenType'][5]}")
    print("="*60 + "\n")


if __name__ == "__main__":
    
    # Verificar si se proporcionó un ID como argumento
    if len(sys.argv) > 1:
        try:
            subhalo_id = int(sys.argv[1])
        except ValueError:
            print("Error: El ID del subhalo debe ser un número entero.")
            print("Uso: python get_subhalo_info.py <subhalo_id>")
            sys.exit(1)
    else:
        # ID por defecto para pruebas
        subhalo_id = 0
        print("No se proporcionó ID. Usando subhalo 0 como ejemplo.")
    
    # Obtener y mostrar información
    info = get_subhalo_info(subhalo_id)
    print_subhalo_info(info)
