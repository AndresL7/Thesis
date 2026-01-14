import numpy as np
import illustris_python as il
import pandas as pd
import os

"""
Encuentra todos los halos FOF en un rango de masas que contienen subhalos con discos,
y guarda un archivo indicando si cada halo FOF en ese rango tiene o no un subhalo con disco.
"""


# ============ CONFIGURATION ============
# Adjust these paths if necessary to match your directory structure
BASE_PATH = "../"
SIMULATION = "TNG50-1"
SNAP_Z0 = 99
HALOS_FILE = "halos.txt"
OUTPUT_FILE = "fof_halos_mass_range_status.txt"
# =======================================

def main():
    # Path to simulation output
    basePath = os.path.join(BASE_PATH, SIMULATION, "output")
    
    if not os.path.exists(basePath):
        print(f"Error: Simulation path not found: {basePath}")
        return

    # 1. Load subhalo IDs with disks
    print(f"Loading {HALOS_FILE}...")
    if not os.path.exists(HALOS_FILE):
        print(f"Error: Could not find {HALOS_FILE} in the current directory.")
        return
        
    try:
        disk_subhalo_ids = np.loadtxt(HALOS_FILE, dtype=int)
        # Handle case where file has only one entry or is empty
        if disk_subhalo_ids.size == 0:
            print("Error: halos.txt is empty.")
            return
        if disk_subhalo_ids.ndim == 0:
            disk_subhalo_ids = np.array([disk_subhalo_ids])
    except Exception as e:
        print(f"Error reading {HALOS_FILE}: {e}")
        return

    print(f"Found {len(disk_subhalo_ids)} disk subhalos.")

    # 2. Load SubhaloGrNr to find which FoF halo each subhalo belongs to
    print("Loading SubhaloGrNr from GroupCat...")
    # subhalo_grnr[i] is the FoF ID of subhalo i
    try:
        subhalo_grnr = il.groupcat.loadSubhalos(basePath, SNAP_Z0, fields=['SubhaloGrNr'])
    except Exception as e:
        print(f"Error loading Subhalo catalog: {e}")
        return

    # 3. Identify FoF halos containing the disk subhalos
    # Check if indices are within bounds
    if np.max(disk_subhalo_ids) >= len(subhalo_grnr):
        print("Error: Some subhalo IDs in halos.txt are out of bounds for the catalog.")
        return

    fof_ids_with_disks = subhalo_grnr[disk_subhalo_ids]
    unique_fof_ids_with_disks = np.unique(fof_ids_with_disks)
    
    print(f"Identified {len(unique_fof_ids_with_disks)} unique FoF halos containing disk subhalos.")

    # 4. Load Halo masses
    print("Loading Halo masses (Group_M_TopHat200)...")
    # using Group_M_TopHat200 as mass definition. Units: 10^10 Msun/h
    try:
        halo_masses = il.groupcat.loadHalos(basePath, SNAP_Z0, fields=['Group_M_TopHat200'])
    except Exception as e:
        print(f"Error loading Halo catalog: {e}")
        return
    
    # 5. Get mass range of the identified FoF halos
    masses_of_interest = halo_masses[unique_fof_ids_with_disks]
    min_mass = np.min(masses_of_interest)
    max_mass = np.max(masses_of_interest)
    
    print(f"Mass range of FoF halos with disks: {min_mass:.6e} to {max_mass:.6e} (10^10 Msun/h)")

    # 6. Find all FoF halos in this mass range
    # Indices of halos within the mass range
    halos_in_range_indices = np.where((halo_masses >= min_mass) & (halo_masses <= max_mass))[0]
    
    print(f"Found {len(halos_in_range_indices)} total FoF halos in this mass range.")

    # 7. Create result columns
    # Column 1: FoF Halo Index (halos_in_range_indices)
    # Column 2: Boolean (True if in unique_fof_ids_with_disks)
    
    # Efficient check using set
    fof_with_disks_set = set(unique_fof_ids_with_disks)
    
    has_disk_column = [idx in fof_with_disks_set for idx in halos_in_range_indices]
    
    # Create DataFrame
    df = pd.DataFrame({
        'FoF_ID': halos_in_range_indices,
        'Has_Disk_Subhalo': has_disk_column,
        'Group_M_TopHat200': halo_masses[halos_in_range_indices]
    })

    # 8. Save to file
    print(f"Saving results to {OUTPUT_FILE}...")
    # Saving as tab-separated values
    df.to_csv(OUTPUT_FILE, index=False, sep='\t')
    print("Done.")

if __name__ == "__main__":
    main()
