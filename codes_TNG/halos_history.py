import illustris_python as il
import numpy as np

basepath = "output/"
snap_z0 = 99 # Snapshot number for z=0

#[0,0.1,0.24,0.5,0.7,1,1.25,1.5,2,2.44,3.01,4.01]
relevant_snaps = [99,91,81,67,59,50,44,40,33,29,25,21]

#load Halos
groups = il.groupcat.loadHalos(basepath, snap_z0, fields=['GroupFirstSub'])

#I get this number from the len of the results in Mdom in z=0
Nhalos = int(input("Number of halos:"))

with open("halos_history.txt","w") as f:
        for H in np.arange(0,Nhalos):
                if (H%1000==0):
                        print(f"Writing history from Halo {H}")
                histID = []
                snaps = []

                #obtain central subhalo id from halo with id H
                central_sub_id = groups[H]

                # Main tree branch of this subhalo
                tree = il.sublink.loadTree(basepath, snap_z0, central_sub_id,
                                           fields=['SnapNum', 'SubhaloGrNr'],
                                           onlyMPB=True)
                
                if tree is not None:
                        for snap, gid in zip(tree['SnapNum'], tree['SubhaloGrNr']):
                                if snap in relevant_snaps:
                                        histID.append(int(snap))
                                        snaps.append(int(gid))

                f.write(f"{snaps}, {histID}\n")
                    
