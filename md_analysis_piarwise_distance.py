#!/usr/bin/env python3
import warnings
warnings.filterwarnings('ignore')

import pandas as pd 
import MDAnalysis as mda


def plot_distance_mdanalysis(un, choice):

    """
    un is the md analysis universe created by reading psf and dcd files.
    choice:  residues names you want to calculate pairwise distance 
    
    """
    sel_choice = un.select_atoms(choice)
    self_dists=[]
    pair_dist_list = []
    for i, y in enumerate(zip(sel_choice.resnames, sel_choice.names, sel_choice.segids, sel_choice.resids)):
        for j in range(i+1, len(sel_choice.resids)):
            name_iden = y[2][2:]+ y[0] + str(y[3])+ "_" + sel_choice.segids[j][2:] + sel_choice.resnames[j]+ str(sel_choice.resids[j])
            pair_dist_list.append(name_iden)

    for ts in un.trajectory:
        time = un.trajectory.time
        self_dists.append(mda.lib.distances.self_distance_array(sel_choice.positions))


    df = pd.DataFrame(self_dists,columns= pair_dist_list)
    return df


if __name__ == "__main__":
    selection_res = '(name CA and resid 225 235 130 163 and segid NQ8) or (name CA and resid 220 and segid NQ8)\
                or (name CA and resid 251 233 and segid NQ8) or (name CA and resid 240 and segid NQ8) \
                or (name CA and resid 50 216 and segid NQ4) or (name CA and resid 139 and segid NQ4)\
                or (name CA and resid 74 and segid NQ7) or (name CA and resid 72 and segid NQ7) \
                or (name CA and resid 60 and segid NQ7) or (name CA and resid 5 and segid NQ10) \
                or (name CA and resid 32 and segid NQ11)'

    file_path = "/Users/umeshkhaniya/Dropbox (City College)/Umesh/paper_umesh_microstate_complexI/MD_file/MD_trajectory"

    u_apo = mda.Universe(file_path + "/apo/ComplexI_popc_ox-red_withWat.psf",
                  file_path + "/apo/ComplexI_popc_ox-red_eq6.dcd")

    plot_distance_mdanalysis(u_apo, selection_res).to_csv("apo_calpha_all_distance.csv")

    
        
    


