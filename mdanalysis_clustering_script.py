import MDAnalysis as mda
import MDAnalysis.analysis.rms
import MDAnalysis.analysis.encore 
import warnings
warnings.filterwarnings('ignore')
# define the name of residues based on you want do clucter
select_res = '(backbone and (resid 128 126 68 173) and chain H) \
                                           or (backbone and (resid 236 240) and chain M)'
name_label = ["brc"]


if __name__ == "__main__":
    file_path = "/Users/umeshkhaniya/Dropbox (City College)/2022_bin_traj"

    u_brc = mda.Universe(file_path + "/step5_input.psf",
                  [file_path + "/step7_50.dcd", file_path + "/step7_51.dcd", file_path + "/step7_52.dcd"])

    cluster_brc= MDAnalysis.analysis.encore.cluster(u_brc, selection= select_res)
    
    universes= [u_brc]
    cluster_set = [cluster_brc]


    for i, cluster_collection in enumerate(cluster_set):
        with open(f"cluster_info_{name_label[i]}", "w") as clust_info:
            clust_info.write('\n')
            clust_info.write(name_label[i]+':\n')
            for cluster in cluster_collection:
                clust_info.write(str(cluster.id) +' ' + '(size: '+ str(cluster.size)+', centroid: '+ str(cluster.centroid)\
                           + ') elements: ' + str(cluster.elements)+'\n')
                time_point=cluster.centroid
                universes[i].trajectory[time_point]
                with mda.Writer(name_label[i]+'_'+str(time_point)+'_frame.pdb') as pdb:
                    pdb.write(universes[i])
            












