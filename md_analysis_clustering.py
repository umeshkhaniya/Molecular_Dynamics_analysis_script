import MDAnalysis as mda
import MDAnalysis.analysis.rms
import MDAnalysis.analysis.encore 
import warnings
warnings.filterwarnings('ignore')

select_res = '(backbone and (resid 225 235 130 163 251 233 240 220) and segid NQ8) \
                                           or (backbone and (resid 50 216 139) and segid NQ4) or\
                                           (backbone and (resid 72 74 60) and segid NQ7) or \
                                           (backbone and (resid 5) and segid NQ10) or \
                                                       (backbone and (resid 32) and segid NQ11)'
name_label = ["apo", "MQ", "Quinol"]


if __name__ == "__main__":
    file_path = "/Users/umeshkhaniya/Dropbox (City College)/Umesh/paper_umesh_microstate_complexI/MD_file/MD_trajectory"

    u_apo = mda.Universe(file_path + "/apo/ComplexI_popc_ox-red_withWat.psf",
                  file_path + "/apo/ComplexI_popc_ox-red_eq6.dcd")

    u_mq = mda.Universe(file_path + "/quinone/qdocked_red.psf",
                  file_path + "/quinone/SMD_eve10_qdocked-red.dcd")

    u_quinol= mda.Universe(file_path + "/quinol/CI-MQ8-OH.psf",
                  file_path + "/quinol/ComplexI_popc_eq23r-SMD-OH-st10.dcd")


    cluster_apo= MDAnalysis.analysis.encore.cluster(u_apo, selection= select_res)
    cluster_mq= MDAnalysis.analysis.encore.cluster(u_mq, selection= select_res)
    cluster_quinol= MDAnalysis.analysis.encore.cluster(u_quinol, selection= select_res)

    universes= [u_apo, u_mq, u_quinol]
    cluster_set = [cluster_apo, cluster_mq, cluster_quinol]


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
            












