#!/usr/bin/env python

from scipy.stats import norm, skewnorm
import sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import math
from matplotlib.colors import LinearSegmentedColormap,ListedColormap,BoundaryNorm


import logging
import ms_analysis as msa 
import weighted_correlation as wc



logging.basicConfig(
    filename='summary_data.log',
    filemode = 'w',
    level=logging.INFO,
    format='%(asctime)s.%(msecs)03d %(levelname)s: %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
)


def histogram_ms_energy(ms_enrg_file):
    energy_lst_count =np.asarray([a for a,f in zip([x[0] for x in ms_enrg_file], [x[1] for x in ms_enrg_file]) for _ in range(f)])
    (mu, sigma) = norm.fit(energy_lst_count)
    skewness, mean, std = skewnorm.fit(energy_lst_count)
    fig = plt.figure(figsize = (10,8))
    graph_hist = plt.hist(energy_lst_count,bins=100, alpha=0.6)
    Y = graph_hist[0]
    y = skewnorm.pdf(np.array(energy_lst_count), skewness, mean, std)
    pdf_data = Y.max()/max(y)*y
    plt.plot(energy_lst_count,pdf_data , label = 'approximated skewnorm', color = 'black')
    plt.title("Skewness:"+str(round(skewness,2))+ " " + "Mean:"+str(round(mean,2))+ " " + "STD:"+str(round(std,2)), fontsize = 15)
    plt.yticks(fontsize = 15)
    plt.xticks(fontsize = 15)
    plt.xlabel(" Microstate Energy (Kcal/Mol)", fontsize = 15)
    plt.ylabel("Count", fontsize = 15)
    plt.tick_params(axis="x", direction="out", length=8, width=2)
    plt.tick_params(axis="y", direction="out", length=8, width=2)
    fig.savefig("enthalpy_dis.pdf", dpi = 300, bbox_inches = "tight")
    logging.info("Histogram figure is saved as name: enthalpy_dis.pdf")

    return 


    

def free_residue():
	free_residues = []
	for res in mc.free_residues:
	    free_residues.append(head3_conf[res[0]].resid)
	ms_free_residues = pd.DataFrame(free_residues,columns = ["Residue"])
	return ms_free_residues


def fixed_residue():
    fixed_res_crg_dict = {}
    for conf in head3_conf:
        if conf.iconf in mc.fixed_iconfs:
            if conf.resid not in fixed_res_crg_dict:
                fixed_res_crg_dict[conf.resid] = conf.crg
    background_charge = sum(fixed_res_crg_dict.values())
    fixed_interested_res = {i:j for i,j in fixed_res_crg_dict.items() if i[:3] in your_interested_res}
    fixed_residues_crg = pd.DataFrame(fixed_interested_res.items(), columns=['Residue', 'crg'])
    return background_charge, fixed_residues_crg



def id_crg_relation():
	id_vs_charge = {}
	for conf in head3_conf:
	    id_vs_charge[conf.iconf] = conf.crg
	return id_vs_charge 
    

def convert_ms_crg(l, d):

    crg_lst =[[y[0], y[1], [convert_ms_crg(x, d) if isinstance(x, list) else d.get(x, x) for x in y[2]]] for y in l]
    return crg_lst



def findUniqueCrgmsCountOrder(crg_list_ms, begin_energy = None, end_energy = None):
    """You need to supply the charge list microstate here. Make sure
    the crg file is sorted with increasing order.If you want to filter 
    charge id based on the energy. Then you need to supply here.
    unique_crg_state_order gives the order of unique charge state based on energy. 
    Lowest energy charge state will give the order 1 and then second unique charge 
    state will give the order 2. This order is based on unique charge ms order.
    """
    if not begin_energy and not end_energy:
        logging.info("All energy microstates are selected.")
        begin_energy = crg_list_ms[0][0]
        end_energy = crg_list_ms[-1][0]
    elif begin_energy and end_energy:
        crg_list_ms = [[x[0], x[1], x[2]] for x in crg_list_ms if x[0] >= begin_energy and x[0] <= end_energy]   
    else:
        logging.critical("No energy bound is defined")
        sys.exit('Give the lower or upper energy bound.')

        
    #unique charge as key and energy, count and order
    crg_all_count = {}
    unique_crg_state_order = 1
    for x, array in enumerate(crg_list_ms):
        if tuple(array[2]) not in crg_all_count.keys():
            crg_all_count[(tuple(array[2]))] = [array[1], [array[0]], [unique_crg_state_order]]
            unique_crg_state_order +=1
        else:
            crg_all_count[(tuple(array[2]))][0] += array[1]
            # add the maximum and minimum energy 
            min_energy = min(min(crg_all_count[(tuple(array[2]))][1]),array[0]) 
            max_energy = max(max(crg_all_count[(tuple(array[2]))][1]),array[0])
            # clear energy list and append minimum and maximum energy
            crg_all_count[(tuple(array[2]))][1].clear()
            crg_all_count[(tuple(array[2]))][1].append(min_energy)
            crg_all_count[(tuple(array[2]))][1].append(max_energy)

    # make a list of count, unique charge microstate, energy difference and order.
    all_crg_ms_unique  = []
    all_count = []
    energy_diff_all = []
    unique_crg_state_order = []
    for u,v in crg_all_count.items():
        all_crg_ms_unique.append(list(u))
        all_count.append(v[0])
        unique_crg_state_order.append(v[2][0])
        if len(v[1]) == 2:
            energy_diff_all.append(round(v[1][1]-v[1][0], 6))
        elif len(v[1]) == 1:
            energy_diff_all.append(0)
        else:
            logging.critical("There is error while creating unique charge state.")
            sys.exit("There is error while creating unique charge state.")
    logging.info(f"Total number of charge ms: {len(crg_list_ms)}")
    logging.info(f"Total number of  Unique charge ms: {len(all_crg_ms_unique)}")
    return all_crg_ms_unique, all_count, unique_crg_state_order, energy_diff_all 


def ConcaCrgMsPandas(unique_crg_ms_list, ms_count, ms_order, free_residues, background_charge, residue_interest_list):
    unique_crg_ms_list_pd = pd.DataFrame(unique_crg_ms_list).T
    ms_count_pd = pd.DataFrame(ms_count,columns = ["Count"]).T
    ms_order_pd = pd.DataFrame(ms_order,columns = ["Order"]).T
    crg_ms_count_pd = pd.concat([unique_crg_ms_list_pd, ms_count_pd, ms_order_pd])
    crg_count_res_1 = pd.concat([free_residues,crg_ms_count_pd], axis=1)
    crg_count_res_1.loc["Count", 'Residue'] = 'Count'
    crg_count_res_1.loc["Order", 'Residue'] = 'Order'
    all_crg_count_res = crg_count_res_1.set_index("Residue")
    # sort based on the count
    all_crg_count_res = all_crg_count_res.sort_values(by = "Count", axis = 1, ascending = False)
    all_crg_count_res.columns = range(all_crg_count_res.shape[1])
    all_crg_count_res = all_crg_count_res.T.set_index("Order")
    all_crg_count_res["Occupancy"] = round(all_crg_count_res["Count"]/ sum(all_crg_count_res["Count"]),3)
    all_crg_count_res['Sum_crg_protein'] =  all_crg_count_res.iloc[:, :-2].sum(axis = 1) + background_charge
    crg_count_res = all_crg_count_res.copy()
    for i in all_crg_count_res.columns:
        if i[:3] not in residue_interest_list and i != "Occupancy" and i != "Count" and i != "Sum_crg_protein":
            crg_count_res.drop([i], axis = 1, inplace = True)

    return crg_count_res


#which tautomer charge state is most populated. This includes the background charge also.
def plots_unique_crg_histogram(charge_ms_file, background_charge):
    x_av = [sum(x) + background_charge for x in charge_ms_file[0]]
    y_av = [math.log10(x) for x in charge_ms_file[1]]
    energy_diff_all_fl = [float(x) for x in charge_ms_file[3]]
    g1 = sns.JointGrid(marginal_ticks=True, height = 6)
    ax = sns.scatterplot(x=x_av, y=y_av, hue = energy_diff_all_fl, palette = 'viridis', size = energy_diff_all_fl, sizes=(10, 200), ax=g1.ax_joint)
    ax.set_xticks(range(int(min(x_av)), int(max(x_av)) + 1))
    ax.set_xlabel("Charge",fontsize=15)
    ax.set_ylabel("log$_{10}$(Count)",fontsize=16)
    ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    ax2 = sns.histplot(x=x_av, linewidth=2,discrete=True, ax=g1.ax_marg_x)

    ax2.set_ylabel(None,fontsize=16)
    g1.ax_marg_y.set_axis_off()
    g1.fig.subplots_adjust(top= 0.9)
    g1.fig.suptitle('All microstate energy', fontsize = 16)
    g1.savefig("all_en_cr_vs_log(count).pdf", dpi = 300, bbox_inches = "tight")

    logging.info("Unique charge distribuion is saved as: all_en_cr_vs_log(count).pdf")

    return 

def combine_free_fixed_residue(fixed_residue_file, all_crg_count_read):

    df_fixed_residues_crg_tar = fixed_residue_file.T
    df_fixed_residues_crg_tar.columns = df_fixed_residues_crg_tar.iloc[0]
    df_fixed_residues_crg_tar = df_fixed_residues_crg_tar.iloc[1:,:].reset_index(drop = True)
    df_fixed_res_dup = pd.concat([df_fixed_residues_crg_tar]*len(all_crg_count_read), ignore_index=True)
    df_join_free_back = all_crg_count_read.join(df_fixed_res_dup)
    return df_join_free_back



def correlation_data_parsing(df):
    all_crg_count =  df.iloc[:,:-2]
    all_crg_count = all_crg_count.T
    all_crg_count["std"] = all_crg_count.std(axis = 1).round(3)
    all_crg_count_std = all_crg_count.loc[all_crg_count['std'] != 0].T[:-1].reset_index(drop = True)
    logging.info(f"Number of residues that change the protonation state: {len(all_crg_count_std.columns)-1}")

    return all_crg_count_std
    


#rename the residues with shorter name and 
# keep acid first, then polar, base residues and then MQ in column
def renameOrderResidues(file_input):
    rename_dict = {}
    acid_list = []
    base_list = []
    polar_rest_list = []
    ub_q_list = []
    non_residue_list = []

    for i in file_input.columns[:-1]:
        rename_dict[i] = i[3] + str("_")+ i[:3] + str(int(i[4:8])) + i[8:]
    rename_dict["Count"] = 'Count'
        
    for x, y in rename_dict.items():
        if y[2:5] == 'ASP':
            rename_dict[x] = y[:1]  +'D' + y[5:]
            acid_list.append(y[:1]  +'D' + y[5:])
        elif y[2:5] == 'GLU':
            rename_dict[x] = y[:1] +'E' + y[5:]
            acid_list.append(y[:1]  +'E' + y[5:])
        elif y[2:5] == 'ARG':
            rename_dict[x] = y[:1] + 'R' + y[5:]
            base_list.append(y[:1]  +'R' + y[5:])
        elif y[2:5] == 'HIS':
            rename_dict[x] = y[:1] +'H' + y[5:]
            base_list.append(y[:1] +'H' + y[5:])
        elif y[2:5] == 'LYS':
            rename_dict[x] = y[:1]+ 'K' + y[5:]
            base_list.append(y[:1]+ 'K' + y[5:])
        elif y[2:5] == 'TYR':
            rename_dict[x] = y[:1] + 'Y' + y[5:]
            polar_rest_list.append(y[:1]+ 'Y' + y[5:])
        elif y[2:5] == 'CYS':
            rename_dict[x] = y[:1] + 'C' + y[5:]
            polar_rest_list.append(y[:1]+ 'C' + y[5:])
        elif y[2:5] == 'MQ8':
            rename_dict[x] ='MQ' + y[5:]
            ub_q_list.append('MQ' + y[5:])
        else:
            non_residue_list.append(y)
        
          
    col_order_list = acid_list + polar_rest_list+ base_list + ub_q_list+ non_residue_list
    file_input = file_input.rename(rename_dict, axis =1)
    file_input = file_input[col_order_list]
    return file_input

def dropCorrCriterion(df, cutoff = 0.0):
    data_frame= wc.WeightedCorr(df=df, wcol='Count')(method='pearson')
    for i in data_frame.columns:
        if list(abs(data_frame[i]) >= cutoff).count(True) == 1:
            data_frame.drop(i, inplace = True)
            data_frame.drop(i, axis =1, inplace = True)
    return data_frame

def corr_heat_map(df_corr):
    plt.figure(figsize=(25, 8))
    cmap = ListedColormap(["darkred", "red", "orange", "lightgray","skyblue", "blue", 'darkblue'])
    bounds= [-1.0, -0.5, -0.3, -0.1, 0.1, 0.3, 0.5, 1.0]
    norm =BoundaryNorm(bounds, cmap.N)
    heatmap = sns.heatmap(df_corr, annot=True, linecolor= "gray", cmap=cmap, norm = norm, square = True, fmt=".2f",linewidths=.01,annot_kws={"fontsize":12})
    plt.ylabel(None)
    plt.xlabel(None)
    plt.yticks(fontsize = 15, rotation = 0)
    plt.xticks(fontsize = 15, rotation = 90)
    cbar = heatmap.collections[0].colorbar
    cbar.ax.tick_params(labelsize=20)
    plt.savefig("corr.pdf", dpi = 300, bbox_inches = 'tight')
    logging.info("Correlation heat map is saved as: corr.pdf")

    return 

# This to rename based on cluster
def renameResidueCluster(df):
    rename_dict_cls = {}
    cls1 = ["HISE0038", "TYRE0087", "ASPE0139", "ARGE0217", "ARGQ0216", "GLUQ0223", "GLUQ0225"]
    cls2 = ["GLUE0050", "GLUE0051", "GLUE0216", "HISP0060", "ARGQ0154", "ASPQ0220", "GLUQ0227", "ARGQ0301", "TYRQ0302", "ASPQ0303", "ARGQ0307"]
    con1_2 = ["GLUE0388", "ASPE0392", "ARGQ0299"]
    cls3 = ["ARGE0042", "HISE0058", "TYRE0061", "GLUP0045", "ASPP0049", "GLUP0053", "ASPQ0071", "LYSQ0146", "TYRQ0147", "GLUQ0235", "TYRQ0236", "LYSQ0240"]
    cls4 = ["ASPP0072", "GLUP0074", "TYRP0109", "GLUP0110", "LYSP0113", "ARGP0117", "TYRQ0124", "GLUQ0130", "TYRQ0134", "TYRQ0162", "GLUQ0163", "TYRQ0206", "GLUQ0213", "HISQ0251", "GLUR0158"]
    con1_4 = ["GLUQ0035", "TYRQ0232", "HISQ0233", "GLUQ0248", "TYRQ0249", "ARGQ0294"]
    cls5 = ["TYRR0043", "TYRR0059", "GLUS0032"]
    cls6 = ["GLUP0006", "TYRP0007", "TYRQ0005", "ASPQ0008", "TYRQ0010", "LYSQ0016", "ASPQ0119", "ASPQ0184", "LYSQ0190", "GLUR0005", "ASPR0047", "ARGR0049", "ASPR0118", "LYSR0120", "TYRS0003", "ARGS0047"]
    rename_dict_cls = {}
    for i in df.columns:
        if i[:-1] in cls1:
            rename_dict_cls[i] = i + "cl1"
        elif i[:-1] in cls2:
            rename_dict_cls[i] = i + "cl2"
        elif i[:-1] in con1_2:
            rename_dict_cls[i] = i + "co1_2"
        elif i[:-1] in cls3:
            rename_dict_cls[i] = i + "cl3"
        elif i[:-1] in cls4:
            rename_dict_cls[i] = i + "cl4"
        elif i[:-1] in con1_4:
            rename_dict_cls[i] = i + "co1_4"
        elif i[:-1] in cls5:
            rename_dict_cls[i] = i + "cl5"

        elif i[:-1] in cls6:
            rename_dict_cls[i] = i + "cl6"
        else:
            rename_dict_cls[i] = i 
    df = df.rename(rename_dict_cls, axis =1)
    return df

def sum_cls_crg(df, cluster_id_list):
    for cluster_id in cluster_id_list:
        col_cl = df.columns[df.columns.str.contains(f"_{cluster_id}")]
        df[f"Crg_{cluster_id}"] = df.loc[:, col_cl].sum(axis = 1)
    return df


if __name__ == "__main__":
    
    # if you have ligands or other residues. Make sure to add here.

    your_interested_res = ["ASP", "GLU", "ARG", "HIS", "LYS", "CYS", "TYR", "NTR", "CTR"] 

   # open the microstate and head3 files
   #STRATAGY: A FOLER GIVEN BY dir0 THIS UNIQULY NAMED FOLDER
   #CURRENTLY THE .ms AND head3.lst FILES ALSO HAVE CONTAIN THIER FOLDER UNIQUE NAME
   #dir0 = '/Users/marilyngunner/Dropbox/Umesh/antiporter-mrg/antiporter/mrg-work/'
   #name = 'dry_apo_final'    #CHANGE THIS FOR UNIQUE INPUT
   #dir = dir0+name
   #mc = msa.MSout(dir+'/'+name+'.ms')
   #head3_conf = msa.read_conformers(dir+'/'+name+'_head3.lst')

#   mc = msa.MSout("ms_out/pH7eH0ms.txt")
    mc = msa.MSout("ms_out.ms")
    head3_conf = msa.read_conformers("head3.lst")

    ms_orig_lst = [[ms.E, ms.count, ms.state] for  ms in list((mc.microstates.values()))]
    ms_orig_lst = sorted(ms_orig_lst, key = lambda x:x[0])
    

    logging.info(f"Total ms steps (total count): {mc.N_ms}")
    logging.info(f"Total unique conformers ms: {mc.N_uniq}")
    histogram_ms_energy(ms_orig_lst)
    ms_free_residues = free_residue()
    ms_fixed_residues = fixed_residue()
    background_charge = ms_fixed_residues[0]
    fixed_residues_crg = ms_fixed_residues[1]
    id_vs_charge =  id_crg_relation()
    crg_orig_lst = convert_ms_crg(ms_orig_lst, id_vs_charge)

    charge_ms_file = findUniqueCrgmsCountOrder(crg_orig_lst)
    all_crg_count_res = ConcaCrgMsPandas(charge_ms_file[0], charge_ms_file[1], charge_ms_file[2], 
                                     ms_free_residues, background_charge,your_interested_res)

    plots_unique_crg_histogram(charge_ms_file, background_charge)

    combine_free_fixed = combine_free_fixed_residue(fixed_residues_crg, all_crg_count_res)
    df_join_free_back_rena = renameResidueCluster(combine_free_fixed)
    cluster_id_list = ["cl1", "cl2", "co1_2","cl3", "cl4", "co1_4", "cl5", "cl6"]
    df_all_crg_cls = sum_cls_crg(df_join_free_back_rena, cluster_id_list = cluster_id_list)
    

    all_crg_count_res.to_csv('all_crg_count_res.csv', header = True)

    all_crg_count_std = correlation_data_parsing(all_crg_count_res)
    df = renameOrderResidues(all_crg_count_std)
#CHANGE THE CORRELATION LEVEL HERE
    df_correlation= dropCorrCriterion(df, cutoff = 0.2) # define the cutoff here
    corr_heat_map(df_correlation)







