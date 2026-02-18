import os
import re
from pathlib import Path
import pandas as pd
import numpy as np

import tribus
from visualization import heatmap_for_median_expression, marker_expression, umap_vis, z_score, cell_type_distribution

# INPUT:
# 1. segm_cells_qced files containing qced segmented cells along with filtered+normalised signal for (at least) channels specified in marker_list_colnames
# should contain at lest the following columns: CellID,sample_id,Sample,X_centroid,Y_centroid, marker_list_colnames
# 2. tribus logic table with depth specified in depth parameter and specific sheets in order specified in logic object (line 70)

# set up paths
pdrive_mount_dir = "F:"

#"/run/user/1001/gvfs/smb-share:server=group3.ad.helsinki.fi,share=h30492/"
#pdrive_mount_dir = "/home/ad/P-drive/h30492/"

segm_cells_qced_dir = os.path.join(pdrive_mount_dir, "farkkilab2//9_EyeMT//Data_analysis//cycif_processing//phenotyping//tribus_multiple_runs_test//input_segm_cells_qced_b3_missing")
segm_cells_qced_list = [os.path.join(dp, f) for dp, dn, fn in os.walk(os.path.expanduser(segm_cells_qced_dir)) for f in fn]


logic_table_path = os.path.join(pdrive_mount_dir, "farkkilab2//9_EyeMT//Data_analysis//cycif_processing//phenotyping//tribus_multiple_runs_test//tribus_logic_tables//global_stroma_immune_single_best//eyemt_tribus_global_neg_stroma_neg_tumor_zero_immune_zero_nkdrop.xlsx")

output_dir = os.path.join(pdrive_mount_dir, "farkkilab2//9_EyeMT//Data_analysis//cycif_processing//phenotyping//tribus_multiple_runs_test//tribus_output_single_best_b3")

logic_table_path = "F://farkkilab2//9_EyeMT//Data_analysis//cycif_processing//phenotyping//tribus_multiple_runs_test//tribus_logic_tables//global_stroma_immune_single_best//eyemt_tribus_global_neg_stroma_neg_tumor_zero_immune_zero_nkdrop.xlsx"
segm_cells_qced_list "F://farkkilab2//9_EyeMT//Data_analysis//cycif_processing//phenotyping//tribus_multiple_runs_test//input_segm_cells_qced_b3_missing"
output_dir = "F://farkkilab2//9_EyeMT//Data_analysis//cycif_processing//phenotyping//tribus_multiple_runs_test//tribus_output_single_best_b3"

print(segm_cells_qced_list)

#####################################################
# set up depth param
depth = 3
# colnames matching logic tables
markers_list = ['PanCK', 'Vimentin', 'Iba1', 'CD11c', 'CD4', 'CD8a'] 

# colnames matching input signal intensity tables in input segmented cells qced csvs - ENSURE ORDER with markers_list
# normalised + filtered signal
markers_list_colnames = ['PanCK_normalized', 'Vimentin_normalized', 'Iba1_normalized', 'CD11c_filt_norm', 
                         'CD4_filt_norm', 'CD8a_filt_norm']

# suffix of segmented cells qced csvs - will be removed to get sample name
segm_cells_qced_suffix = '_segm_qced_selected_markers'

################################

if not os.path.exists(output_dir):
    os.makedirs(output_dir)

################################

# iterate through patients qc-ed segmented cells and run tribus
for segm_cells_pt_path in segm_cells_qced_list:

    pt_name = Path(segm_cells_pt_path).stem
    pt_name = pt_name.replace(segm_cells_qced_suffix, "")

    print(pt_name)
   
    # load cell csv
    sample_data = pd.read_csv(segm_cells_pt_path)
    sample_data_selected = sample_data.loc[:, 
    ['CellID', 'Y_centroid', 'X_centroid', 'Area', 'Eccentricity'] + markers_list_colnames]
        
    # change colnames to match logic
    sample_data_selected.columns = ['CellID', 'Y_centroid', 'X_centroid', 'Area', 'Eccentricity'] + markers_list

    print("qced cells loaded")

    # load logic xls
    df = pd.ExcelFile(logic_table_path)

    logic = {"Global": pd.read_excel(df, df.sheet_names[0], index_col=0),
             "Stroma": pd.read_excel(df, df.sheet_names[1], index_col=0),
             "Tumor": pd.read_excel(df, df.sheet_names[2], index_col=0),
             "Immune": pd.read_excel(df, df.sheet_names[3], index_col=0),
             "Stroma_Immune": pd.read_excel(df, df.sheet_names[4], index_col=0),
             "Tumor_Immune": pd.read_excel(df, df.sheet_names[5], index_col=0)}

    print("logic table loaded")

    # set random seed to ensure reproducing
    labels, scores = tribus.run_tribus(np.arcsinh(sample_data_selected/5.0), logic, depth=depth, normalization=z_score, 
                                       tuning=0, sigma=1, learning_rate=1, 
                                       clustering_threshold=100, undefined_threshold=0.0005, other_threshold=0.4, random_state=42)
        
    result_data = sample_data_selected.join(labels)
    labels_new = labels.join(result_data["CellID"])
    labels_new.to_csv(os.path.join(output_dir, pt_name + '_tribus_annotation.csv'))
    result_data.to_csv(os.path.join(output_dir, pt_name + '_raw_tribus_annotated.csv'))
    scores.to_csv(os.path.join(output_dir, pt_name + '_tribus_scores.csv'))
        
    print('output saved')
    print('%%%%')
