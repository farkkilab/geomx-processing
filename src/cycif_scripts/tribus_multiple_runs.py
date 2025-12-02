import os
import re
from pathlib import Path
import pandas as pd
import numpy as np


import tribus
from visualization import heatmap_for_median_expression, marker_expression, umap_vis, z_score, cell_type_distribution


# set up paths
pdrive_mount_dir = "/run/user/1001/gvfs/smb-share:server=group3.ad.helsinki.fi,share=h30492/"

segm_cells_dir = os.path.join(pdrive_mount_dir, "farkkilab2/7_TLS/Data/Exp_4/csv")
segm_cells_qced_dir = os.path.join(pdrive_mount_dir, "farkkilab2/9_EyeMT/Data_analysis/cycif_processing/phenotyping/tribus_multiple_runs_test/input_segm_cells_qced")
pt_names = ["7-S026_iOme1_01253_4A", "7-S026_iOme2_01253_4A", "7-S046_iOme_01253_4A"]

logic_table_dir = os.path.join(pdrive_mount_dir, "farkkilab2/9_EyeMT/Data_analysis/cycif_processing/phenotyping/tribus_multiple_runs_test/tribus_logic_tables/global_immune")
logic_table_list = [os.path.join(dp, f) for dp, dn, fn in os.walk(os.path.expanduser(logic_table_dir)) for f in fn]

output_dir = os.path.join(pdrive_mount_dir, "farkkilab2/9_EyeMT/Data_analysis/cycif_processing/phenotyping/tribus_multiple_runs_test/tribus_output_global_immune_norm_filt_markers2")
#segm_cells_pt_list = [os.path.join(segm_cells_dir, pt + '.csv') for pt in pt_names]

# set up depth param
depth = 2

markers_list = ['PanCK', 'Vimentin', 'Iba1', 'CD11c', 'CD4', 'CD8a', 'NKG2a'] # matching logic tables
#markers_list_colnames = ['PanCK_2', 'Vimentin_1', 'Iba1_1', 'CD11c_1', 'CD4_1', 'CD8a_2', 'NKG2A_2'] # raw signal
# normalised + filtered signal
markers_list_colnames = ['PanCK_normalized', 'Vimentin_normalized', 'Iba1_filt_norm', 'CD11c_normalized', 
                         'CD4_filt_norm', 'CD8a_filt_norm', 'NKGA_filt_norm']

################################

for pt_name in pt_names:
    for logic_table_path in logic_table_list:
        # raw signal csv
        # segm_cells_pt_path = os.path.join(segm_cells_dir, pt_name + '.csv')
        # filtered from norm+filt h5ad object
        segm_cells_pt_path = os.path.join(segm_cells_qced_dir, pt_name + '_segm_qced2.csv')
        output_name = pt_name + Path(logic_table_path).stem
        
        #print(segm_cells_pt_path)
        #print(logic_table_path)
        print(output_name)
        
        #load cell csv
        sample_data = pd.read_csv(segm_cells_pt_path)
        sample_data_selected = sample_data.loc[:, 
        ['CellID', 'Y_centroid', 'X_centroid', 'Area', 'Eccentricity'] + markers_list_colnames]
        # change colnames to match logic
        sample_data_selected.columns = ['CellID', 'Y_centroid', 'X_centroid', 'Area', 'Eccentricity'] + markers_list
        
        print("cell csv loaded")
        #load logic xls
        #logic1_name = re.findall(re.compile(r'global_[a-z]*'), Path(logic_table_path).stem)[0]
        #logic2_name = re.findall(re.compile(r'immune_.*'), Path(logic_table_path).stem)[0]
        
        #logic1 = pd.read_excel(logic_table_path, sheet_name= logic1_name)
        #logic2 = pd.read_excel(logic_table_path, sheet_name= logic2_name)
        df = pd.ExcelFile(logic_table_path)
        logic = {"Global": pd.read_excel(df, df.sheet_names[0], index_col=0), "Immune": pd.read_excel(df, df.sheet_names[1], index_col=0)}

        print("logic table loaded")

        # set random seed to ensure reproducing
        labels, scores = tribus.run_tribus(np.arcsinh(sample_data_selected/5.0), logic, depth=depth, normalization=z_score, 
                                    tuning=0, sigma=1, learning_rate=1, 
                                    clustering_threshold=100, undefined_threshold=0.0005, other_threshold=0.4, random_state=42)
        
        result_data = sample_data_selected.join(labels)
        labels_new = labels.join(result_data["CellID"])
        labels_new.to_csv(os.path.join(output_dir, output_name + '_tribus_annotation.csv'))
        result_data.to_csv(os.path.join(output_dir, output_name + '_raw_tribus_annotated.csv'))
        
        print('output saved')
        print('%%%%')
        
