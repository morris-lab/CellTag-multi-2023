import scanpy as sc
import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.gridspec import GridSpec

import plotly.express as px

sc.settings.set_verbosity=3
sc.set_figure_params(figsize=(3,3))

import sf_utils

from spectra import spectra as spc
import spectra
import json


adata = sc.read_h5ad("../06-integrate_kd_oe/final_objs/cc_regressed_intobj_cr_aggr.h5ad")
adata.X = adata.raw.X

with open("./spectra_dicts/spectra_rd2_markers.json") as f:
    annots = json.load(f)
    
adata.obs.seurat_clusters_old = adata.obs.seurat_clusters_old.astype('str')
model = spc.est_spectra(adata = adata,
                        gene_set_dictionary = annots,
                        cell_type_key = "seurat_clusters_old",
                        use_highly_variable = True, lam = 0.01)

model.save("./spectra_models/spectra_intobj_rd2_model")
adata.write_h5ad("./sc_objects/cc_regressed_intobj_cr_aggr_spectra_rd2_cellchat.h5ad")

##sbatch -c 24 --mem-per-cpu=3G --mail-type=ALL --mail-user=$USER@wustl.edu -J spectra_cellchat_rd2 --output logs/spectra-%J.out --wrap=". /opt/apps/labs/smlab/software/py_env/sc_env/bin/activate; python spectra_rd2_cellchat.py"