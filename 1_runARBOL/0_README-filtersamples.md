To filter each 10X sample file individually we used batch .wdl pipline, The pipeline
can be found https://dockstore.org/workflows/github.com/BenjaminDoran/Filter10XData/Filter10XData:0.0.2dev

And through https://firecloud.terra.bio

We used it because we had a lot of samples and we were testing out the newer clouds resources available at the time. Work flow The basic workflow distills to calling these functions to filter cells based on standard qc thresholds used in scanpy and seurat.

```python
import os
import yaml
from pprint import pprint

import scanpy.api as sc
import pandas as pd
import deepdish as dd

def filter_data_samples_00(raw_paths, out_dir, meta_file, user_thresholds={}, verbose=True):
    """Filters cells and genes from provided sample paths
    Args:
        raw_paths (list): list of filepaths to raw h5ad samples
        out_dir (str): filepath to folder to write filtered samples
        meta_file (str): filepath write h5 meta data file.
        user_thresholds (dict, optional): Defaults to {}. thresholds for each
            sample in format {<smp_id>: {<param>: value}}. See help in
            `utils.prc.filter_adata` for more details
        verbose (bool, optional): Defaults to True. prints status while running
    Returns:
        bool: True if all samples filtered successfully
    """
    ad = {}
    for f in sorted(raw_paths):
        smp = f.split("/")[-1][:-5]
        verbose and print(f"\rReading {smp}          ", end="")
        ad[smp] = sc.read(f)
        ad[smp].var_names_make_unique()
    verbose and print("\rAll samples read in!          ")

    ad_filt = {}
    meta_filt = {}
    for smp in ad:
        verbose and print(f"\rFiltering {smp}           ", end="")
        ad_filt[smp], meta_filt[smp] = filter_adata(
            ad[smp], user_thresholds.get(smp, {}))
        ad_filt[smp].write(f"{out_dir}/{smp}.h5ad")

    verbose and print("\rFiltered all samples!          ")
    os.makedirs(os.path.dirname(meta_file), exist_ok=True)
    dd.io.save(meta_file, meta_filt)
    return True

def filter_adata(adata_raw, user_thresholds={}):
    """Filters low quality cells and genes from raw AnnData object.
    Args:
        adata_raw (AnnData): raw cell X gene dataset with metadata
        user_thresholds (dict): thresholds to filter cells and genes.
            Defaults to `utils.prc.default_filter_thresholds`
                {
                    "max_count": 15000,
                    "max_genes": 3000,
                    "min_cells": 3,
                    "min_genes": 300
                }
            Supported keys: 'min_cells',
                            'max_cells',
                            'max_count',
                            'max_genes'
    Returns:
        Tuple:
            0 (AnnData): filtered AnnData object
            1 (dict): meta data associated with filtering
                    {
                        "n_cell_raw": # num cells raw
                        "n_gene_raw": # num genes raw
                        "n_cell_min": # num cells after min filter
                        "n_gene_min": # num genes after min filter
                        "n_cell_max": # num cells after max filter
                        "n_gene_max": # num genes after max filter
                        "pct_cells_min": # % cells filtered in min filter
                        "pct_cells_max": # % cells filtered in max filter
                        "min_barcodes": # cell barcodes kept in min filter
                        "max_barcodes": # cell barcodes kept in max filter
                    }
    """

    params = default_filter_thresholds.copy()
    params.update(user_thresholds)

    # min filter
    adata_min = adata_raw.copy()
    print(adata_min)
    adata_min = add_qc_stats(adata_min)
    sc.pp.filter_cells(adata_min, min_genes=params['min_genes'])
    sc.pp.filter_genes(adata_min, min_cells=params['min_cells'])

    # max filter
    adata_max = adata_min.copy()
    sc.pp.filter_cells(adata_max, max_counts=params['max_count'])
    sc.pp.filter_cells(adata_max, max_genes=params['max_genes'])


    meta = {
        "min_barcodes": adata_min.obs.index,
        "max_barcodes": adata_max.obs.index
    }
    return adata_max, meta
```


The filters for the samples used in this paper are as follow:

Exact thresholds used for each sample:
| Sample_ID | Diagnosis | min_cells | min_genes | max_genes | max_count | N_cells_post_filter |
| --- | --- | --- | --- |--- | --- | --- |
| p011_T0D_ILE_LPS_3p | CD | 3 | 200 | 3000 | 20000 | 11087 |
| p014_T0D_ILE_LPS_3p | CD | 3 | 200 | 3000 | 20000 | 7285 |
| p016_T0D_ILE_LPS_3p | CD | 3 | 200 | 3000 | 20000 | 7713 |
| p018_T0D_ILE_LPS_3p | CD | 3 | 300 | 4000 | 20000 | 15469 |
| p019_T0D_ILE_LPS_3p | CD | 3 | 200 | 5000 | 20000 | 11200 |
| p022_T0D_ILE_LPS_3p | CD | 3 | 200 | 3000 | 20000 | 7798 |
| p024_T0D_ILE_LPS_3p | CD | 3 | 200 | 3000 | 20000 | 9705 |
| p026_T0D_ILE_LPS_3p | CD | 3 | 300 | 4000 | 30000 | 10847 |
| p027_T0D_ILE_LPS_3p | CD | 3 | 200 | 4000 | 40000 | 8124 |
| p031_T0D_ILE_LPS_3p | CD | 3 | 200 | 3000 | 20000 | 12177 |
| p032_T0D_ILE_LPS_3p | CD | 3 | 200 | 3000 | 20000 | 9052 |
| p041_T0D_ILE_LPS_3p | CD | 3 | 200 | 4000 | 40000 | 7508 |
| p042_T0D_ILE_LPS_3p | CD | 3 | 200 | 4000 | 30000 | 12219 |
| p048_T0D_ILE_LPS_3p | CD | 3 | 200 | 4000 | 40000 | 9158 |
| p009_T0D_ILE_LPS_3p | FGID | 3 | 200 | 3000 | 20000 | 11426 |
| p017_T0D_ILE_LPS_3p | FGID | 3 | 200 | 3000 | 20000 | 11107 |
| p023_T0D_ILE_LPS_3p | FGID | 3 | 200 | 3000 | 20000 | 11262 |
| p028_T0D_ILE_LPS_3p | FGID | 3 | 200 | 4000 | 20000 | 10293 |
| p030_T0D_ILE_LPS_3p | FGID | 3 | 200 | 4000 | 30000 | 9590 |
| p033_T0D_ILE_LPS_3p | FGID | 3 | 200 | 3000 | 20000 | 9514 |
| p034_T0D_ILE_LPS_3p | FGID | 3 | 200 | 3000 | 20000 | 8142 |
| p035_T0D_ILE_LPS_3p | FGID | 3 | 200 | 3000 | 20000 | 5668 |
| p037_T0D_ILE_LPS_3p | FGID | 3 | 200 | 3000 | 20000 | 4210 |
| p043_T0D_ILE_LPS_3p | FGID | 3 | 200 | 5000 | 30000 | 5045 |
| p044_T0D_ILE_LPS_3p | FGID | 3 | 350 | 4000 | 30000 | 12727 |
| p049_T0D_ILE_LPS_3p | FGID | 3 | 200 | 3000 | 20000 | 8003 |
| p050_T0D_ILE_LPS_3p | FGID | 3 | 200 | 4000 | 30000 | 8582 |



merging was then performed with this code:

```R
library(Seurat)
filtered.FGpatient.h5ad.files = sampledata$h5adfiles[
    sampledata$patent %in% paste0("p0", c("09", 17, 23, 28, 30, 33, 34, 35, 37, 43, 44, 49, 50))
]
srobjs = lapply(filtered.FGpatient.h5ad.files, readH5AD)
srobj = merge(x = srobjs[[1]],
              y = unlist(srobjs[2:length(srobjs)]),
              add.cell.ids = smp_ids,
              project = ANLS_LBL)
srobj$orig.ident <- Idents(srobj)
saveRDS(srobj, Sys.getenv("FGSROBJfile"))
```

```R
library(Seurat)
filtered.CDpatient.h5ad.files = sampledata$h5adfiles[
    sampledata$patent %in% paste0("p0", c(11, 14, 16, 18,19, 22, 24, 26, 27, 31, 32, 41, 42, 48))
]
srobjs = lapply(filtered.CDpatient.h5ad.files, readH5AD)
srobj = merge(x = srobjs[[1]],
              y = unlist(srobjs[2:length(srobjs)]),
              add.cell.ids = smp_ids,
              project = ANLS_LBL)
srobj$orig.ident <- Idents(srobj)
saveRDS(srobj, Sys.getenv("CDSROBJfile"))
```