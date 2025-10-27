# pangeapy

## Installation

### create conda environment
conda create -n pangea python=3.10 

conda activate pangea   

### install python package through pip/github
pip install --user git+https:<area>//github.com/srkim727/pangeapy.git   

## Usage
you can find an example of usage from the following link
cell and meta annotation: https://github.com/srkim727/pangeapy/blob/main/docs/vignette_annotation.ipynb   
calling missing cell types: https://github.com/srkim727/pangeapy/blob/main/docs/vignette_identifying_missing_cells.ipynb

1) cell annotation   
   
  loading cell annotation modules   
  ```cell_anno = CellAnnotator()```   
   
cell annotation modules will label every cell in the scRNA-seq data in AnnData format   
(counts should be 10^4 normalized and log1p-transformed)   
```pred = cell_anno.annotate(adata, sample_key='Sample')```      

generated is a dataframe of cell indexes with hierarchical cell annotations and prediction probabilities

|index|Level1&#124;predicted_labels   |Level1&#124;over_clustering|Level1&#124;majority_voting|Level1&#124;conf_score|Level1&#124;cert_score|Level2&#124;predicted_labels|Level2&#124;over_clustering|Level2&#124;majority_voting|Level2&#124;conf_score|Level2&#124;cert_score|PG_annotations|PG_combined_score|Sample  |
|--------------------------|----------------------|----------------------|-----------------|-----------------|-----------------------|----------------------|----------------------|-----------------|-----------------|--------------|-----------------|--------|------------------|
|AAACCTGAGAAAGTGG-1-gPlexA1|T&NK                  |4                     |T&NK             |0.999999         |0.275238               |NK_CD16               |36                    |NK_CD16          |0.999932         |0.216732      |T&NK&#124;NK_CD16     |0.999966|S00109-Ja001E-PBCa|
|AAAGATGGTTCCACTC-1-gPlexA1|T&NK                  |29                    |T&NK             |0.999993         |0.503419               |NK_CD16               |35                    |NK_CD16          |0.992871         |0.344663      |T&NK&#124;NK_CD16     |0.996425|S00109-Ja001E-PBCa|


2) metadata prediction

loading metadata annotation module   
```meta_anno = MetaAnnotator()```

metadata annotation modules predict each sample's organ identity and disease-associated phenotypes   
```meta = meta_anno.annotate(pred, sample_key='Sample')```   
```merged = meta.integrate()```

this process will generate a dataframe of samples with predicted organ identities and phenotypes   

|index                    |Organ_pred   |Organ_prob|Pheno_pred|Pheno_prob|
|--------------------------|-------------|----------|----------|----------|
|S00052-Ja005E-PBCa        |Blood        |0.949866  |BP2       |0.535271  |
|H00067-Ha001E-PBGa        |Blood        |0.939722  |BP4       |0.997693  |
|S00028-Ja001E-PBCa        |Blood        |0.996022  |BP7       |0.751714  |
|N00032-Ja001E-PBGa        |Blood        |0.845979  |BP7       |0.778266  |
|H00054-Ha001E-PBGa        |Blood        |0.889262  |BP1       |0.852949  |


3) identification of missing cell types

use compute_uncertainty = True option in cell_anno model
please refer to https://github.com/srkim727/pangeapy/blob/main/docs/vignette_identifying_missing_cells.ipynb

## Citation
Kim, unpublished
