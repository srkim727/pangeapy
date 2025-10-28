# 🌏 pangeapy

`pangeapy` is a Python package designed for **automated cell type annotation** and **metadata prediction** using the **PANGEA reference atlas**.
It enables researchers to easily label single-cell transcriptomic data and predict higher-level phenotypic contexts such as organ or disease state, with minimal preprocessing.

## Explore PANGEA through a web interface

You can explore the interactive documentation at **[https://srkim727.github.io/](https://srkim727.github.io/)**, where you can:

1. **Annotate cells** — Automatically assign fine-grained cell type labels using the PANGEA reference model.  
2. **Predict metadata** — Infer higher-level attributes such as organ, tissue, or multicellular phenotype.  
3. **Identify missing cell types** — Detect potential novel or unrepresented cell populations through uncertainty estimation.

## Installation

### 1) Create a Conda environment

```bash
conda create -n pangea python=3.10
conda activate pangea
```

### 2) Install via GitHub

```bash
pip install --user git+https://github.com/srkim727/pangeapy.git
```

## Usage

Below are example workflows for major modules in `pangeapy`.

### 1) Cell Annotation

The **CellAnnotator** module assigns cell-type labels to every cell in your single-cell dataset (`AnnData` format).

#### Example

```python
from pangeapy import CellAnnotator

cell_anno = CellAnnotator()
pred = cell_anno.annotate(adata, sample_key='Sample')
```

#### Input requirements

* Input should be an **AnnData** object (`adata`)
* Counts must be **normalized to 1e4 and log1p-transformed**

#### Output

The function returns a DataFrame where each cell is annotated across multiple hierarchical levels:

|index|Level1&#124;predicted_labels   |Level1&#124;over_clustering|Level1&#124;majority_voting|Level1&#124;conf_score|Level1&#124;cert_score|Level2&#124;predicted_labels|Level2&#124;over_clustering|Level2&#124;majority_voting|Level2&#124;conf_score|Level2&#124;cert_score|PG_annotations|PG_combined_score|Sample  |
|--------------------------|----------------------|----------------------|-----------------|-----------------|-----------------------|----------------------|----------------------|-----------------|-----------------|--------------|-----------------|--------|------------------|
|AAACCTGAGAAAGTGG-1-gPlexA1|T&NK                  |4                     |T&NK             |0.999999         |0.275238               |NK_CD16               |36                    |NK_CD16          |0.999932         |0.216732      |T&NK&#124;NK_CD16     |0.999966|S00109-Ja001E-PBCa|
|AAAGATGGTTCCACTC-1-gPlexA1|T&NK                  |29                    |T&NK             |0.999993         |0.503419               |NK_CD16               |35                    |NK_CD16          |0.992871         |0.344663      |T&NK&#124;NK_CD16     |0.996425|S00109-Ja001E-PBCa|

> Each row corresponds to one cell and includes confidence and certainty scores for hierarchical annotations.

### 2) Metadata Prediction

The **MetaAnnotator** module predicts higher-level sample attributes such as **organ identity** and **multicellular blood phenotypes**.

#### Example

```python
from pangeapy import MetaAnnotator

meta_anno = MetaAnnotator()
meta = meta_anno.annotate(pred, sample_key='Sample')
merged = meta.integrate()
```

#### Output

| index              | Organ_pred | Organ_prob | Pheno_pred | Pheno_prob |
| ------------------ | ---------- | ---------- | ---------- | ---------- |
| S00052-Ja005E-PBCa | Blood      | 0.949866   | BP2        | 0.535271   |
| H00067-Ha001E-PBGa | Blood      | 0.939722   | BP4        | 0.997693   |
| S00028-Ja001E-PBCa | Blood      | 0.996022   | BP7        | 0.751714   |

> The results can be merged with cell annotation outputs for downstream analysis or visualization.

### 3) Identifying Missing Cell Types

To detect potentially **unrepresented or missing cell types**,
enable uncertainty estimation in the annotation step
by setting `compute_uncertainty=True`:

```python
cell_anno = CellAnnotator(compute_uncertainty=True)
pred = cell_anno.annotate(adata)
```

## Tutorials

| Task                           | Notebook                                                                                                                                 |
| ------------------------------ | ---------------------------------------------------------------------------------------------------------------------------------------- |
| 🔹 Cell & Metadata Annotation  | [vignette_annotation.ipynb](https://github.com/srkim727/pangeapy/blob/main/docs/vignette_annotation.ipynb)                               |
| 🔹 Missing Cell Type Detection | [vignette_identifying_missing_cells.ipynb](https://github.com/srkim727/pangeapy/blob/main/docs/vignette_identifying_missing_cells.ipynb) |

## Citation

Kim, *unpublished.*

