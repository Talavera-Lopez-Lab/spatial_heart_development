# Environments Installation
## YML File
Export Environment
```bash
mamba env export --no-builds > path/to/environment.yml
```
Create environment
```bash
mamba env create -f path/to/environment.yml
```
## Manual
### cellxgene-env
```bash
mamba create -n cellxgene-env python=3.10
mamba activate cellxgene-env
pip install -U cellxgene-census
```
Notebooks:
- 0-240423-Downloading_Data/0-240423-Cellxgene_api_querying.ipynb
- 0-240423-Downloading_Data/1-240423-concat_anndata.ipynb