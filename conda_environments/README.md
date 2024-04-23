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