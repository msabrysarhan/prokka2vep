# prokka2vep

## Introduction

prokka (https://github.com/tseemann/prokka)
VEP (https://www.ensembl.org/info/docs/tools/vep/index.html)
This script converts prokka GFF files to VEP-friendly GFF format. 

## Installation

### Dependencies

#### Mandatory dependencies:
1. Python >= 3.7
2. Pandas >= 2.0
3. python-csv == 0.0.13


## Usage

```bash 
python3 prokka2vep.py --gff SGB4837.gff --out vep_gff.gff
```

```
Reading GFF file as pandas dataframe ...

Creating transcript records ...

Merging GFF dataframes ... 

Reordering GFF rows ... 

Processing the non-coding RNA ... 

Writing GFF file to vep_gff.gff 

GFF conversion is done. Thanks for using prokka2vep!
```
