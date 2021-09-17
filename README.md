*This README document describes the script get_OG_fastas.py. A generic workflow description is provided in workflow_description.docx*

# INPUT
- input_table

CSV table with header and at least the columns "gene" and "locus_tag".\
Example:
```
gene,locus_tag
gene1,locus_tag1
gene2,locus_tag2
gene2,locus_tag3
gene3,locus_tag4
```

- gbk_dir

Folder with assemblies in genbank format (.gbk or .gbff).\
Example folder structure:
```
gbk_dir
   |-assembly1.gbk
   |-assembly2.gbk
   |-assembly3.gbff
```

- og_matrix

Tab-delimited matrix with OGs in rows, and assemblies in columns. Each cell lists which locus_tag from which assemblies
belong to which OG. If multiple locus tags in an assembly belong to the same OG, they must be whitespace-separated in the cell.\
Example:
```
    assembly1   assembly2
OG1 locus11     locus12 locus22
OG2             locus21
OG3 locus31
```
# Installation:

1) Clone the git repo (or simply download / unpackage the repo archive) and cd to to project folder
2) Recreate the environment (with conda or mamba), then activate the environment:

- mamba env create -f envs/biopandas.yaml
- conda activate biopandas

OR

- conda env create -f envs/biopandas.yaml
- conda activate biopandas

# Test

Requires the input files specified in Input (input table csv, folder of gbk's, OG matrix); file can be place in a 'data' folder such that a test command is:

./scripts/get_OG_fastas.py -i data/input_table.csv  -g data/gbks -m data/matrix_annotated.txt -o results/test01

# Notes:

-The output files are named AFTER the gene names. This means that if multiple gene names in the input table correspond to the same OG, these output files will be identical, but have different file names.
-If there are duplicate genenames, they are renamed with the suffix _DUP_1, _DUP_2, etc.
