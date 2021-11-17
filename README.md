# TITLE
All scripts used for analysis in the inter-pangenome study.
## pangenome_splitter
A script to split a pangenome's representative cluster sequences into separate core and accessory cluster sequence lists based on a given core threshold.
### Required input files
This script requires the `gene_presence_absence_roary.csv` and `pan_genome_reference.fa` files that are output from 
[Panaroo](https://gtonkinhill.github.io/panaroo/#/).
(Also works with the `gene_presence_absence.csv` and `pan_genome_reference.fa` output files from [Roary](https://sanger-pathogens.github.io/Roary/).)
### Usage
`pangenome_splitter.py -m gene_presence_absence_roary.csv -f pan_genome_reference.fa -o output -core 0.95`
