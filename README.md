# Interpangenome *E.coli* analysis
Scripts used for analysis in the inter-pangenome study.
## pangenome_splitter
A script to split a pangenome's representative cluster sequences into separate core and accessory cluster sequence lists based on a specified core threshold.
### Required input files
This script requires the `gene_presence_absence_roary.csv` and `pan_genome_reference.fa` files that are output from 
[Panaroo](https://gtonkinhill.github.io/panaroo/#/).
(Also works with the `gene_presence_absence.csv` and `pan_genome_reference.fa` output files from [Roary](https://sanger-pathogens.github.io/Roary/)).
### Usage
`pangenome_splitter.py -m gene_presence_absence_roary.csv -f pan_genome_reference.fa -o output -core 0.95`
## process_partial_hits
A script to process an ABRicate summary output file to account for partial gene hits.
### Required input files
This script requires an [ABRicate](https://github.com/tseemann/abricate) summary file in .tsv format.
### Usage
`process_partial_hits.py -i summary_file.tsv -o output`
## function_distribution_clustermap
A script to visualise the percentage of genes per COG functional category between multiple lineages and output .csv containing the breakdown of percentages. 
### Required input files
This script requires the path to a single directory containing annotation files output from [EggNOG-mapper](https://github.com/eggnogdb/eggnog-mapper) in .tsv format.
### Usage
`function_distribution_clustermap.py -p path/to/your/directory -o output`
