# EnrichmentatorList
Enrichmentator version, for gene list.

This series of scripts perform enrichment analysis on a directory with several files containing lists of genes. 
The enrichment is performed on database KEGG and GeneOntology (biological process, molecular function, and cellular component).
The output is a heatmap showing the processes of each database as rows and columns as each gene list (labeled by its filename).

For use run the shell script Enrichment.sh, and as arguments first the directory containing the lists of genes (genelist.txt) and second p-value for the analysis.

There must be at least two gene lists and such files must have the extension ".txt".
If there is no process that enrisquesca any comundiad, there will be an error.
