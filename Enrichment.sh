#!/bin/bash

# params: DirList pvalue (0.001)

mkdir RESULTADOS
d=$1

for j in $(ls $d/*.txt)
  do
    Rscript lib/Enrichmentator.R $j $2
  done

  python lib/EnrichmentMatrix.py $d
  
  Rscript lib/heatmap4.r $d'_GS_GO_BP.csv'
  Rscript lib/heatmap4.r $d'_GS_GO_CC.csv'
  Rscript lib/heatmap4.r $d'_GS_GO_MF.csv'
  Rscript lib/heatmap4.r $d'_KEGG.csv'
   
  mkdir $d'EnrichmentMatrices'
  mv *.csv $d'EnrichmentMatrices'
  mkdir RESULTADOS/$d
  mv *.pdf RESULTADOS/$d
