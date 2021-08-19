# PREDICT (2021) Code Repo

This repository contains code and analysis to replicate the findings of [Zheng et al]()

## Contents

Folders and files of the repo are numbered roughly in order of dependency, such that a higher number requires the results of lower numbered analyses. Numbering restarts within each folder so files should be run from lowest number to highest. Non-numbered notebooks can be run in any order after the numbered files. Other non-numbered files are outputs or direct dependencies of numbered batch scripts.

This repo shows how we:

* Run interactive tiered clustering with ARBOL (see external tool for updated instructions for latest version)
* Make standardized names for all end clusters
* Generate a binary tree representing the transcriptomic relations between clusters
* Map between FGID and CD clusters using RandomForest model predictions
* Correlate clinical metadata to transcriptomic information
* Use GSEA to isolate pathways connected clinical disease severity
* Analyze transcriptomic landscape of Myeloid and T cells for cell state differences connected to clinical disease severity.
* Create marker gene lists for difference levels of the cell-type hierarchy
* Compare to traditional clustering methods

## Links to external tutorials and interactive tools

ARBOL tutorial: https://shaleklab.github.io/SCTieredClustering/

Annotation Rank Score shiny app: https://kylekimler.shinyapps.io/shinyrank/

GSEA analysis: https://jo-m-lab.github.io/3p-PREDICT-Paper/4_GSEA/PREDICT_GSEA_final.html
