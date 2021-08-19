#! /usr/bin/Rscript
install.packages("pbmcapply")
library(gtools)
library(parallel)
library(pbmcapply)
library(tidyverse)
library(Seurat)
library(data.table)

message("starting commands with ", detectCores(), " cores available")

metafile = Sys.getenv("METAFILE")
srobjfile = Sys.getenv("SROBJFILE")
column = Sys.getenv("COLUMN")
outfile = Sys.getenv("OUTRDSFILE")

meta = read.table(metafile, header=1) %>%
    column_to_rownames("CellID")
srobj = readRDS(srobjfile)
srobj@meta.data = meta[colnames(srobj), ]

message("read in data")

pairs <- t(combn(levels(srobj@meta.data[[column]]), 2)) %>%
    data.frame(X3 = paste0(.[,1], "__VS__", .[,2])) %>%
    setNames(c('ident1', 'ident2', 'comparison_name'))
rownames(pairs) <- pairs$comparison_name
pairs.t <- as.data.frame(t(pairs))

message("created Pairs list")

Idents(srobj) = srobj@meta.data[[column]]
pw.de.lst.rna = pbmclapply(pairs.t,
                       function(r) {
                           print(r)
                           rownames_to_column(FindMarkers(srobj, assay="RNA",
                                ident.1 = r['ident1'],
                                ident.2 = r['ident2'],
                                test.use = "wilcox", min.pct=0.1,
                                max.cells.per.ident = 500), "gene")
                       }, mc.cores=detectCores(), ignore.interactive = T)

message("finished differential expression tests!")

pw.de.df.rna = data.table::rbindlist(pw.de.lst.rna, idcol="comparison_name")

message("finished converting to table")

saveRDS(pw.de.df.rna, outfile)

message("saved file")


