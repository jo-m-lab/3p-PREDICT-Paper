#! /usr/bin/Rscript
require(Seurat)
require(tidyverse)
# source("ITCcode_v0.1.R")
source("http://gist.githubusercontent.com/BenjaminDoran/c5d31abec78d4529575b3c5b676e2129/raw/9b91709582dfc543ea2b214483b3ca0d9690d872/GenerateTieredClusters.R")
srobj <- readRDS(Sys.getenv("STARTING_SROBJ"))

bcells1 <- subset(srobj, idents="0")
tiers <- GenTieredClusters(srobj,
                           saveSROBJdir = sprintf("%s/srobjs_clust0", Sys.getenv("OUT")),
                           figdir = sprintf("%s/figs_clust0", Sys.getenv("OUT")),
                           SaveEndNamesDir = sprintf("%s/endclusts_clust0", Sys.getenv("OUT")))

bcells2 <- subset(srobj, idents="1")
tiers <- GenTieredClusters(srobj,
                           saveSROBJdir = sprintf("%s/srobjs_clust1", Sys.getenv("OUT")),
                           figdir = sprintf("%s/figs_clust1", Sys.getenv("OUT")),
                           SaveEndNamesDir = sprintf("%s/endclusts_clust1", Sys.getenv("OUT")))


