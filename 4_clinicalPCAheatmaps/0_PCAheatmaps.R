# R 4.0
# title:  PCA code for PREDICT paper 2021
# author: Kyle Kimler
# date:   2 Aug 2021
# ---

library(tidyverse)
library(data.table)
library(fdrtool)
library(ggnewscale)
library(ggfortify)

#For the PCA heatmaps, one only needs the per-cell metadata from the Alexandria single-cell-portal study
#This can also be reconstructed with delimiter parsing for major type from Supplemental Table 11
meta <- read_tsv('CD_full_metadata.tsv')

#if they are included, remove doublets and low quality cells
metaCleared <- meta %>% slice(-(grep("(Doub)|(LowQual)", meta$type)))

#Edit Mast cell type from plasma to Mstcl
metaCleared <- metaCleared %>% mutate(type=ifelse(curatedname=='CD.Mstcl.TPSAB1.HPGDS','Mstcl',type))

#Group together all plasma subtypes into IgA,IgM,IgG only.
metaCleared <- metaCleared %>% mutate(curatedname=ifelse(type=='Plsma',gsub("^([^\\.]*\\.[^\\.]*)\\..*", "\\1", curatedname),curatedname))

#Now to reproduce Figure 5, filter only T, Myeloid, and Epithelial cells
meta <- metaCleared %>% filter(type=='Tclls'|type=='Mloid'|type=='Epith')

meta$curatedname = factor(meta$curatedname)
levels(meta$curatedname) = levels(meta$curatedname) %>%
    gsub("(Hs\\.)|(CD\\.)", "", .)

###########################################################
##### Produce the CPM per major celltype table
###########################################################

tbl = as.data.frame(table(factor(meta$curatedname),
                        paste0(meta$patient,"-",meta$antiTNF_response))) %>%
        separate("Var2", c("Var2", "Var3")) %>%
        setNames(c("curatedtype", "patient", "ATR", "count")) %>%
        group_by(patient) %>%
        mutate(CPM = ((count / sum(count)) * 1e6)) %>%
        group_by(curatedtype)

tblwide = tbl %>% dplyr::select(-count) %>% pivot_wider(names_from=curatedtype, values_from= CPM)
tblwide[1:10,]

tbl <- tbl %>% data.frame %>% arrange(desc(CPM))

#

PCmeta <- tblwide %>% data.frame
rownames(PCmeta)=paste0(PCmeta$patient,"-",PCmeta$ATR)
PCmeta <- PCmeta %>% dplyr::select(-patient,-ATR)

pcaCD <- prcomp(PCmeta, scale. = TRUE)

varianceExplained <- pcaCD$sdev^2/sum(pcaCD$sdev^2)

pcaCDloadings <- pcaCD$rotation %>% data.frame %>% rownames_to_column('celltype') %>% arrange(desc(PC1))

#############################################################################
########### Optional plots for variance explained per PC and patient PCA viz (PCs 1 and 2)
#############################################################################

varianceExplained <- varianceExplained %>% data.frame
plotPCVariance <- varianceExplained
plotPCVariance$pc <- paste0('pc_',seq(length(plotPCVariance[,1])))
colnames(plotPCVariance) <- c('exp','pc')
celltype = 'Combo.Tclls.Mloid.Epith'

#var explained per PC
ggplot(data=plotPCVariance,aes(x=reorder(pc, -exp),y=exp)) + geom_col() + ggtitle('Percent variance explained per PC')
ggsave(plot=last_plot(),sprintf('CD_PCA_Percent_variance_per_PC_CelltypesTotalCPM_%s.pdf',celltype))

#PCA viz
autoplot(pcaCD,data=tblwide,colour='patient',shape='ATR')
ggsave(plot=last_plot(),sprintf('CD_PCA_CelltypesTotalCPM_%s.pdf',celltype))

##############################################################################
##############################################################################
# At this point, PC 2 is correlated with patient metadata
# This code is omitted to comply with HIPAA
# We used prcomp on the pcaCDloadings as above.
##############################################################################
##############################################################################
##############################################################################


#############################################################################
########### Loading PCA correlations table
#############################################################################

t <- as.matrix(read_csv('~/Projects/PREDICT/celltypePCA/cortable.csv') %>% distinct %>% select(-Initial_Visit,-Timepoint1,-Timepoint2))
t <- head(t,-1)
t <- t %>% data.frame %>% select(-Patient_ID)

#############################################################################
########### Create heatmap Fig. 5a
#############################################################################

m <- ncol(t)

#cluster correlations and re-order data into clustered order
ord <- hclust(dist(t(t),method='euclidean'),method='ward.D')$order
#methods?
#ord <- hclust(dist(t(t),method='euclidean'))$order
t <- t[,ord]

#Calculate correlations, producing i,j coordinate df
dt <- CJ(i=seq_len(m), j=seq_len(m))[i<j]
dt[, c("p.value"):=(cor.test(t[,i],t[,j],method='spearman')$p.value), by=.(i,j)]
dt[, c("corr"):=(cor(t[,i],t[,j],method='spearman')), by=.(i,j)]

#calculate FDR using fdrtool
dt[,lfdr:=fdrtool::fdrtool(p.value, statistic="pvalue")$lfdr]

#calculate p-value significance stars
dt$stars <- ifelse(dt$p.value < 0.05,ifelse(dt$p.value <0.01,ifelse(dt$p.value<0.001,'***','**'),'*'),'')

textcol <- "grey40"

#double the dataframe to make symmetrical heatmap
dt2 <- dt %>% rename(i=j,j=i)
#add middle row
middles <- data.frame(i=seq(1:length(t)),j=seq(1:length(t)),corr=1,lfdr=0,stars='',rectx=NA,recty=NA,p.value=0,lfdr=0)
dt <- bind_rows(dt,dt2,middles)

#add rectangle coordinates depending on FDR significance
dt$fdrtest <- ifelse(dt$lfdr < 0.05,TRUE,FALSE)
dt$rectx <- ifelse(dt$fdrtest,dt$i,NA)
dt$recty <- ifelse(dt$fdrtest,dt$j,NA)



#plot heatmap
h1 <- ggplot(dt) +
  geom_tile(aes(x=i,y=j, fill=corr,  height=corr,  width=corr)) +
  scale_fill_distiller(palette = 'RdYlBu', direction=-1, limits=c(-1,1),name="Correlation") +
  scale_color_distiller(palette = 'RdYlBu', direction=-1, limits=c(-1,1),name="Correlation") +
  scale_x_continuous('',breaks = seq_len(m),expand=c(0,0), labels = colnames(t),position='top') +
  scale_y_continuous('',breaks = seq_len(m),expand=c(0,0), labels = colnames(t),trans='reverse') +
  theme(axis.text.x=element_text(angle=85, colour=textcol,hjust=-0.1,vjust=0.2),
        panel.background=element_blank(),
        panel.grid.minor=element_blank(),
        panel.grid.major.x=element_blank(),
        panel.grid.major.y=element_blank(),
        #panel.ontop=TRUE,
        #panel.background=element_rect(fill='transparent',colour=textcol),
        axis.text.y=element_text(vjust=0.2,colour=textcol),
        axis.ticks=element_line(size=0.4),
        #panel.border=element_rect(fill='transparent',size=1,colour=textcol),
        plot.margin=margin(0.7,0.4,0.1,0.2,"cm"),
        plot.title=element_text(colour=textcol,hjust=0,size=14,face="bold")) +
  coord_fixed() +
  geom_tile(colour='gray80',size=0.15,fill='transparent',aes(x=i,y=j)) +
  new_scale_color() +
  scale_colour_manual(name='FDR',values='black',labels='< 0.05') +
  geom_rect(size=0.3,fill=NA,aes(colour='black',xmin=rectx-.5,xmax=rectx+.5,ymin=recty-0.5,ymax=recty+0.5)) +
  new_scale_color() +
  scale_colour_manual(name='p-values',values=c('black','grey0','grey1','grey2'),labels=c('< 0.05','< 0.01','<0.001')) +
  geom_text(aes(x=i,y=j+0.125,label=stars,colour='black'),size=2,na.rm=TRUE)

pdf(width=12,height=12,'PCA_clinicalVariables_celltypes_PREDICT3p_Paper.pdf')
print(h1)
dev.off()

#############################################################################
###########  Load pc2 loadings per cell state (Supplemental Table 13)
#############################################################################

pc2 <- read_csv('~/Projects/PREDICT/celltypePCA/pc2cellsubsetsforheatmap.csv') %>% data.frame %>% filter(!is.na(PC2.most.positive),PC2.most.positive!='PC2 most negative')

celltypes <- read_csv('~/Projects/PREDICT/celltypePCA/pc2cellsubsetsforheatmap.csv') %>% data.frame %>% select(PC2.most.positive) %>% filter(!is.na(PC2.most.positive),PC2.most.positive!='PC2 most negative')

#Obtain per patient top and bottom cell frequencies for PC2_EPI_MYE_T
freqs <- read_csv('~/Projects/PREDICT/celltypePCA/CPMtableCD14_perPatient_PREDICT_3p_ILE_LPS.csv')

#reformat subset names
#remove CD. from subset names
freqs <- freqs %>% mutate(subset=sub('...','',freqs$subset))
#change / to .
freqs <- freqs %>% mutate(subset=sub('/','.',freqs$subset))
#change - to .
freqs <- freqs %>% mutate(subset=sub('-','.',freqs$subset))
#change another - to .
freqs <- freqs %>% mutate(subset=sub('-','.',freqs$subset))

#############################################################################
########### Create heatmap Fig. 5a
#############################################################################


#filter celltypes we want to plot
freqs <- freqs %>% filter(subset %in% celltypes[,1])


pxc <- freqs %>% select(-ATR,-count) %>% pivot_wider(names_from=patient,values_from=CPM)
pxc <- pxc %>% column_to_rownames('subset')


#cluster correlations and re-order data into clustered order
ord <- hclust(dist(t(pxc),method='euclidean'),method='complete')$order
#manual order
ord = c(14,8,9,5,6,13,2,1,3,4,7,10,11,12)

pxc <- pxc[,ord]

#Order celltypes by PC2 levels
pxc <- pxc[celltypes$PC2.most.positive,]
#methods?
#ord <- hclust(dist(t(t),method='euclidean'))$order

freqs$i <- match(freqs$subset,rownames(pxc))
freqs$j <- match(freqs$patient,colnames(pxc))

#freqs <- freqs[order(match(freqs$j,ord)),]

freqs <- freqs %>% mutate(ATR=str_replace_all(freqs$ATR,'RESP','FR'))

#Add box for groups
groups <- data.frame(j=unique(freqs$j),i=51,count=0,ATR=count(freqs,patient,ATR)$ATR,CPM=0,patient=unique(freqs$patient),subset='')

dend<-hclust(dist(t(pxc),method='euclidean'),method='complete')

h4 <- ggplot(freqs) +
  geom_tile(aes(x=i,y=j, fill=log(CPM)
    #height=CPM/max(CPM),  width=CPM/max(CPM)
    )) +
  scale_fill_distiller(palette = 'RdGy', direction=-1, limits=c(0,max(log(freqs$CPM))),name=bquote(log ~ "CPM")) +
  scale_color_distiller(palette = 'RdGy', direction=-1, limits=c(0,max(log(freqs$CPM))),name=bquote(log ~ "CPM")) +
  scale_x_continuous('',breaks = seq(1:length(freqs$i %>% unique)),expand=c(0,0), labels = rownames(pxc),position='top') +
  scale_y_continuous('',breaks = seq(1:length(freqs$j %>% unique)),expand=c(0,0), labels = colnames(pxc),trans='reverse') +
  #facet_wrap(factor(ATR,levels=c('NOA','FR','PR')) ~ ., strip.position = "left") +
  theme(axis.text.x=element_text(angle=85, colour=textcol,hjust=-0.1,vjust=0.2,size=6),
      panel.spacing = unit(0, "lines"),
        strip.background = element_blank(),
        strip.placement = "outside",
        panel.background=element_blank(),
        panel.grid.minor=element_blank(),
        panel.grid.major.x=element_blank(),
        panel.grid.major.y=element_blank(),
        aspect.ratio=.1,
        #panel.ontop=TRUE,
        #panel.background=element_rect(fill='transparent',colour=textcol),
        axis.text.y=element_text(vjust=0.2,colour=textcol),
        axis.ticks=element_line(size=0.4),
        #panel.border=element_rect(fill='transparent',size=1,colour=textcol),
        plot.margin=margin(0.7,0.4,0.1,0.2,"cm"),
        plot.title=element_text(colour=textcol,hjust=0,size=14,face="bold")) +
  #coord_fixed() +
  new_scale_color() +
  geom_tile(colour='gray80',size=0.15,fill='transparent',aes(x=i,y=j)) +
  geom_vline(xintercept=25.5, linetype="dashed",
                color = "black", size=.4) +
  facet_grid(fct_relevel(ATR,c('NOA','FR','PR'))~., scales="free_y")


pdf(width=12,height=4.5,'heatmap_PC2_patients_celltypes_PREDICT3p_Paper_facets.pdf')
print(h4)
dev.off()
#library(cowplot)
#library(ggplotify)
#dendro <- as.grob(~plot(dend,hang=-1,cex=0.6))
#plot_grid(dendro,h3,align='v',rel_widths=c(.5,1))



########################################################################################
########### Spearman correlations between cell states: heatmap
########################################################################################


celltypes <- read_csv('~/Projects/PREDICT/celltypePCA/pc2cellsubsetsforheatmap.csv') %>% data.frame %>% select(PC2.most.positive) %>% filter(!is.na(PC2.most.positive),PC2.most.positive!='PC2 most negative')
pc2 <- read_csv('~/Projects/PREDICT/celltypePCA/pc2cellsubsetsforheatmap.csv') %>% data.frame %>% filter(!is.na(PC2.most.positive),PC2.most.positive!='PC2 most negative')
freqs <- read_csv('~/Projects/PREDICT/celltypePCA/CPMtableCD14_perPatient_PREDICT_3p_ILE_LPS.csv')


freqs <- freqs %>% mutate(subset=sub('...','',freqs$subset))
#change / to .
freqs <- freqs %>% mutate(subset=sub('/','.',freqs$subset))
#change - to .
freqs <- freqs %>% mutate(subset=sub('-','.',freqs$subset))
#change another - to .
freqs <- freqs %>% mutate(subset=sub('-','.',freqs$subset))

#filter celltypes we want to plot
freqs <- freqs %>% filter(subset %in% celltypes[,1])


pxc <- freqs %>% select(-ATR,-count) %>% pivot_wider(names_from=patient,values_from=CPM)
pxc <- pxc %>% column_to_rownames('subset')

#cluster correlations and re-order data into clustered order


#methods?
#ord <- hclust(dist(t(t),method='euclidean'))$order
pxc <- pxc[celltypes$PC2.most.positive,]
ord <- hclust(dist(pxc,method='euclidean'),method='ward.D')$order
pxc <- pxc[ord,]

m <- nrow(pxc)

#Calculate correlations, producing i,j coordinate df
dt <- CJ(i=seq_len(m), j=seq_len(m))[i<j]
dt[, c("p.value"):=(cor.test(as.numeric(pxc[i,]),as.numeric(pxc[j,]),method='spearman')$p.value), by=.(i,j)]
dt[, c("corr"):=(cor(as.numeric(pxc[i,]),as.numeric(pxc[j,]),method='spearman')), by=.(i,j)]

#calculate FDR using fdrtool
dt[,lfdr:=fdrtool::fdrtool(p.value, statistic="pvalue")$lfdr]

#calculate p-value significance stars
dt$stars <- ifelse(dt$p.value < 0.05,ifelse(dt$p.value <0.01,ifelse(dt$p.value<0.001,'***','**'),'*'),'')

textcol <- "grey40"

#double the dataframe to make symmetrical heatmap
dt2 <- dt %>% rename(i=j,j=i)
#add middle row
middles <- data.frame(i=seq(1:m),j=seq(1:m),corr=1,lfdr=0,stars='',rectx=NA,recty=NA,p.value=0,lfdr=0)
dt <- bind_rows(dt,dt2,middles)

#add rectangle coordinates depending on FDR significance
dt$fdrtest <- ifelse(dt$lfdr < 0.05,TRUE,FALSE)
dt$rectx <- ifelse(dt$fdrtest,dt$i,NA)
dt$recty <- ifelse(dt$fdrtest,dt$j,NA)

h5 <- ggplot(dt) +
  geom_tile(aes(x=i,y=j, fill=corr,  height=corr,  width=corr)) +
  scale_fill_distiller(palette = 'RdYlBu', direction=-1, limits=c(-1,1),name="Correlation") +
  scale_color_distiller(palette = 'RdYlBu', direction=-1, limits=c(-1,1),name="Correlation") +
  scale_x_continuous('',breaks = seq_len(m),expand=c(0,0), labels = rownames(pxc),position='top') +
  scale_y_continuous('',breaks = seq_len(m),expand=c(0,0), labels = rownames(pxc),trans='reverse') +
  theme(axis.text.x=element_text(angle=85, colour=textcol,hjust=-0.1,vjust=0.2),
        panel.background=element_blank(),
        panel.grid.minor=element_blank(),
        panel.grid.major.x=element_blank(),
        panel.grid.major.y=element_blank(),
        #panel.ontop=TRUE,
        #panel.background=element_rect(fill='transparent',colour=textcol),
        axis.text.y=element_text(vjust=0.2,colour=textcol),
        axis.ticks=element_line(size=0.4),
        #panel.border=element_rect(fill='transparent',size=1,colour=textcol),
        plot.margin=margin(0.7,0.4,0.1,0.2,"cm"),
        plot.title=element_text(colour=textcol,hjust=0,size=14,face="bold")) +
  coord_fixed() +
  geom_tile(colour='gray80',size=0.15,fill='transparent',aes(x=i,y=j)) +
  new_scale_color() +
  scale_colour_manual(name='FDR',values='black',labels='< 0.05') +
  geom_rect(size=0.3,fill=NA,aes(colour='black',xmin=rectx-.5,xmax=rectx+.5,ymin=recty-0.5,ymax=recty+0.5)) +
  new_scale_color() +
  scale_colour_manual(name='p-values',values=c('black','grey0','grey1','grey2'),labels=c('< 0.05','< 0.01','<0.001')) +
  geom_text(aes(x=i,y=j+0.125,label=stars,colour='black'),size=2,na.rm=TRUE)


pdf(width=12,height=12,'heatmap_PC2_celltype_loadings_by_celltypes_Spearman_PREDICT3p_Paper.pdf')
print(h5)
dev.off()

########################################################################################
########### hclus-ordered heatmap with dendrogram correlating the cell states included in the clinical analysis above.
########################################################################################

library(cowplot)
library(ggplotify)
library(ape)
library(rdist)
library(phyloseq)
library(ggdendro)
dend<-hclust(dist(pxc,method='euclidean'),method='complete')
require(dendextend)
gg.dend <- as.ggdend(as.dendrogram(dend))
ddata_x <- dendro_data(as.dendrogram(dend))

plt_dendr <- ggplot(segment(ddata_x)) +
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
  theme_void() +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

sample_names <- rownames(pxc)

# Obtain the dendrogram
#dend <- as.dendrogram(hclust(dist(mat)))
dend_data <- dendro_data(dend)

# Setup the data, so that the layout is inverted (this is more
# "clear" than simply using coord_flip())
segment_data <- with(
    segment(dend_data),
    data.frame(x = x, y = y, xend = xend, yend = yend))
# Use the dendrogram label data to position the gene labels
gene_pos_table <- with(
    dend_data$labels,
    data.frame(x_center = x, sample = as.character(label), height = 1))

# Table to position the samples
sample_pos_table <- data.frame(sample = sample_names) %>%
    mutate(x_center = (1:n()),
           width = 1)

# Neglecting the gap parameters
heatmap_data <- pxc %>% t %>%
    reshape2::melt(value.name = "expr", varnames = c("gene", "sample")) %>%
    left_join(gene_pos_table) %>%
    left_join(sample_pos_table)

# Limits for the vertical axes
gene_axis_limits <- with(
    gene_pos_table,
    c(min(x_center - 0.5 * height), max(x_center + 0.5 * height))
) +
    0.1 * c(-1, 1) # extra spacing: 0.1

# Dendrogram plot
plt_dendr <- ggplot(segment_data) +
    geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
    scale_y_continuous(expand = c(0, 0.5)) +
    scale_x_continuous(breaks = gene_pos_table$x_center,
                       labels = rep("",length(gene_pos_table$sample)),
                       limits = gene_axis_limits,
                       expand = c(0, 0)
                       ) +
    labs(y = "Distance", colour = "", size = "",x="") +
    theme_bw() +
    theme(panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 90))

#Heatmap with labels on lower end

h6 <- ggplot(dt) +
  geom_tile(aes(x=i,y=j, fill=corr,  height=corr,  width=corr)) +
  scale_fill_distiller(palette = 'RdYlBu', direction=-1, limits=c(-1,1),name="Correlation") +
  scale_color_distiller(palette = 'RdYlBu', direction=-1, limits=c(-1,1),name="Correlation") +
  scale_x_continuous('',breaks = seq_len(m),expand=c(0,0), labels = rownames(pxc), position='bottom') +
  scale_y_continuous('',breaks = seq_len(m),expand=c(0,0), labels = rownames(pxc),trans='reverse') +
  theme(axis.text.x=element_text(angle=85, colour=textcol,hjust=-0.1,vjust=0.2),
        panel.background=element_blank(),
        panel.grid.minor=element_blank(),
        panel.grid.major.x=element_blank(),
        panel.grid.major.y=element_blank(),
        #panel.ontop=TRUE,
        #panel.background=element_rect(fill='transparent',colour=textcol),
        axis.text.y=element_text(vjust=0.2,colour=textcol),
        axis.ticks=element_line(size=0.4),
        #panel.border=element_rect(fill='transparent',size=1,colour=textcol),
        plot.margin=margin(0.7,0.4,0.1,0.2,"cm"),
        plot.title=element_text(colour=textcol,hjust=0,size=14,face="bold")) +
  coord_fixed() +
  geom_tile(colour='gray80',size=0.15,fill='transparent',aes(x=i,y=j)) +
  new_scale_color() +
  scale_colour_manual(name='FDR',values='black',labels='< 0.05') +
  geom_rect(size=0.3,fill=NA,aes(colour='black',xmin=rectx-.5,xmax=rectx+.5,ymin=recty-0.5,ymax=recty+0.5)) +
  new_scale_color() +
  scale_colour_manual(name='p-values',values=c('black','grey0','grey1','grey2'),labels=c('< 0.05','< 0.01','<0.001')) +
  geom_text(aes(x=i,y=j+0.125,label=stars,colour='black'),size=2,na.rm=TRUE)


library(cowplot)
pdf(height=20,width=12,'~/Projects/PREDICT/celltypePCA/testdend.pdf')
plot_grid(plt_dendr, h6, align = 'v', nrow=2, ncol=1, rel_heights=c(.1,1), axis='lr', vjust=c(6,1.5))
dev.off()


########################################################################################
########### Pearson correlations between cell states: heatmap
########################################################################################



pxc <- freqs %>% select(-ATR,-count) %>% pivot_wider(names_from=patient,values_from=CPM)
pxc <- pxc %>% column_to_rownames('subset')

ord <- hclust(dist(pxc,method='euclidean'),method='ward.D')$order
#methods?
#ord <- hclust(dist(t(t),method='euclidean'))$order
pxc <- pxc[celltypes$PC2.most.positive,]

#optionally log transform CPM to perform Pearson correlation
pxc <- log(1+pxc)

pxc <- pxc[ord,]

#Calculate correlations, producing i,j coordinate df
dt <- CJ(i=seq_len(m), j=seq_len(m))[i<j]
dt[, c("p.value"):=(cor.test(as.numeric(pxc[i,]),as.numeric(pxc[j,]),method='pearson')$p.value), by=.(i,j)]
dt[, c("corr"):=(cor(as.numeric(pxc[i,]),as.numeric(pxc[j,]),method='pearson')), by=.(i,j)]

#calculate FDR using fdrtool
dt[,lfdr:=fdrtool::fdrtool(p.value, statistic="pvalue")$lfdr]

#calculate p-value significance stars
dt$stars <- ifelse(dt$p.value < 0.05,ifelse(dt$p.value <0.01,ifelse(dt$p.value<0.001,'***','**'),'*'),'')

textcol <- "grey40"

#double the dataframe to make symmetrical heatmap
dt2 <- dt %>% rename(i=j,j=i)
#add middle row
middles <- data.frame(i=seq(1:m),j=seq(1:m),corr=1,lfdr=0,stars='',rectx=NA,recty=NA,p.value=0,lfdr=0)
dt <- bind_rows(dt,dt2,middles)

#add rectangle coordinates depending on FDR significance
dt$fdrtest <- ifelse(dt$lfdr < 0.05,TRUE,FALSE)
dt$rectx <- ifelse(dt$fdrtest,dt$i,NA)
dt$recty <- ifelse(dt$fdrtest,dt$j,NA)

m <- nrow(pxc)


h7 <- ggplot(dt) +
  geom_tile(aes(x=i,y=j, fill=corr,  height=corr,  width=corr)) +
  scale_fill_distiller(palette = 'RdYlBu', direction=-1, limits=c(-1,1),name="Correlation") +
  scale_color_distiller(palette = 'RdYlBu', direction=-1, limits=c(-1,1),name="Correlation") +
  scale_x_continuous('',breaks = seq_len(m),expand=c(0,0), labels = rownames(pxc),position='top') +
  scale_y_continuous('',breaks = seq_len(m),expand=c(0,0), labels = rownames(pxc),trans='reverse') +
  theme(axis.text.x=element_text(angle=85, colour=textcol,hjust=-0.1,vjust=0.2),
        panel.background=element_blank(),
        panel.grid.minor=element_blank(),
        panel.grid.major.x=element_blank(),
        panel.grid.major.y=element_blank(),
        #panel.ontop=TRUE,
        #panel.background=element_rect(fill='transparent',colour=textcol),
        axis.text.y=element_text(vjust=0.2,colour=textcol),
        axis.ticks=element_line(size=0.4),
        #panel.border=element_rect(fill='transparent',size=1,colour=textcol),
        plot.margin=margin(0.7,0.4,0.1,0.2,"cm"),
        plot.title=element_text(colour=textcol,hjust=0,size=14,face="bold")) +
  coord_fixed() +
  geom_tile(colour='gray80',size=0.15,fill='transparent',aes(x=i,y=j)) +
  new_scale_color() +
  scale_colour_manual(name='FDR',values='black',labels='< 0.05') +
  geom_rect(size=0.3,fill=NA,aes(colour='black',xmin=rectx-.5,xmax=rectx+.5,ymin=recty-0.5,ymax=recty+0.5)) +
  new_scale_color() +
  scale_colour_manual(name='p-values',values=c('black','grey0','grey1','grey2'),labels=c('< 0.05','< 0.01','<0.001')) +
  geom_text(aes(x=i,y=j+0.125,label=stars,colour='black'),size=2,na.rm=TRUE)

pdf(width=12,height=12,'heatmap_log1pPC2_celltype_loadings_by_celltypes_Pearson_PREDICT3p_Paper.pdf')
print(h7)
dev.off()
