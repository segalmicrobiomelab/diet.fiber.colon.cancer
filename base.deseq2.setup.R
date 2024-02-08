### DESeq2 Work 
# DESEq2 Calculations 

library(phyloseq)
library(ade4)
library(vegan)
library(phyloseq)
library(ade4)
library(vegan)
library(biomformat)
library(devtools)
library(readr)
library(readtext)
library(qiime2R)
library("ape")
library("ggpubr")
library("decontam")
library("ggplot2")
library("phyloseq")
library("reshape2")
library("tidyr")
library("matrixStats")
library(DESeq2)
library(edgeR)
library(limma)
library(Glimma)
library(gplots)
library(RColorBrewer)
library(pheatmap)
library(ggplot2)
library(ggrepel)
library(scales)
library(data.table)
library(fBasics)
library(forcats)
library(omu)
# library(maptools)
library(phyloseq)
library("vegan")
library(ade4)
library("reshape2")
library(dplyr)
library(gridExtra)
library("decontam")
library("ggplot2")
library("phyloseq")
library(gridExtra)
library(gtable)

physeq
subset.cecum.physeq = subset_samples(physeq, Source  %in% c("cecum"))
subset.fecal.physeq = subset_samples(physeq, Source  %in% c("fecal"))

# Subset as needed
subset.cecum.physeq.HF.v.NC = subset_samples(subset.cecum.physeq, Treatment  %in% c("HF", "NC"))
subset.cecum.physeq.LF.v.NC = subset_samples(subset.cecum.physeq, Treatment  %in% c("LF", "NC"))
subset.cecum.physeq.HF.v.LF = subset_samples(subset.cecum.physeq, Treatment  %in% c("HF", "LF"))
subset.fecal.physeq.HF.v.NC = subset_samples(subset.fecal.physeq, Treatment  %in% c("HF", "NC"))
subset.fecal.physeq.LF.v.NC = subset_samples(subset.fecal.physeq, Treatment  %in% c("LF", "NC"))
subset.fecal.physeq.HF.v.LF = subset_samples(subset.fecal.physeq, Treatment  %in% c("HF", "LF"))

subset.cecum.physeq
subset.cecum.physeq

sample_data(subset.cecum.physeq)
tax_table(subset.cecum.physeq)

sample_data(subset.fecal.physeq)
tax_table(subset.fecal.physeq)

#Set Theme For Figures
theme <-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black"),
              axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"),
              legend.position="none")

# Set your alpha level - can be 0.05, 0.1, or 0.2 as it is here
alpha = 0.2

sample_data(subset.cecum.physeq)$Treatment
sample_data(subset.fecal.physeq)$Treatment


# Subset as needed
rownames(sample_data(subset.cecum.physeq.HF.v.NC))
rownames(sample_data(subset.cecum.physeq.LF.v.NC))
rownames(sample_data(subset.cecum.physeq.HF.v.LF))
rownames(sample_data(subset.fecal.physeq.HF.v.NC))
rownames(sample_data(subset.fecal.physeq.LF.v.NC))
rownames(sample_data(subset.fecal.physeq.HF.v.LF))

subset.genus.cecum.physeq.HF.v.NC = tax_glom(subset.cecum.physeq.HF.v.NC, taxrank = "Genus")
subset.genus.cecum.physeq.LF.v.NC = tax_glom(subset.cecum.physeq.LF.v.NC, taxrank = "Genus")
subset.genus.cecum.physeq.HF.v.LF = tax_glom(subset.cecum.physeq.HF.v.LF, taxrank = "Genus")
subset.genus.fecal.physeq.HF.v.NC = tax_glom(subset.fecal.physeq.HF.v.NC, taxrank = "Genus")
subset.genus.fecal.physeq.LF.v.NC = tax_glom(subset.fecal.physeq.LF.v.NC, taxrank = "Genus")
subset.genus.fecal.physeq.HF.v.LF = tax_glom(subset.fecal.physeq.HF.v.LF, taxrank = "Genus")

# Change filter strategy
filtered.genus.cecum.physeq.HF.v.NC = genefilter_sample(subset.genus.cecum.physeq.HF.v.NC, filterfun_sample(function(x) x > 0.001), A = 0.001 * nsamples(subset.genus.cecum.physeq.HF.v.NC))
pruned.genus.cecum.physeq.HF.v.NC = prune_taxa(filtered.genus.cecum.physeq.HF.v.NC, subset.genus.cecum.physeq.HF.v.NC)
ntaxa(subset.genus.cecum.physeq.HF.v.NC)
ntaxa(pruned.genus.cecum.physeq.HF.v.NC)

filtered.genus.cecum.physeq.LF.v.NC = genefilter_sample(subset.genus.cecum.physeq.LF.v.NC, filterfun_sample(function(x) x > 0.001), A = 0.001 * nsamples(subset.genus.cecum.physeq.LF.v.NC))
pruned.genus.cecum.physeq.LF.v.NC = prune_taxa(filtered.genus.cecum.physeq.LF.v.NC, subset.genus.cecum.physeq.LF.v.NC)
ntaxa(subset.genus.cecum.physeq.LF.v.NC)
ntaxa(pruned.genus.cecum.physeq.LF.v.NC)

filtered.genus.cecum.physeq.HF.v.LF = genefilter_sample(subset.genus.cecum.physeq.HF.v.LF, filterfun_sample(function(x) x > 0.001), A = 0.001 * nsamples(subset.genus.cecum.physeq.HF.v.LF))
pruned.genus.cecum.physeq.HF.v.LF = prune_taxa(filtered.genus.cecum.physeq.HF.v.LF, subset.genus.cecum.physeq.HF.v.LF)
ntaxa(subset.genus.cecum.physeq.HF.v.LF)
ntaxa(pruned.genus.cecum.physeq.HF.v.LF)

filtered.genus.fecal.physeq.HF.v.NC = genefilter_sample(subset.genus.fecal.physeq.HF.v.NC, filterfun_sample(function(x) x > 0.001), A = 0.001 * nsamples(subset.genus.fecal.physeq.HF.v.NC))
pruned.genus.fecal.physeq.HF.v.NC = prune_taxa(filtered.genus.fecal.physeq.HF.v.NC, subset.genus.fecal.physeq.HF.v.NC)
ntaxa(subset.genus.fecal.physeq.HF.v.NC)
ntaxa(pruned.genus.fecal.physeq.HF.v.NC)

filtered.genus.fecal.physeq.LF.v.NC = genefilter_sample(subset.genus.fecal.physeq.LF.v.NC, filterfun_sample(function(x) x > 0.001), A = 0.001 * nsamples(subset.genus.fecal.physeq.LF.v.NC))
pruned.genus.fecal.physeq.LF.v.NC = prune_taxa(filtered.genus.fecal.physeq.LF.v.NC, subset.genus.fecal.physeq.LF.v.NC)
ntaxa(subset.genus.fecal.physeq.LF.v.NC)
ntaxa(pruned.genus.fecal.physeq.LF.v.NC)

filtered.genus.fecal.physeq.HF.v.LF = genefilter_sample(subset.genus.fecal.physeq.HF.v.LF, filterfun_sample(function(x) x > 0.001), A = 0.001 * nsamples(subset.genus.fecal.physeq.HF.v.LF))
pruned.genus.fecal.physeq.HF.v.LF = prune_taxa(filtered.genus.fecal.physeq.HF.v.LF, subset.genus.fecal.physeq.HF.v.LF)
ntaxa(subset.genus.fecal.physeq.HF.v.LF)
ntaxa(pruned.genus.fecal.physeq.HF.v.LF)

# Prune the samples 
pruned.genus.cecum.physeq.HF.v.NC
pruned.genus.cecum.physeq.LF.v.NC
pruned.genus.cecum.physeq.HF.v.LF
pruned.genus.fecal.physeq.HF.v.NC
pruned.genus.fecal.physeq.LF.v.NC
pruned.genus.fecal.physeq.HF.v.LF

# Normalizing the sample here
normalizeSample = function(x) {
  x/sum(x)
}

pruned.genus.rel.table.cecum.physeq.HF.v.NC = transformSampleCounts(pruned.genus.cecum.physeq.HF.v.NC, normalizeSample)
pruned.genus.rel.table.cecum.physeq.LF.v.NC = transformSampleCounts(pruned.genus.cecum.physeq.LF.v.NC, normalizeSample)
pruned.genus.rel.table.cecum.physeq.HF.v.LF = transformSampleCounts(pruned.genus.cecum.physeq.HF.v.LF, normalizeSample)
pruned.genus.rel.table.fecal.physeq.HF.v.NC = transformSampleCounts(pruned.genus.fecal.physeq.HF.v.NC, normalizeSample)
pruned.genus.rel.table.fecal.physeq.LF.v.NC = transformSampleCounts(pruned.genus.fecal.physeq.LF.v.NC, normalizeSample)
pruned.genus.rel.table.fecal.physeq.HF.v.LF = transformSampleCounts(pruned.genus.fecal.physeq.HF.v.LF, normalizeSample)

# Change variable of interest

# Change the phyloseq object to DESeq2 Object
diagdds.pruned.genus.cecum.physeq.HF.v.NC <- phyloseq_to_deseq2(pruned.genus.cecum.physeq.HF.v.NC, ~ Treatment)
diagdds.pruned.genus.cecum.physeq.LF.v.NC <- phyloseq_to_deseq2(pruned.genus.cecum.physeq.LF.v.NC, ~ Treatment)
diagdds.pruned.genus.cecum.physeq.HF.v.LF <- phyloseq_to_deseq2(pruned.genus.cecum.physeq.HF.v.LF, ~ Treatment)
diagdds.pruned.genus.fecal.physeq.HF.v.NC <- phyloseq_to_deseq2(pruned.genus.fecal.physeq.HF.v.NC, ~ Treatment)
diagdds.pruned.genus.fecal.physeq.LF.v.NC <- phyloseq_to_deseq2(pruned.genus.fecal.physeq.LF.v.NC, ~ Treatment)
diagdds.pruned.genus.fecal.physeq.HF.v.LF <- phyloseq_to_deseq2(pruned.genus.fecal.physeq.HF.v.LF, ~ Treatment)

otu_table(subset.cecum.physeq.HF.v.NC)

gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}

geoMeans.diagdds.pruned.genus.cecum.physeq.HF.v.NC = apply(counts(diagdds.pruned.genus.cecum.physeq.HF.v.NC), 1, gm_mean)
geoMeans.diagdds.pruned.genus.cecum.physeq.LF.v.NC = apply(counts(diagdds.pruned.genus.cecum.physeq.LF.v.NC), 1, gm_mean)
geoMeans.diagdds.pruned.genus.cecum.physeq.HF.v.LF = apply(counts(diagdds.pruned.genus.cecum.physeq.HF.v.LF), 1, gm_mean)
geoMeans.diagdds.pruned.genus.fecal.physeq.HF.v.NC = apply(counts(diagdds.pruned.genus.fecal.physeq.HF.v.NC), 1, gm_mean)
geoMeans.diagdds.pruned.genus.fecal.physeq.LF.v.NC = apply(counts(diagdds.pruned.genus.fecal.physeq.LF.v.NC), 1, gm_mean)
geoMeans.diagdds.pruned.genus.fecal.physeq.HF.v.LF = apply(counts(diagdds.pruned.genus.fecal.physeq.HF.v.LF), 1, gm_mean)

diagdds.pruned.genus.cecum.physeq.HF.v.NC = estimateSizeFactors(diagdds.pruned.genus.cecum.physeq.HF.v.NC, geoMeans = geoMeans.diagdds.pruned.genus.cecum.physeq.HF.v.NC)
diagdds.pruned.genus.cecum.physeq.HF.v.NC = estimateDispersions(diagdds.pruned.genus.cecum.physeq.HF.v.NC)

diagdds.pruned.genus.cecum.physeq.LF.v.NC = estimateSizeFactors(diagdds.pruned.genus.cecum.physeq.LF.v.NC, geoMeans = geoMeans.diagdds.pruned.genus.cecum.physeq.LF.v.NC)
diagdds.pruned.genus.cecum.physeq.LF.v.NC = estimateDispersions(diagdds.pruned.genus.cecum.physeq.LF.v.NC)

diagdds.pruned.genus.cecum.physeq.HF.v.LF = estimateSizeFactors(diagdds.pruned.genus.cecum.physeq.HF.v.LF, geoMeans = geoMeans.diagdds.pruned.genus.cecum.physeq.HF.v.LF)
diagdds.pruned.genus.cecum.physeq.HF.v.LF = estimateDispersions(diagdds.pruned.genus.cecum.physeq.HF.v.LF)

diagdds.pruned.genus.fecal.physeq.HF.v.NC = estimateSizeFactors(diagdds.pruned.genus.fecal.physeq.HF.v.NC, geoMeans = geoMeans.diagdds.pruned.genus.fecal.physeq.HF.v.NC)
diagdds.pruned.genus.fecal.physeq.HF.v.NC = estimateDispersions(diagdds.pruned.genus.fecal.physeq.HF.v.NC)

diagdds.pruned.genus.fecal.physeq.LF.v.NC = estimateSizeFactors(diagdds.pruned.genus.fecal.physeq.LF.v.NC, geoMeans = geoMeans.diagdds.pruned.genus.fecal.physeq.LF.v.NC)
diagdds.pruned.genus.fecal.physeq.LF.v.NC = estimateDispersions(diagdds.pruned.genus.fecal.physeq.LF.v.NC)

diagdds.pruned.genus.fecal.physeq.HF.v.LF = estimateSizeFactors(diagdds.pruned.genus.fecal.physeq.HF.v.LF, geoMeans = geoMeans.diagdds.pruned.genus.fecal.physeq.HF.v.LF)
diagdds.pruned.genus.fecal.physeq.HF.v.LF = estimateDispersions(diagdds.pruned.genus.fecal.physeq.HF.v.LF)

# change variable of interest
diagdds.pruned.genus.cecum.physeq.HF.v.NC$Treatment <- droplevels(diagdds.pruned.genus.cecum.physeq.HF.v.NC$Treatment)
diagdds.pruned.genus.cecum.physeq.LF.v.NC$Treatment <- droplevels(diagdds.pruned.genus.cecum.physeq.LF.v.NC$Treatment)
diagdds.pruned.genus.cecum.physeq.HF.v.LF$Treatment <- droplevels(diagdds.pruned.genus.cecum.physeq.HF.v.LF$Treatment)
diagdds.pruned.genus.fecal.physeq.HF.v.NC$Treatment <- droplevels(diagdds.pruned.genus.fecal.physeq.HF.v.NC$Treatment)
diagdds.pruned.genus.fecal.physeq.LF.v.NC$Treatment <- droplevels(diagdds.pruned.genus.fecal.physeq.LF.v.NC$Treatment)
diagdds.pruned.genus.cecum.physeq.HF.v.NC$Treatment <- droplevels(diagdds.pruned.genus.cecum.physeq.HF.v.NC$Treatment)

# choose variable of interest and reference
diagdds.pruned.genus.cecum.physeq.HF.v.NC$Treatment <- relevel(diagdds.pruned.genus.cecum.physeq.HF.v.NC$Treatment, ref ="NC")
diagdds.pruned.genus.cecum.physeq.LF.v.NC$Treatment <- relevel(diagdds.pruned.genus.cecum.physeq.LF.v.NC$Treatment, ref ="NC")
diagdds.pruned.genus.cecum.physeq.HF.v.LF$Treatment <- relevel(diagdds.pruned.genus.cecum.physeq.HF.v.LF$Treatment, ref ="LF")
diagdds.pruned.genus.fecal.physeq.HF.v.NC$Treatment <- relevel(diagdds.pruned.genus.fecal.physeq.HF.v.NC$Treatment, ref ="NC")
diagdds.pruned.genus.fecal.physeq.LF.v.NC$Treatment <- relevel(diagdds.pruned.genus.fecal.physeq.LF.v.NC$Treatment, ref ="NC")
diagdds.pruned.genus.fecal.physeq.HF.v.LF$Treatment <- relevel(diagdds.pruned.genus.fecal.physeq.HF.v.LF$Treatment, ref ="LF")

diagdds.pruned.genus.cecum.physeq.HF.v.NC <- DESeq(diagdds.pruned.genus.cecum.physeq.HF.v.NC)
diagdds.pruned.genus.cecum.physeq.LF.v.NC <- DESeq(diagdds.pruned.genus.cecum.physeq.LF.v.NC)
diagdds.pruned.genus.cecum.physeq.HF.v.LF <- DESeq(diagdds.pruned.genus.cecum.physeq.HF.v.LF)
diagdds.pruned.genus.fecal.physeq.HF.v.NC <- DESeq(diagdds.pruned.genus.fecal.physeq.HF.v.NC)
diagdds.pruned.genus.fecal.physeq.LF.v.NC <- DESeq(diagdds.pruned.genus.fecal.physeq.LF.v.NC)
diagdds.pruned.genus.fecal.physeq.HF.v.LF <- DESeq(diagdds.pruned.genus.fecal.physeq.HF.v.LF)

res.diagdds.pruned.genus.cecum.physeq.HF.v.NC <- results(diagdds.pruned.genus.cecum.physeq.HF.v.NC)
res.diagdds.pruned.genus.cecum.physeq.LF.v.NC <- results(diagdds.pruned.genus.cecum.physeq.LF.v.NC)
res.diagdds.pruned.genus.cecum.physeq.HF.v.LF <- results(diagdds.pruned.genus.cecum.physeq.HF.v.LF)
res.diagdds.pruned.genus.fecal.physeq.HF.v.NC <- results(diagdds.pruned.genus.fecal.physeq.HF.v.NC)
res.diagdds.pruned.genus.fecal.physeq.LF.v.NC <- results(diagdds.pruned.genus.fecal.physeq.LF.v.NC)
res.diagdds.pruned.genus.fecal.physeq.HF.v.LF <- results(diagdds.pruned.genus.fecal.physeq.HF.v.LF)

res.diagdds.pruned.genus.cecum.physeq.HF.v.NC = res.diagdds.pruned.genus.cecum.physeq.HF.v.NC[order(res.diagdds.pruned.genus.cecum.physeq.HF.v.NC$padj, na.last = NA), ]
res.diagdds.pruned.genus.cecum.physeq.LF.v.NC = res.diagdds.pruned.genus.cecum.physeq.LF.v.NC[order(res.diagdds.pruned.genus.cecum.physeq.LF.v.NC$padj, na.last = NA), ]
res.diagdds.pruned.genus.cecum.physeq.HF.v.LF = res.diagdds.pruned.genus.cecum.physeq.HF.v.LF[order(res.diagdds.pruned.genus.cecum.physeq.HF.v.LF$padj, na.last = NA), ]
res.diagdds.pruned.genus.fecal.physeq.HF.v.NC = res.diagdds.pruned.genus.fecal.physeq.HF.v.NC[order(res.diagdds.pruned.genus.fecal.physeq.HF.v.NC$padj, na.last = NA), ]
res.diagdds.pruned.genus.fecal.physeq.LF.v.NC = res.diagdds.pruned.genus.fecal.physeq.LF.v.NC[order(res.diagdds.pruned.genus.fecal.physeq.LF.v.NC$padj, na.last = NA), ]
res.diagdds.pruned.genus.fecal.physeq.HF.v.LF = res.diagdds.pruned.genus.fecal.physeq.HF.v.LF[order(res.diagdds.pruned.genus.fecal.physeq.HF.v.LF$padj, na.last = NA), ]

select_genes.res.diagdds.pruned.genus.cecum.physeq.HF.v.NC = rownames(res.diagdds.pruned.genus.cecum.physeq.HF.v.NC[res.diagdds.pruned.genus.cecum.physeq.HF.v.NC$padj < alpha & !is.na(res.diagdds.pruned.genus.cecum.physeq.HF.v.NC$padj), ])[1:50]
select_genes.res.diagdds.pruned.genus.cecum.physeq.LF.v.NC = rownames(res.diagdds.pruned.genus.cecum.physeq.LF.v.NC[res.diagdds.pruned.genus.cecum.physeq.LF.v.NC$padj < alpha & !is.na(res.diagdds.pruned.genus.cecum.physeq.LF.v.NC$padj), ])[1:50]
select_genes.res.diagdds.pruned.genus.cecum.physeq.HF.v.LF = rownames(res.diagdds.pruned.genus.cecum.physeq.HF.v.LF[res.diagdds.pruned.genus.cecum.physeq.HF.v.LF$padj < alpha & !is.na(res.diagdds.pruned.genus.cecum.physeq.HF.v.LF$padj), ])[1:50]
select_genes.res.diagdds.pruned.genus.fecal.physeq.HF.v.NC = rownames(res.diagdds.pruned.genus.fecal.physeq.HF.v.NC[res.diagdds.pruned.genus.fecal.physeq.HF.v.NC$padj < alpha & !is.na(res.diagdds.pruned.genus.fecal.physeq.HF.v.NC$padj), ])[1:50]
select_genes.res.diagdds.pruned.genus.fecal.physeq.LF.v.NC = rownames(res.diagdds.pruned.genus.fecal.physeq.LF.v.NC[res.diagdds.pruned.genus.fecal.physeq.LF.v.NC$padj < alpha & !is.na(res.diagdds.pruned.genus.fecal.physeq.LF.v.NC$padj), ])[1:50]
select_genes.res.diagdds.pruned.genus.fecal.physeq.HF.v.LF = rownames(res.diagdds.pruned.genus.fecal.physeq.HF.v.LF[res.diagdds.pruned.genus.fecal.physeq.HF.v.LF$padj < alpha & !is.na(res.diagdds.pruned.genus.fecal.physeq.HF.v.LF$padj), ])[1:50]

# from above
pruned.genus.cecum.physeq.HF.v.NC
pruned.genus.cecum.physeq.LF.v.NC
pruned.genus.cecum.physeq.HF.v.LF
pruned.genus.fecal.physeq.HF.v.NC
pruned.genus.fecal.physeq.LF.v.NC
pruned.genus.fecal.physeq.HF.v.LF

pruned.genus.rel.table.cecum.physeq.HF.v.NC 
pruned.genus.rel.table.cecum.physeq.LF.v.NC 
pruned.genus.rel.table.cecum.physeq.HF.v.LF
pruned.genus.rel.table.fecal.physeq.HF.v.NC
pruned.genus.rel.table.fecal.physeq.LF.v.NC
pruned.genus.rel.table.fecal.physeq.HF.v.LF

res.diagdds.pruned.genus.cecum.physeq.HF.v.NC
res.diagdds.pruned.genus.cecum.physeq.LF.v.NC
res.diagdds.pruned.genus.cecum.physeq.HF.v.LF
res.diagdds.pruned.genus.fecal.physeq.HF.v.NC
res.diagdds.pruned.genus.fecal.physeq.LF.v.NC
res.diagdds.pruned.genus.fecal.physeq.HF.v.LF

subset.genus.cecum.physeq.HF.v.NC
subset.genus.cecum.physeq.LF.v.NC
subset.genus.cecum.physeq.HF.v.LF
subset.genus.fecal.physeq.HF.v.NC
subset.genus.fecal.physeq.LF.v.NC
subset.genus.fecal.physeq.HF.v.LF

pruned.genus.rel.table.cecum.physeq.HF.v.NC = transformSampleCounts(pruned.genus.cecum.physeq.HF.v.NC, normalizeSample)
pruned.genus.rel.table.cecum.physeq.LF.v.NC = transformSampleCounts(pruned.genus.cecum.physeq.LF.v.NC, normalizeSample)
pruned.genus.rel.table.cecum.physeq.HF.v.LF = transformSampleCounts(pruned.genus.cecum.physeq.HF.v.LF, normalizeSample)
pruned.genus.rel.table.fecal.physeq.HF.v.NC = transformSampleCounts(pruned.genus.fecal.physeq.HF.v.NC, normalizeSample)
pruned.genus.rel.table.fecal.physeq.LF.v.NC = transformSampleCounts(pruned.genus.fecal.physeq.LF.v.NC, normalizeSample)
pruned.genus.rel.table.fecal.physeq.HF.v.LF = transformSampleCounts(pruned.genus.fecal.physeq.HF.v.LF, normalizeSample)

res.diagdds.pruned.genus.cecum.physeq.HF.v.NC <- merge(as(res.diagdds.pruned.genus.cecum.physeq.HF.v.NC, "data.frame"), as(tax_table(pruned.genus.rel.table.cecum.physeq.HF.v.NC), "matrix"), by.x=0, by.y=0)
res.diagdds.pruned.genus.cecum.physeq.LF.v.NC <- merge(as(res.diagdds.pruned.genus.cecum.physeq.LF.v.NC, "data.frame"), as(tax_table(pruned.genus.rel.table.cecum.physeq.LF.v.NC), "matrix"), by.x=0, by.y=0)
res.diagdds.pruned.genus.cecum.physeq.HF.v.LF <- merge(as(res.diagdds.pruned.genus.cecum.physeq.HF.v.LF, "data.frame"), as(tax_table(pruned.genus.rel.table.cecum.physeq.HF.v.LF), "matrix"), by.x=0, by.y=0)
res.diagdds.pruned.genus.fecal.physeq.HF.v.NC <- merge(as(res.diagdds.pruned.genus.fecal.physeq.HF.v.NC, "data.frame"), as(tax_table(pruned.genus.rel.table.fecal.physeq.HF.v.NC), "matrix"), by.x=0, by.y=0)
res.diagdds.pruned.genus.fecal.physeq.LF.v.NC <- merge(as(res.diagdds.pruned.genus.fecal.physeq.LF.v.NC, "data.frame"), as(tax_table(pruned.genus.rel.table.fecal.physeq.LF.v.NC), "matrix"), by.x=0, by.y=0)
res.diagdds.pruned.genus.fecal.physeq.HF.v.LF <- merge(as(res.diagdds.pruned.genus.fecal.physeq.HF.v.LF, "data.frame"), as(tax_table(pruned.genus.rel.table.fecal.physeq.HF.v.LF), "matrix"), by.x=0, by.y=0)

# res$row2 <- paste(res$Kingdom,res$Phylum,res$Class,res$Order,res$Family,res$Genus,res$Species,res$Row.names)
res.diagdds.pruned.genus.cecum.physeq.HF.v.NC$row2 <- paste(res.diagdds.pruned.genus.cecum.physeq.HF.v.NC$Kingdom,res.diagdds.pruned.genus.cecum.physeq.HF.v.NC$Phylum,res.diagdds.pruned.genus.cecum.physeq.HF.v.NC$Class,res.diagdds.pruned.genus.cecum.physeq.HF.v.NC$Order,res.diagdds.pruned.genus.cecum.physeq.HF.v.NC$Family,res.diagdds.pruned.genus.cecum.physeq.HF.v.NC$Genus,res.diagdds.pruned.genus.cecum.physeq.HF.v.NC$Species,res.diagdds.pruned.genus.cecum.physeq.HF.v.NC$Row.names)
res.diagdds.pruned.genus.cecum.physeq.LF.v.NC$row2 <- paste(res.diagdds.pruned.genus.cecum.physeq.LF.v.NC$Kingdom,res.diagdds.pruned.genus.cecum.physeq.LF.v.NC$Phylum,res.diagdds.pruned.genus.cecum.physeq.LF.v.NC$Class,res.diagdds.pruned.genus.cecum.physeq.LF.v.NC$Order,res.diagdds.pruned.genus.cecum.physeq.LF.v.NC$Family,res.diagdds.pruned.genus.cecum.physeq.LF.v.NC$Genus,res.diagdds.pruned.genus.cecum.physeq.LF.v.NC$Species,res.diagdds.pruned.genus.cecum.physeq.LF.v.NC$Row.names)
res.diagdds.pruned.genus.cecum.physeq.HF.v.LF$row2 <- paste(res.diagdds.pruned.genus.cecum.physeq.HF.v.LF$Kingdom,res.diagdds.pruned.genus.cecum.physeq.HF.v.LF$Phylum,res.diagdds.pruned.genus.cecum.physeq.HF.v.LF$Class,res.diagdds.pruned.genus.cecum.physeq.HF.v.LF$Order,res.diagdds.pruned.genus.cecum.physeq.HF.v.LF$Family,res.diagdds.pruned.genus.cecum.physeq.HF.v.LF$Genus,res.diagdds.pruned.genus.cecum.physeq.HF.v.LF$Species,res.diagdds.pruned.genus.cecum.physeq.HF.v.LF$Row.names)
res.diagdds.pruned.genus.fecal.physeq.HF.v.NC$row2 <- paste(res.diagdds.pruned.genus.fecal.physeq.HF.v.NC$Kingdom,res.diagdds.pruned.genus.fecal.physeq.HF.v.NC$Phylum,res.diagdds.pruned.genus.fecal.physeq.HF.v.NC$Class,res.diagdds.pruned.genus.fecal.physeq.HF.v.NC$Order,res.diagdds.pruned.genus.fecal.physeq.HF.v.NC$Family,res.diagdds.pruned.genus.fecal.physeq.HF.v.NC$Genus,res.diagdds.pruned.genus.fecal.physeq.HF.v.NC$Species,res.diagdds.pruned.genus.fecal.physeq.HF.v.NC$Row.names)
res.diagdds.pruned.genus.fecal.physeq.LF.v.NC$row2 <- paste(res.diagdds.pruned.genus.fecal.physeq.LF.v.NC$Kingdom,res.diagdds.pruned.genus.fecal.physeq.LF.v.NC$Phylum,res.diagdds.pruned.genus.fecal.physeq.LF.v.NC$Class,res.diagdds.pruned.genus.fecal.physeq.LF.v.NC$Order,res.diagdds.pruned.genus.fecal.physeq.LF.v.NC$Family,res.diagdds.pruned.genus.fecal.physeq.LF.v.NC$Genus,res.diagdds.pruned.genus.fecal.physeq.LF.v.NC$Species,res.diagdds.pruned.genus.fecal.physeq.LF.v.NC$Row.names)
res.diagdds.pruned.genus.fecal.physeq.HF.v.LF$row2 <- paste(res.diagdds.pruned.genus.fecal.physeq.HF.v.LF$Kingdom,res.diagdds.pruned.genus.fecal.physeq.HF.v.LF$Phylum,res.diagdds.pruned.genus.fecal.physeq.HF.v.LF$Class,res.diagdds.pruned.genus.fecal.physeq.HF.v.LF$Order,res.diagdds.pruned.genus.fecal.physeq.HF.v.LF$Family,res.diagdds.pruned.genus.fecal.physeq.HF.v.LF$Genus,res.diagdds.pruned.genus.fecal.physeq.HF.v.LF$Species,res.diagdds.pruned.genus.fecal.physeq.HF.v.LF$Row.names)

# res$row2 <- gsub('\\s+', '|', res$row2)
res.diagdds.pruned.genus.cecum.physeq.HF.v.NC$row2 <- gsub('\\s+', '|', res.diagdds.pruned.genus.cecum.physeq.HF.v.NC$row2)
res.diagdds.pruned.genus.cecum.physeq.LF.v.NC$row2 <- gsub('\\s+', '|', res.diagdds.pruned.genus.cecum.physeq.LF.v.NC$row2)
res.diagdds.pruned.genus.cecum.physeq.HF.v.LF$row2 <- gsub('\\s+', '|', res.diagdds.pruned.genus.cecum.physeq.HF.v.LF$row2)
res.diagdds.pruned.genus.fecal.physeq.HF.v.NC$row2 <- gsub('\\s+', '|', res.diagdds.pruned.genus.fecal.physeq.HF.v.NC$row2)
res.diagdds.pruned.genus.fecal.physeq.LF.v.NC$row2 <- gsub('\\s+', '|', res.diagdds.pruned.genus.fecal.physeq.LF.v.NC$row2)
res.diagdds.pruned.genus.fecal.physeq.HF.v.LF$row2 <- gsub('\\s+', '|', res.diagdds.pruned.genus.fecal.physeq.HF.v.LF$row2)

# res <- as.data.frame(res)
res.diagdds.pruned.genus.cecum.physeq.HF.v.NC <- as.data.frame(res.diagdds.pruned.genus.cecum.physeq.HF.v.NC)
res.diagdds.pruned.genus.cecum.physeq.LF.v.NC <- as.data.frame(res.diagdds.pruned.genus.cecum.physeq.LF.v.NC)
res.diagdds.pruned.genus.cecum.physeq.HF.v.LF <- as.data.frame(res.diagdds.pruned.genus.cecum.physeq.HF.v.LF)
res.diagdds.pruned.genus.fecal.physeq.HF.v.NC <- as.data.frame(res.diagdds.pruned.genus.fecal.physeq.HF.v.NC)
res.diagdds.pruned.genus.fecal.physeq.LF.v.NC <- as.data.frame(res.diagdds.pruned.genus.fecal.physeq.LF.v.NC)
res.diagdds.pruned.genus.fecal.physeq.HF.v.LF <- as.data.frame(res.diagdds.pruned.genus.fecal.physeq.HF.v.LF)

rownames(res.diagdds.subset.otu.table.10.8.v.PBS.con)

# res$names <- res$Taxa
# res$Taxa <- res$row2
# otu.to.save.res <-as.character(res$names)

res.diagdds.pruned.genus.cecum.physeq.HF.v.NC$names <- res.diagdds.pruned.genus.cecum.physeq.HF.v.NC$Row.names
res.diagdds.pruned.genus.cecum.physeq.HF.v.NC$Row.names <- res.diagdds.pruned.genus.cecum.physeq.HF.v.NC$row2
otu.to.save.res.diagdds.pruned.genus.cecum.physeq.HF.v.NC <-as.character(res.diagdds.pruned.genus.cecum.physeq.HF.v.NC$names)

res.diagdds.pruned.genus.cecum.physeq.LF.v.NC$names <- res.diagdds.pruned.genus.cecum.physeq.LF.v.NC$Row.names
res.diagdds.pruned.genus.cecum.physeq.LF.v.NC$Row.names <- res.diagdds.pruned.genus.cecum.physeq.LF.v.NC$row2
otu.to.save.res.diagdds.pruned.genus.cecum.physeq.LF.v.NC <-as.character(res.diagdds.pruned.genus.cecum.physeq.LF.v.NC$names)

res.diagdds.pruned.genus.cecum.physeq.HF.v.LF$names <- res.diagdds.pruned.genus.cecum.physeq.HF.v.LF$Row.names
res.diagdds.pruned.genus.cecum.physeq.HF.v.LF$Row.names <- res.diagdds.pruned.genus.cecum.physeq.HF.v.LF$row2
otu.to.save.res.diagdds.pruned.genus.cecum.physeq.HF.v.LF <-as.character(res.diagdds.pruned.genus.cecum.physeq.HF.v.LF$names)

res.diagdds.pruned.genus.fecal.physeq.HF.v.NC$names <- res.diagdds.pruned.genus.fecal.physeq.HF.v.NC$Row.names
res.diagdds.pruned.genus.fecal.physeq.HF.v.NC$Row.names <- res.diagdds.pruned.genus.fecal.physeq.HF.v.NC$row2
otu.to.save.res.diagdds.pruned.genus.fecal.physeq.HF.v.NC <-as.character(res.diagdds.pruned.genus.fecal.physeq.HF.v.NC$names)

res.diagdds.pruned.genus.fecal.physeq.LF.v.NC$names <- res.diagdds.pruned.genus.fecal.physeq.LF.v.NC$Row.names
res.diagdds.pruned.genus.fecal.physeq.LF.v.NC$Row.names <- res.diagdds.pruned.genus.fecal.physeq.LF.v.NC$row2
otu.to.save.res.diagdds.pruned.genus.fecal.physeq.LF.v.NC <-as.character(res.diagdds.pruned.genus.fecal.physeq.LF.v.NC$names)

res.diagdds.pruned.genus.fecal.physeq.HF.v.LF$names <- res.diagdds.pruned.genus.fecal.physeq.HF.v.LF$Row.names
res.diagdds.pruned.genus.fecal.physeq.HF.v.LF$Row.names <- res.diagdds.pruned.genus.fecal.physeq.HF.v.LF$row2
otu.to.save.res.diagdds.pruned.genus.fecal.physeq.HF.v.LF <-as.character(res.diagdds.pruned.genus.fecal.physeq.HF.v.LF$names)

pruned.genus.rel.table.cecum.physeq.HF.v.NC
pruned.genus.rel.table.cecum.physeq.LF.v.NC
pruned.genus.rel.table.cecum.physeq.HF.v.LF
pruned.genus.rel.table.fecal.physeq.HF.v.NC
pruned.genus.rel.table.fecal.physeq.LF.v.NC
pruned.genus.rel.table.fecal.physeq.HF.v.LF

res.diagdds.pruned.genus.cecum.physeq.HF.v.NC
res.diagdds.pruned.genus.cecum.physeq.LF.v.NC
res.diagdds.pruned.genus.cecum.physeq.HF.v.LF
res.diagdds.pruned.genus.fecal.physeq.HF.v.NC
res.diagdds.pruned.genus.fecal.physeq.LF.v.NC
res.diagdds.pruned.genus.fecal.physeq.HF.v.LF

otu.to.save.res.diagdds.pruned.genus.cecum.physeq.HF.v.NC
otu.to.save.res.diagdds.pruned.genus.cecum.physeq.LF.v.NC
otu.to.save.res.diagdds.pruned.genus.cecum.physeq.HF.v.LF
otu.to.save.res.diagdds.pruned.genus.fecal.physeq.HF.v.NC
otu.to.save.res.diagdds.pruned.genus.fecal.physeq.LF.v.NC
otu.to.save.res.diagdds.pruned.genus.fecal.physeq.HF.v.LF

# change variable of interest
experiment.pruned.genus.rel.table = subset_samples(pruned.genus.rel.table.cecum.physeq.HF.v.NC, Treatment %in% c("HF"))
control.pruned.genus.rel.table = subset_samples(pruned.genus.rel.table.cecum.physeq.HF.v.NC, Treatment %in% c("NC"))
experiment.pruned.genus.rel.table.df <- data.frame(otu_table(experiment.pruned.genus.rel.table))
experiment.pruned.genus.rel.table.df.meanRA <- rowMeans(experiment.pruned.genus.rel.table.df)
control.pruned.genus.rel.table.df <- data.frame(otu_table(control.pruned.genus.rel.table))
control.pruned.genus.rel.table.df.meanRA <- rowMeans(control.pruned.genus.rel.table.df)
experiment.pruned.genus.rel.table.df.meanRA.save <- experiment.pruned.genus.rel.table.df.meanRA[otu.to.save.res.diagdds.pruned.genus.cecum.physeq.HF.v.NC]
control.pruned.genus.rel.table.df.meanRA.save <- control.pruned.genus.rel.table.df.meanRA[otu.to.save.res.diagdds.pruned.genus.cecum.physeq.HF.v.NC]
res.diagdds.pruned.genus.cecum.physeq.HF.v.NC$abundance.experiment <- experiment.pruned.genus.rel.table.df.meanRA.save
res.diagdds.pruned.genus.cecum.physeq.HF.v.NC$abundance.control <- control.pruned.genus.rel.table.df.meanRA.save

experiment.pruned.genus.rel.table = subset_samples(pruned.genus.rel.table.cecum.physeq.LF.v.NC, Treatment %in% c("LF"))
control.pruned.genus.rel.table = subset_samples(pruned.genus.rel.table.cecum.physeq.LF.v.NC, Treatment %in% c("NC"))
experiment.pruned.genus.rel.table.df <- data.frame(otu_table(experiment.pruned.genus.rel.table))
experiment.pruned.genus.rel.table.df.meanRA <- rowMeans(experiment.pruned.genus.rel.table.df)
control.pruned.genus.rel.table.df <- data.frame(otu_table(control.pruned.genus.rel.table))
control.pruned.genus.rel.table.df.meanRA <- rowMeans(control.pruned.genus.rel.table.df)
experiment.pruned.genus.rel.table.df.meanRA.save <- experiment.pruned.genus.rel.table.df.meanRA[otu.to.save.res.diagdds.pruned.genus.cecum.physeq.LF.v.NC]
control.pruned.genus.rel.table.df.meanRA.save <- control.pruned.genus.rel.table.df.meanRA[otu.to.save.res.diagdds.pruned.genus.cecum.physeq.LF.v.NC]
res.diagdds.pruned.genus.cecum.physeq.LF.v.NC$abundance.experiment <- experiment.pruned.genus.rel.table.df.meanRA.save
res.diagdds.pruned.genus.cecum.physeq.LF.v.NC$abundance.control <- control.pruned.genus.rel.table.df.meanRA.save

experiment.pruned.genus.rel.table = subset_samples(pruned.genus.rel.table.cecum.physeq.HF.v.LF, Treatment %in% c("HF"))
control.pruned.genus.rel.table = subset_samples(pruned.genus.rel.table.cecum.physeq.HF.v.LF, Treatment %in% c("LF"))
experiment.pruned.genus.rel.table.df <- data.frame(otu_table(experiment.pruned.genus.rel.table))
experiment.pruned.genus.rel.table.df.meanRA <- rowMeans(experiment.pruned.genus.rel.table.df)
control.pruned.genus.rel.table.df <- data.frame(otu_table(control.pruned.genus.rel.table))
control.pruned.genus.rel.table.df.meanRA <- rowMeans(control.pruned.genus.rel.table.df)
experiment.pruned.genus.rel.table.df.meanRA.save <- experiment.pruned.genus.rel.table.df.meanRA[otu.to.save.res.diagdds.pruned.genus.cecum.physeq.HF.v.LF]
control.pruned.genus.rel.table.df.meanRA.save <- control.pruned.genus.rel.table.df.meanRA[otu.to.save.res.diagdds.pruned.genus.cecum.physeq.HF.v.LF]
res.diagdds.pruned.genus.cecum.physeq.HF.v.LF$abundance.experiment <- experiment.pruned.genus.rel.table.df.meanRA.save
res.diagdds.pruned.genus.cecum.physeq.HF.v.LF$abundance.control <- control.pruned.genus.rel.table.df.meanRA.save

experiment.pruned.genus.rel.table = subset_samples(pruned.genus.rel.table.fecal.physeq.HF.v.NC, Treatment %in% c("HF"))
control.pruned.genus.rel.table = subset_samples(pruned.genus.rel.table.fecal.physeq.HF.v.NC, Treatment %in% c("NC"))
experiment.pruned.genus.rel.table.df <- data.frame(otu_table(experiment.pruned.genus.rel.table))
experiment.pruned.genus.rel.table.df.meanRA <- rowMeans(experiment.pruned.genus.rel.table.df)
control.pruned.genus.rel.table.df <- data.frame(otu_table(control.pruned.genus.rel.table))
control.pruned.genus.rel.table.df.meanRA <- rowMeans(control.pruned.genus.rel.table.df)
experiment.pruned.genus.rel.table.df.meanRA.save <- experiment.pruned.genus.rel.table.df.meanRA[otu.to.save.res.diagdds.pruned.genus.fecal.physeq.HF.v.NC]
control.pruned.genus.rel.table.df.meanRA.save <- control.pruned.genus.rel.table.df.meanRA[otu.to.save.res.diagdds.pruned.genus.fecal.physeq.HF.v.NC]
res.diagdds.pruned.genus.fecal.physeq.HF.v.NC$abundance.experiment <- experiment.pruned.genus.rel.table.df.meanRA.save
res.diagdds.pruned.genus.fecal.physeq.HF.v.NC$abundance.control <- control.pruned.genus.rel.table.df.meanRA.save

experiment.pruned.genus.rel.table = subset_samples(pruned.genus.rel.table.fecal.physeq.LF.v.NC, Treatment %in% c("LF"))
control.pruned.genus.rel.table = subset_samples(pruned.genus.rel.table.fecal.physeq.LF.v.NC, Treatment %in% c("NC"))
experiment.pruned.genus.rel.table.df <- data.frame(otu_table(experiment.pruned.genus.rel.table))
experiment.pruned.genus.rel.table.df.meanRA <- rowMeans(experiment.pruned.genus.rel.table.df)
control.pruned.genus.rel.table.df <- data.frame(otu_table(control.pruned.genus.rel.table))
control.pruned.genus.rel.table.df.meanRA <- rowMeans(control.pruned.genus.rel.table.df)
experiment.pruned.genus.rel.table.df.meanRA.save <- experiment.pruned.genus.rel.table.df.meanRA[otu.to.save.res.diagdds.pruned.genus.fecal.physeq.LF.v.NC]
control.pruned.genus.rel.table.df.meanRA.save <- control.pruned.genus.rel.table.df.meanRA[otu.to.save.res.diagdds.pruned.genus.fecal.physeq.LF.v.NC]
res.diagdds.pruned.genus.fecal.physeq.LF.v.NC$abundance.experiment <- experiment.pruned.genus.rel.table.df.meanRA.save
res.diagdds.pruned.genus.fecal.physeq.LF.v.NC$abundance.control <- control.pruned.genus.rel.table.df.meanRA.save

experiment.pruned.genus.rel.table = subset_samples(pruned.genus.rel.table.fecal.physeq.HF.v.LF, Treatment %in% c("HF"))
control.pruned.genus.rel.table = subset_samples(pruned.genus.rel.table.fecal.physeq.HF.v.LF, Treatment %in% c("LF"))
experiment.pruned.genus.rel.table.df <- data.frame(otu_table(experiment.pruned.genus.rel.table))
experiment.pruned.genus.rel.table.df.meanRA <- rowMeans(experiment.pruned.genus.rel.table.df)
control.pruned.genus.rel.table.df <- data.frame(otu_table(control.pruned.genus.rel.table))
control.pruned.genus.rel.table.df.meanRA <- rowMeans(control.pruned.genus.rel.table.df)
experiment.pruned.genus.rel.table.df.meanRA.save <- experiment.pruned.genus.rel.table.df.meanRA[otu.to.save.res.diagdds.pruned.genus.fecal.physeq.HF.v.LF]
control.pruned.genus.rel.table.df.meanRA.save <- control.pruned.genus.rel.table.df.meanRA[otu.to.save.res.diagdds.pruned.genus.fecal.physeq.HF.v.LF]
res.diagdds.pruned.genus.fecal.physeq.HF.v.LF$abundance.experiment <- experiment.pruned.genus.rel.table.df.meanRA.save
res.diagdds.pruned.genus.fecal.physeq.HF.v.LF$abundance.control <- control.pruned.genus.rel.table.df.meanRA.save

res.diagdds.pruned.genus.cecum.physeq.HF.v.NC
res.diagdds.pruned.genus.cecum.physeq.LF.v.NC
res.diagdds.pruned.genus.cecum.physeq.HF.v.LF
res.diagdds.pruned.genus.fecal.physeq.HF.v.NC
res.diagdds.pruned.genus.fecal.physeq.LF.v.NC
res.diagdds.pruned.genus.fecal.physeq.HF.v.LF

names(res.diagdds.pruned.genus.cecum.physeq.HF.v.NC)
names(res.diagdds.pruned.genus.cecum.physeq.LF.v.NC)
names(res.diagdds.pruned.genus.cecum.physeq.HF.v.LF)
names(res.diagdds.pruned.genus.fecal.physeq.HF.v.NC)
names(res.diagdds.pruned.genus.fecal.physeq.LF.v.NC)
names(res.diagdds.pruned.genus.fecal.physeq.HF.v.LF)

# Change abundance column headers
# res.1 <- res[,c("Taxa", "abundance.experiment", "abundance.control", "baseMean", "logFC", "pvalue", "adj.P.Val")]
# 
# change table name
# write.table(res.1,file="20220311.All_skin_vs_sheath_taxa_abundance.txt", sep="\t", col.names = NA, row.names = TRUE, quote=FALSE)

res.diagdds.pruned.genus.cecum.physeq.HF.v.NC$sig <- -log10(res.diagdds.pruned.genus.cecum.physeq.HF.v.NC$padj)
res.diagdds.pruned.genus.cecum.physeq.LF.v.NC$sig <- -log10(res.diagdds.pruned.genus.cecum.physeq.LF.v.NC$padj)
res.diagdds.pruned.genus.cecum.physeq.HF.v.LF$sig <- -log10(res.diagdds.pruned.genus.cecum.physeq.HF.v.LF$padj)
res.diagdds.pruned.genus.fecal.physeq.HF.v.NC$sig <- -log10(res.diagdds.pruned.genus.fecal.physeq.HF.v.NC$padj)
res.diagdds.pruned.genus.fecal.physeq.LF.v.NC$sig <- -log10(res.diagdds.pruned.genus.fecal.physeq.LF.v.NC$padj)
res.diagdds.pruned.genus.fecal.physeq.HF.v.LF$sig <- -log10(res.diagdds.pruned.genus.fecal.physeq.HF.v.LF$padj)

sum(is.infinite(res.diagdds.pruned.genus.cecum.physeq.HF.v.NC$sig))
sum(is.infinite(res.diagdds.pruned.genus.cecum.physeq.LF.v.NC$sig))
sum(is.infinite(res.diagdds.pruned.genus.cecum.physeq.HF.v.LF$sig))
sum(is.infinite(res.diagdds.pruned.genus.fecal.physeq.HF.v.NC$sig))
sum(is.infinite(res.diagdds.pruned.genus.fecal.physeq.LF.v.NC$sig))
sum(is.infinite(res.diagdds.pruned.genus.fecal.physeq.HF.v.LF$sig))

cols.res.diagdds.pruned.genus.cecum.physeq.HF.v.NC <- densCols(res.diagdds.pruned.genus.cecum.physeq.HF.v.NC$log2FoldChange, res.diagdds.pruned.genus.cecum.physeq.HF.v.NC$sig)
cols.res.diagdds.pruned.genus.cecum.physeq.LF.v.NC <- densCols(res.diagdds.pruned.genus.cecum.physeq.LF.v.NC$log2FoldChange, res.diagdds.pruned.genus.cecum.physeq.LF.v.NC$sig)
cols.res.diagdds.pruned.genus.cecum.physeq.HF.v.LF <- densCols(res.diagdds.pruned.genus.cecum.physeq.HF.v.LF$log2FoldChange, res.diagdds.pruned.genus.cecum.physeq.HF.v.LF$sig)
cols.res.diagdds.pruned.genus.fecal.physeq.HF.v.NC <- densCols(res.diagdds.pruned.genus.fecal.physeq.HF.v.NC$log2FoldChange, res.diagdds.pruned.genus.fecal.physeq.HF.v.NC$sig)
cols.res.diagdds.pruned.genus.fecal.physeq.LF.v.NC <- densCols(res.diagdds.pruned.genus.fecal.physeq.LF.v.NC$log2FoldChange, res.diagdds.pruned.genus.fecal.physeq.LF.v.NC$sig)
cols.res.diagdds.pruned.genus.fecal.physeq.HF.v.LF <- densCols(res.diagdds.pruned.genus.fecal.physeq.HF.v.LF$log2FoldChange, res.diagdds.pruned.genus.fecal.physeq.HF.v.LF$sig)

cols.res.diagdds.pruned.genus.cecum.physeq.HF.v.NC[res.diagdds.pruned.genus.cecum.physeq.HF.v.NC$pvalue==0] <- "purple"
cols.res.diagdds.pruned.genus.cecum.physeq.LF.v.NC[res.diagdds.pruned.genus.cecum.physeq.LF.v.NC$pvalue==0] <- "purple"
cols.res.diagdds.pruned.genus.cecum.physeq.HF.v.LF[res.diagdds.pruned.genus.cecum.physeq.HF.v.LF$pvalue==0] <- "purple"
cols.res.diagdds.pruned.genus.fecal.physeq.HF.v.NC[res.diagdds.pruned.genus.fecal.physeq.HF.v.NC$pvalue==0] <- "purple"
cols.res.diagdds.pruned.genus.fecal.physeq.LF.v.NC[res.diagdds.pruned.genus.fecal.physeq.LF.v.NC$pvalue==0] <- "purple"
cols.res.diagdds.pruned.genus.fecal.physeq.HF.v.LF[res.diagdds.pruned.genus.fecal.physeq.HF.v.LF$pvalue==0] <- "purple"

write.table(res.diagdds.pruned.genus.cecum.physeq.HF.v.NC, file="res.diagdds.pruned.genus.cecum.physeq.HF.v.NC.txt", sep="\t")
write.table(res.diagdds.pruned.genus.cecum.physeq.LF.v.NC, file="res.diagdds.pruned.genus.cecum.physeq.LF.v.NC.txt", sep="\t")
write.table(res.diagdds.pruned.genus.cecum.physeq.HF.v.LF, file="res.diagdds.pruned.genus.cecum.physeq.HF.v.LF.txt", sep="\t")
write.table(res.diagdds.pruned.genus.fecal.physeq.HF.v.NC, file="res.diagdds.pruned.genus.fecal.physeq.HF.v.NC.txt", sep="\t")
write.table(res.diagdds.pruned.genus.fecal.physeq.LF.v.NC, file="res.diagdds.pruned.genus.fecal.physeq.LF.v.NC.txt", sep="\t")
write.table(res.diagdds.pruned.genus.fecal.physeq.HF.v.LF, file="res.diagdds.pruned.genus.fecal.physeq.HF.v.LF.txt", sep="\t")

# Change colours printing out
cols.res.diagdds.pruned.genus.cecum.physeq.HF.v.NC[res.diagdds.pruned.genus.cecum.physeq.HF.v.NC$log2FoldChange > 0 & res.diagdds.pruned.genus.cecum.physeq.HF.v.NC$padj < alpha ] <- "brown"
cols.res.diagdds.pruned.genus.cecum.physeq.HF.v.NC[res.diagdds.pruned.genus.cecum.physeq.HF.v.NC$log2FoldChange < 0 & res.diagdds.pruned.genus.cecum.physeq.HF.v.NC$padj < alpha ] <- "mediumseagreen"

res.diagdds.pruned.genus.cecum.physeq.HF.v.NC$pch <- 19
res.diagdds.pruned.genus.cecum.physeq.HF.v.NC$pch[res.diagdds.pruned.genus.cecum.physeq.HF.v.NC$pvalue ==0] <- 6

pdf(file="DESEQ2_HF_v_NC_FDR_0.2_cec.pdf", width=8, height=8)
ggplot(res.diagdds.pruned.genus.cecum.physeq.HF.v.NC, aes(x = log2FoldChange, y = sig,label=Row.names)) +
  geom_point(color=cols.res.diagdds.pruned.genus.cecum.physeq.HF.v.NC, size=ifelse(res.diagdds.pruned.genus.cecum.physeq.HF.v.NC$log2FoldChange>=1 & res.diagdds.pruned.genus.cecum.physeq.HF.v.NC$padj < alpha, 100 * res.diagdds.pruned.genus.cecum.physeq.HF.v.NC$abundance.experiment, ifelse(res.diagdds.pruned.genus.cecum.physeq.HF.v.NC$log2FoldChange<=-1 & res.diagdds.pruned.genus.cecum.physeq.HF.v.NC$padj < alpha, 100 * res.diagdds.pruned.genus.cecum.physeq.HF.v.NC$abundance.control, 2)),alpha=0.7) +
  geom_text_repel(aes(label=ifelse(res.diagdds.pruned.genus.cecum.physeq.HF.v.NC$padj < alpha & res.diagdds.pruned.genus.cecum.physeq.HF.v.NC$Genus!="g__", as.character(res.diagdds.pruned.genus.cecum.physeq.HF.v.NC$Genus),'')),size=4,force=25, segment.colour="grey",segment.alpha=0.2) +
  geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") +
  xlab("Effect size: log2(fold-change)") +
  ylab("-log10(adjusted p-value)") +
  theme
dev.off()

# Change colours printing out
cols.res.diagdds.pruned.genus.cecum.physeq.LF.v.NC[res.diagdds.pruned.genus.cecum.physeq.LF.v.NC$log2FoldChange > 0 & res.diagdds.pruned.genus.cecum.physeq.LF.v.NC$padj < alpha ] <- "brown"
cols.res.diagdds.pruned.genus.cecum.physeq.LF.v.NC[res.diagdds.pruned.genus.cecum.physeq.LF.v.NC$log2FoldChange < 0 & res.diagdds.pruned.genus.cecum.physeq.LF.v.NC$padj < alpha ] <- "mediumseagreen"

res.diagdds.pruned.genus.cecum.physeq.LF.v.NC$pch <- 19
res.diagdds.pruned.genus.cecum.physeq.LF.v.NC$pch[res.diagdds.pruned.genus.cecum.physeq.LF.v.NC$pvalue ==0] <- 6

pdf(file="DESE2_LF_v_NC_FDR_0.2_cec.pdf", width=8, height=8)
ggplot(res.diagdds.pruned.genus.cecum.physeq.LF.v.NC, aes(x = log2FoldChange, y = sig,label=Row.names)) +
  geom_point(color=cols.res.diagdds.pruned.genus.cecum.physeq.LF.v.NC, size=ifelse(res.diagdds.pruned.genus.cecum.physeq.LF.v.NC$log2FoldChange>=1 & res.diagdds.pruned.genus.cecum.physeq.LF.v.NC$padj < alpha, 100 * res.diagdds.pruned.genus.cecum.physeq.LF.v.NC$abundance.experiment, ifelse(res.diagdds.pruned.genus.cecum.physeq.LF.v.NC$log2FoldChange<=-1 & res.diagdds.pruned.genus.cecum.physeq.LF.v.NC$padj < alpha, 100 * res.diagdds.pruned.genus.cecum.physeq.LF.v.NC$abundance.control, 2)),alpha=0.7) +
  geom_text_repel(aes(label=ifelse(res.diagdds.pruned.genus.cecum.physeq.LF.v.NC$padj < alpha & res.diagdds.pruned.genus.cecum.physeq.LF.v.NC$Genus!="g__", as.character(res.diagdds.pruned.genus.cecum.physeq.LF.v.NC$Genus),'')),size=4,force=25, segment.colour="grey",segment.alpha=0.2) +
  geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") +
  xlab("Effect size: log2(fold-change)") +
  ylab("-log10(adjusted p-value)") +
  theme
dev.off()

# Change colours printing out
cols.res.diagdds.pruned.genus.cecum.physeq.HF.v.LF[res.diagdds.pruned.genus.cecum.physeq.HF.v.LF$log2FoldChange > 0 & res.diagdds.pruned.genus.cecum.physeq.HF.v.LF$padj < alpha ] <- "brown"
cols.res.diagdds.pruned.genus.cecum.physeq.HF.v.LF[res.diagdds.pruned.genus.cecum.physeq.HF.v.LF$log2FoldChange < 0 & res.diagdds.pruned.genus.cecum.physeq.HF.v.LF$padj < alpha ] <- "mediumseagreen"

res.diagdds.pruned.genus.cecum.physeq.HF.v.LF$pch <- 19
res.diagdds.pruned.genus.cecum.physeq.HF.v.LF$pch[res.diagdds.pruned.genus.cecum.physeq.HF.v.LF$pvalue ==0] <- 6

pdf(file="DESEQ2_HF_v_LF_FDR_0.2_cec.pdf", width=8, height=8)
ggplot(res.diagdds.pruned.genus.cecum.physeq.HF.v.LF, aes(x = log2FoldChange, y = sig,label=Row.names)) +
  geom_point(color=cols.res.diagdds.pruned.genus.cecum.physeq.HF.v.LF, size=ifelse(res.diagdds.pruned.genus.cecum.physeq.HF.v.LF$log2FoldChange>=1 & res.diagdds.pruned.genus.cecum.physeq.HF.v.LF$padj < alpha, 100 * res.diagdds.pruned.genus.cecum.physeq.HF.v.LF$abundance.experiment, ifelse(res.diagdds.pruned.genus.cecum.physeq.HF.v.LF$log2FoldChange<=-1 & res.diagdds.pruned.genus.cecum.physeq.HF.v.LF$padj < alpha, 100 * res.diagdds.pruned.genus.cecum.physeq.HF.v.LF$abundance.control, 2)),alpha=0.7) +
  geom_text_repel(aes(label=ifelse(res.diagdds.pruned.genus.cecum.physeq.HF.v.LF$padj < alpha & res.diagdds.pruned.genus.cecum.physeq.HF.v.LF$Genus!="g__", as.character(res.diagdds.pruned.genus.cecum.physeq.HF.v.LF$Genus),'')),size=4,force=25, segment.colour="grey",segment.alpha=0.2) +
  geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") +
  xlab("Effect size: log2(fold-change)") +
  ylab("-log10(adjusted p-value)") +
  theme
dev.off()

# Change colours printing out
cols.res.diagdds.pruned.genus.fecal.physeq.HF.v.NC[res.diagdds.pruned.genus.fecal.physeq.HF.v.NC$log2FoldChange > 0 & res.diagdds.pruned.genus.fecal.physeq.HF.v.NC$padj < alpha ] <- "brown"
cols.res.diagdds.pruned.genus.fecal.physeq.HF.v.NC[res.diagdds.pruned.genus.fecal.physeq.HF.v.NC$log2FoldChange < 0 & res.diagdds.pruned.genus.fecal.physeq.HF.v.NC$padj < alpha ] <- "mediumseagreen"

res.diagdds.pruned.genus.fecal.physeq.HF.v.NC$pch <- 19
res.diagdds.pruned.genus.fecal.physeq.HF.v.NC$pch[res.diagdds.pruned.genus.fecal.physeq.HF.v.NC$pvalue ==0] <- 6

pdf(file="DESEQ2_HF_v_NC_FDR_0.2_fec.pdf", width=8, height=8)
ggplot(res.diagdds.pruned.genus.fecal.physeq.HF.v.NC, aes(x = log2FoldChange, y = sig,label=Row.names)) +
  geom_point(color=cols.res.diagdds.pruned.genus.fecal.physeq.HF.v.NC, size=ifelse(res.diagdds.pruned.genus.fecal.physeq.HF.v.NC$log2FoldChange>=1 & res.diagdds.pruned.genus.fecal.physeq.HF.v.NC$padj < alpha, 100 * res.diagdds.pruned.genus.fecal.physeq.HF.v.NC$abundance.experiment, ifelse(res.diagdds.pruned.genus.fecal.physeq.HF.v.NC$log2FoldChange<=-1 & res.diagdds.pruned.genus.fecal.physeq.HF.v.NC$padj < alpha, 100 * res.diagdds.pruned.genus.fecal.physeq.HF.v.NC$abundance.control, 2)),alpha=0.7) +
  geom_text_repel(aes(label=ifelse(res.diagdds.pruned.genus.fecal.physeq.HF.v.NC$padj < alpha & res.diagdds.pruned.genus.fecal.physeq.HF.v.NC$Genus!="g__", as.character(res.diagdds.pruned.genus.fecal.physeq.HF.v.NC$Genus),'')),size=4,force=25, segment.colour="grey",segment.alpha=0.2) +
  geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") +
  xlab("Effect size: log2(fold-change)") +
  ylab("-log10(adjusted p-value)") +
  theme
dev.off()

# Change colours printing out
cols.res.diagdds.pruned.genus.fecal.physeq.LF.v.NC[res.diagdds.pruned.genus.fecal.physeq.LF.v.NC$log2FoldChange > 0 & res.diagdds.pruned.genus.fecal.physeq.LF.v.NC$padj < alpha ] <- "brown"
cols.res.diagdds.pruned.genus.fecal.physeq.LF.v.NC[res.diagdds.pruned.genus.fecal.physeq.LF.v.NC$log2FoldChange < 0 & res.diagdds.pruned.genus.fecal.physeq.LF.v.NC$padj < alpha ] <- "mediumseagreen"

res.diagdds.pruned.genus.fecal.physeq.LF.v.NC$pch <- 19
res.diagdds.pruned.genus.fecal.physeq.LF.v.NC$pch[res.diagdds.pruned.genus.fecal.physeq.LF.v.NC$pvalue ==0] <- 6

pdf(file="DESEQ2_LF_v_NC_FDR_0.2_fec.pdf", width=8, height=8)
ggplot(res.diagdds.pruned.genus.fecal.physeq.LF.v.NC, aes(x = log2FoldChange, y = sig,label=Row.names)) +
  geom_point(color=cols.res.diagdds.pruned.genus.fecal.physeq.LF.v.NC, size=ifelse(res.diagdds.pruned.genus.fecal.physeq.LF.v.NC$log2FoldChange>=1 & res.diagdds.pruned.genus.fecal.physeq.LF.v.NC$padj < alpha, 100 * res.diagdds.pruned.genus.fecal.physeq.LF.v.NC$abundance.experiment, ifelse(res.diagdds.pruned.genus.fecal.physeq.LF.v.NC$log2FoldChange<=-1 & res.diagdds.pruned.genus.fecal.physeq.LF.v.NC$padj < alpha, 100 * res.diagdds.pruned.genus.fecal.physeq.LF.v.NC$abundance.control, 2)),alpha=0.7) +
  geom_text_repel(aes(label=ifelse(res.diagdds.pruned.genus.fecal.physeq.LF.v.NC$padj < alpha & res.diagdds.pruned.genus.fecal.physeq.LF.v.NC$Genus!="g__", as.character(res.diagdds.pruned.genus.fecal.physeq.LF.v.NC$Genus),'')),size=4,force=25, segment.colour="grey",segment.alpha=0.2) +
  geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") +
  xlab("Effect size: log2(fold-change)") +
  ylab("-log10(adjusted p-value)") +
  theme
dev.off()

# Change colours printing out
cols.res.diagdds.pruned.genus.fecal.physeq.HF.v.LF[res.diagdds.pruned.genus.fecal.physeq.HF.v.LF$log2FoldChange > 0 & res.diagdds.pruned.genus.fecal.physeq.HF.v.LF$padj < alpha ] <- "brown"
cols.res.diagdds.pruned.genus.fecal.physeq.HF.v.LF[res.diagdds.pruned.genus.fecal.physeq.HF.v.LF$log2FoldChange < 0 & res.diagdds.pruned.genus.fecal.physeq.HF.v.LF$padj < alpha ] <- "mediumseagreen"

res.diagdds.pruned.genus.fecal.physeq.HF.v.LF$pch <- 19
res.diagdds.pruned.genus.fecal.physeq.HF.v.LF$pch[res.diagdds.pruned.genus.fecal.physeq.HF.v.LF$pvalue ==0] <- 6

pdf(file="DESEQ2_HF_v_LF_FDR_0.2_fec.pdf", width=8, height=8)
ggplot(res.diagdds.pruned.genus.fecal.physeq.HF.v.LF, aes(x = log2FoldChange, y = sig,label=Row.names)) +
  geom_point(color=cols.res.diagdds.pruned.genus.fecal.physeq.HF.v.LF, size=ifelse(res.diagdds.pruned.genus.fecal.physeq.HF.v.LF$log2FoldChange>=1 & res.diagdds.pruned.genus.fecal.physeq.HF.v.LF$padj < alpha, 100 * res.diagdds.pruned.genus.fecal.physeq.HF.v.LF$abundance.experiment, ifelse(res.diagdds.pruned.genus.fecal.physeq.HF.v.LF$log2FoldChange<=-1 & res.diagdds.pruned.genus.fecal.physeq.HF.v.LF$padj < alpha, 100 * res.diagdds.pruned.genus.fecal.physeq.HF.v.LF$abundance.control, 2)),alpha=0.7) +
  geom_text_repel(aes(label=ifelse(res.diagdds.pruned.genus.fecal.physeq.HF.v.LF$padj < alpha & res.diagdds.pruned.genus.fecal.physeq.HF.v.LF$Genus!="g__", as.character(res.diagdds.pruned.genus.fecal.physeq.HF.v.LF$Genus),'')),size=4,force=25, segment.colour="grey",segment.alpha=0.2) +
  geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") +
  xlab("Effect size: log2(fold-change)") +
  ylab("-log10(adjusted p-value)") +
  theme
dev.off()
