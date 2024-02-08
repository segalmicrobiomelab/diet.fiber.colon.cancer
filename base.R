# Session.info()
sessionInfo()
getwd()
setwd("/set/your/working/directory")

#Load Phyloseq
library(phyloseq)
library(ade4)
library(vegan)
library(biomformat)
library(devtools)
library(readr)
library(readtext)
library(qiime2R)
library("ape")
library()
library(ggplot2)
library(ggrepel)
library(DESeq2)
library(microbiomeMarker)

# install packages if needed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
Biocmanager::install()
BiocManager::install("phyloseq")
install.packages('devtools')
install.packages('tidyverse')
install.packages('readr')
install.packages('readtext')
install.packages('vegan')
install.packages('ade4')
install.packages('biomformat')
install.packages('qiime2R')
install.packages('phyloseq')
BiocManager::install("DESeq2")

# # install packages if needed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
Biocmanager::install()
BiocManager::install("pathfindR")
BiocManager::install("MatrixModels")
BiocManager::install("Matrix", force = TRUE)
BiocManager::install("sf", force = TRUE)
BiocManager::install("phangorn")
if (!requireNamespace("devtools", quietly = TRUE)){install.packages("devtools")}
devtools::install_github("jbisanz/qiime2R")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("microbiomeMarker")

# Load the QIIME2 packages 
# QIIME2 objects from run 
# metadata from ... your metadata 
physeq<-qza_to_phyloseq(
  features="table.qza",
  tree="rooted-tree_quality.qza",
  taxonomy="taxonomy.qza",
  metadata = "msq.fiber.map.txt"
)

rownames(sample_data(physeq))
colnames(sample_data(physeq))

# Remove taxa with 0 abundance
physeq = subset_taxa(physeq, rowSums(otu_table(physeq)) != 0)

##If you want to nomalize OTU table before
## To normalize data you need to set a function
normalizeSample = function(x) {
  x/sum(x)
}
otu.relative.table = transformSampleCounts(physeq, normalizeSample)

rownames(sample_data(otu.relative.table))
colnames(sample_data(otu.relative.table))

Phylum.rel.table = tax_glom(otu.relative.table, taxrank = "Phylum")
Class.rel.table = tax_glom(otu.relative.table, taxrank = "Class")
Order.rel.table = tax_glom(otu.relative.table, taxrank = "Order")
Family.rel.table = tax_glom(otu.relative.table, taxrank = "Family")
Genus.rel.table = tax_glom(otu.relative.table, taxrank = "Genus")
OTU.rel.table = tax_glom(otu.relative.table, taxrank = "Species")

rownames(sample_data(otu.relative.table))
colnames(sample_data(otu.relative.table))
colnames(sample_data(otu.relative.table))