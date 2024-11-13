# LEFSe 
# Old packages 

#Load Packages
library(DESeq2)
library(edgeR)
library(limma)
library(Glimma)
library(gplots)
library(RColorBrewer)
library(pheatmap)
library(ggplot2)
library(ggrepel)
library(pathfindR) #
library(scales)
library("quantreg")
library(data.table)
library(fBasics)
library(forcats)
library(omu)
# library(maptools)
library(phyloseq)
library(vegan)
library(ggpmisc) #
library(dplyr)
library(tibble)
library(formattable) #
library("htmltools")
library("webshot")  #  
library(splitstackshape) #
library(decontam)
library(dplyr)
library(grid)
library(cowplot)
library(wesanderson) # 
library(colorspace)
library(lefser) # 
library(microbiomeMarker) # 
library(knitr)
library(dplyr)


# # install packages if needed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
Biocmanager::install()
BiocManager::install("pathfindR")
BiocManager::install("MatrixModels")
BiocManager::install("Matrix")
BiocManager::install("sf")
BiocManager::install("quantreg")
BiocManager::install("maptools")
BiocManager::install("ggpmisc")
BiocManager::install("formattable")
BiocManager::install("webshot")
BiocManager::install("splitstackshape")
BiocManager::install("decontam")
BiocManager::install("cowplot")
BiocManager::install("wesanderson")
BiocManager::install("lefser")
BiocManager::install("microbiomeMarker")

install.packages("Matrix")
install.packages("maptools")

####

subset.cecum.physeq.eig.HF.v.NC
subset.cecum.physeq.eig.LF.v.NC
subset.cecum.physeq.eig.HF.v.LF
subset.cecum.physeq.nin.HF.v.NC
subset.cecum.physeq.nin.LF.v.NC
subset.cecum.physeq.nin.HF.v.LF
subset.cecum.physeq.ten.HF.v.NC
subset.cecum.physeq.ten.LF.v.NC
subset.cecum.physeq.ten.HF.v.LF
subset.cecum.physeq.twe.HF.v.NC
subset.cecum.physeq.twe.LF.v.NC
subset.cecum.physeq.twe.HF.v.LF
subset.fecal.physeq.HF.v.NC
subset.fecal.physeq.LF.v.NC
subset.fecal.physeq.HF.v.LF

####

physeq
sample_data(physeq)
otu.relative.table
OTU.rel.table

marker_table(lef_out)
sample_data(lef_out)
marker_table(lef_out)
tax_table(lef_out) 
otu_table(lef_out)

# subset.cecum.physeq.eig.HF.v.NC

lefse_results <- run_lefse(subset.cecum.physeq.eig.HF.v.NC,   group = 'Treatment',   
                           subgroup = NULL,  taxa_rank = "all", transform = c("log10"), norm = "CPM", 
                           norm_para = list(), kw_cutoff = 0.05, lda_cutoff = 2,  bootstrap_n = 30, 
                           bootstrap_fraction = 2/3, wilcoxon_cutoff = 0.05, multigrp_strat = FALSE, strict = c("2"), 
                           sample_min = 10, only_same_subgrp = FALSE, curv = FALSE)

marker_table(lefse_results)

lef_out <- run_lefse(subset.cecum.physeq.eig.HF.v.NC, group = "Treatment", taxa_rank = "Genus", norm = "CPM",
                     kw_cutoff = 0.05, lda_cutoff = 2, strict = c("2"))

plot_cladogram(lefse_results, color=c("orange", "blue"),only_marker = FALSE,
               branch_size = 0.2, alpha = 0.2, node_size_scale = 1, node_size_offset = 1, 
               clade_label_level = 3,clade_label_font_size = 4,  annotation_shape = 22,
               annotation_shape_size = 5,  group_legend_param = list(),  
               marker_legend_param = list())
ggsave("LEFSE_clado_subset.cecum.physeq.eig.HF.v.NC.pdf", height=8, width=10, device="pdf")

plot_cladogram(lef_out, color=c("orange", "blue"),only_marker = FALSE,
               branch_size = 0.2, alpha = 0.2, node_size_scale = 1, node_size_offset = 1, 
               clade_label_level = 3,clade_label_font_size = 4,  annotation_shape = 22,
               annotation_shape_size = 5,  group_legend_param = list(),  
               marker_legend_param = list())

quartz()
plot_ef_bar(lef_out)
ggsave("LEFSE_EFbar_subset.cecum.physeq.eig.HF.v.NC.pdf", height=8, width=10, device="pdf")

plot_abundance(lef_out, label_level=6, max_label_len = 200, group = "Treatment")
ggsave("LEFSE_abundbar_subset.cecum.physeq.eig.HF.v.NC.pdf", height=8, width=10, device="pdf")

# subset.cecum.physeq.eig.LF.v.NC

lefse_results <- run_lefse(subset.cecum.physeq.eig.LF.v.NC,   group = 'Treatment',   
                           subgroup = NULL,  taxa_rank = "all", transform = c("log10"), norm = "CPM", 
                           norm_para = list(), kw_cutoff = 0.05, lda_cutoff = 2,  bootstrap_n = 30, 
                           bootstrap_fraction = 2/3, wilcoxon_cutoff = 0.05, multigrp_strat = FALSE, strict = c("2"), 
                           sample_min = 10, only_same_subgrp = FALSE, curv = FALSE)

lef_out <- run_lefse(subset.cecum.physeq.eig.LF.v.NC, group = "Treatment", taxa_rank = "Genus", norm = "CPM",
                     kw_cutoff = 0.05, lda_cutoff = 2)

plot_cladogram(lefse_results, color=c("orange", "blue"),only_marker = FALSE,
               branch_size = 0.2, alpha = 0.2, node_size_scale = 1, node_size_offset = 1, 
               clade_label_level = 3,clade_label_font_size = 4,  annotation_shape = 22,
               annotation_shape_size = 5,  group_legend_param = list(),  
               marker_legend_param = list())
ggsave("LEFSE_clado_subset.cecum.physeq.eig.LF.v.NC.pdf", height=8, width=20, device="pdf")

plot_ef_bar(lef_out)
ggsave("LEFSE_EFbar_subset.cecum.physeq.eig.LF.v.NC.pdf", height=8, width=10, device="pdf")

plot_abundance(lef_out, label_level=6, max_label_len = 200, group = "Treatment")
ggsave("LEFSE_abundbar_subset.cecum.physeq.eig.LF.v.NC.pdf", height=8, width=10, device="pdf")

# subset.cecum.physeq.eig.HF.v.LF

lefse_results <- run_lefse(subset.cecum.physeq.eig.HF.v.LF,   group = 'Treatment',   
                           subgroup = NULL,  taxa_rank = "all", transform = c("log10"), norm = "CPM", 
                           norm_para = list(), kw_cutoff = 0.05, lda_cutoff = 2,  bootstrap_n = 30, 
                           bootstrap_fraction = 2/3, wilcoxon_cutoff = 0.05, multigrp_strat = FALSE, strict = c("2"), 
                           sample_min = 10, only_same_subgrp = FALSE, curv = FALSE)

lef_out <- run_lefse(subset.cecum.physeq.eig.HF.v.LF, group = "Treatment", taxa_rank = "Genus", norm = "CPM",
                     kw_cutoff = 0.05, lda_cutoff = 2)

plot_cladogram(lefse_results, color=c("orange", "blue"),only_marker = FALSE,
               branch_size = 0.2, alpha = 0.2, node_size_scale = 1, node_size_offset = 1, 
               clade_label_level = 3,clade_label_font_size = 4,  annotation_shape = 22,
               annotation_shape_size = 5,  group_legend_param = list(),  
               marker_legend_param = list())
ggsave("LEFSE_clado_subset.cecum.physeq.eig.HF.v.LF.pdf", height=8, width=20, device="pdf")

plot_ef_bar(lef_out)
ggsave("LEFSE_EFbar_subset.cecum.physeq.eig.HF.v.LF.pdf", height=8, width=10, device="pdf")

plot_abundance(lef_out, label_level=6, max_label_len = 200, group = "Treatment")
ggsave("LEFSE_abundbar_subset.cecum.physeq.eig.HF.v.LF.pdf", height=8, width=10, device="pdf")

# subset.cecum.physeq.nin.HF.v.NC

lefse_results <- run_lefse(subset.cecum.physeq.nin.HF.v.NC,   group = 'Treatment',   
                           subgroup = NULL,  taxa_rank = "all", transform = c("log10"), norm = "CPM", 
                           norm_para = list(), kw_cutoff = 0.05, lda_cutoff = 2,  bootstrap_n = 30, 
                           bootstrap_fraction = 2/3, wilcoxon_cutoff = 0.05, multigrp_strat = FALSE, strict = c("2"), 
                           sample_min = 10, only_same_subgrp = FALSE, curv = FALSE)

lef_out <- run_lefse(subset.cecum.physeq.nin.HF.v.NC, group = "Treatment", taxa_rank = "Genus", norm = "CPM",
                     kw_cutoff = 0.05, lda_cutoff = 2)

plot_cladogram(lefse_results, color=c("orange", "blue"),only_marker = FALSE,
               branch_size = 0.2, alpha = 0.2, node_size_scale = 1, node_size_offset = 1, 
               clade_label_level = 3,clade_label_font_size = 4,  annotation_shape = 22,
               annotation_shape_size = 5,  group_legend_param = list(),  
               marker_legend_param = list())
ggsave("LEFSE_clado_subset.cecum.physeq.nin.HF.v.NC.pdf", height=8, width=10, device="pdf")

plot_ef_bar(lef_out)
ggsave("LEFSE_EFbar_subset.cecum.physeq.nin.HF.v.NC.pdf", height=8, width=10, device="pdf")

plot_abundance(lef_out, label_level=6, max_label_len = 200, group = "Treatment")
ggsave("LEFSE_abundbar_subset.cecum.physeq.nin.HF.v.NC.pdf", height=8, width=10, device="pdf")

# subset.cecum.physeq.nin.LF.v.NC

lefse_results <- run_lefse(subset.cecum.physeq.nin.LF.v.NC,   group = 'Treatment',   
                           subgroup = NULL,  taxa_rank = "all", transform = c("log10"), norm = "CPM", 
                           norm_para = list(), kw_cutoff = 0.05, lda_cutoff = 2,  bootstrap_n = 30, 
                           bootstrap_fraction = 2/3, wilcoxon_cutoff = 0.05, multigrp_strat = FALSE, strict = c("2"), 
                           sample_min = 10, only_same_subgrp = FALSE, curv = FALSE)

lef_out <- run_lefse(subset.cecum.physeq.nin.LF.v.NC, group = "Treatment", taxa_rank = "Genus", norm = "CPM",
                     kw_cutoff = 0.05, lda_cutoff = 2)

plot_cladogram(lefse_results, color=c("orange", "blue"),only_marker = FALSE,
               branch_size = 0.2, alpha = 0.2, node_size_scale = 1, node_size_offset = 1, 
               clade_label_level = 3,clade_label_font_size = 4,  annotation_shape = 22,
               annotation_shape_size = 5,  group_legend_param = list(),  
               marker_legend_param = list())
ggsave("LEFSE_clado_subset.cecum.physeq.nin.LF.v.NC.pdf", height=8, width=20, device="pdf")

plot_ef_bar(lef_out)
ggsave("LEFSE_EFbar_subset.cecum.physeq.nin.LF.v.NC.pdf", height=8, width=10, device="pdf")

plot_abundance(lef_out, label_level=6, max_label_len = 200, group = "Treatment")
ggsave("LEFSE_abundbar_subset.cecum.physeq.nin.LF.v.NC.pdf", height=8, width=10, device="pdf")

# subset.cecum.physeq.nin.HF.v.LF

lefse_results <- run_lefse(subset.cecum.physeq.nin.HF.v.LF,   group = 'Treatment',   
                           subgroup = NULL,  taxa_rank = "all", transform = c("log10"), norm = "CPM", 
                           norm_para = list(), kw_cutoff = 0.05, lda_cutoff = 2,  bootstrap_n = 30, 
                           bootstrap_fraction = 2/3, wilcoxon_cutoff = 0.05, multigrp_strat = FALSE, strict = c("2"), 
                           sample_min = 10, only_same_subgrp = FALSE, curv = FALSE)

lef_out <- run_lefse(subset.cecum.physeq.nin.HF.v.LF, group = "Treatment", taxa_rank = "Genus", norm = "CPM",
                     kw_cutoff = 0.05, lda_cutoff = 2)

plot_cladogram(lefse_results, color=c("orange", "blue"),only_marker = FALSE,
               branch_size = 0.2, alpha = 0.2, node_size_scale = 1, node_size_offset = 1, 
               clade_label_level = 3,clade_label_font_size = 4,  annotation_shape = 22,
               annotation_shape_size = 5,  group_legend_param = list(),  
               marker_legend_param = list())
ggsave("LEFSE_clado_subset.cecum.physeq.nin.HF.v.LF.pdf", height=8, width=20, device="pdf")

plot_ef_bar(lef_out)
ggsave("LEFSE_EFbar_subset.cecum.physeq.nin.HF.v.LF.pdf", height=8, width=10, device="pdf")

plot_abundance(lef_out, label_level=6, max_label_len = 200, group = "Treatment")
ggsave("LEFSE_abundbar_subset.cecum.physeq.nin.HF.v.LF.pdf", height=8, width=10, device="pdf")

# subset.cecum.physeq.ten.HF.v.NC

lefse_results <- run_lefse(subset.cecum.physeq.ten.HF.v.NC,   group = 'Treatment',   
                           subgroup = NULL,  taxa_rank = "all", transform = c("log10"), norm = "CPM", 
                           norm_para = list(), kw_cutoff = 0.05, lda_cutoff = 2,  bootstrap_n = 30, 
                           bootstrap_fraction = 2/3, wilcoxon_cutoff = 0.05, multigrp_strat = FALSE, strict = c("2"), 
                           sample_min = 10, only_same_subgrp = FALSE, curv = FALSE)

lef_out <- run_lefse(subset.cecum.physeq.ten.HF.v.NC, group = "Treatment", taxa_rank = "Genus", norm = "CPM",
                     kw_cutoff = 0.05, lda_cutoff = 2)

plot_cladogram(lefse_results, color=c("orange", "blue"),only_marker = FALSE,
               branch_size = 0.2, alpha = 0.2, node_size_scale = 1, node_size_offset = 1, 
               clade_label_level = 3,clade_label_font_size = 4,  annotation_shape = 22,
               annotation_shape_size = 5,  group_legend_param = list(),  
               marker_legend_param = list())
ggsave("LEFSE_clado_subset.cecum.physeq.ten.HF.v.NC.pdf", height=8, width=10, device="pdf")

plot_ef_bar(lef_out)
ggsave("LEFSE_EFbar_subset.cecum.physeq.ten.HF.v.NC.pdf", height=8, width=10, device="pdf")

plot_abundance(lef_out, label_level=6, max_label_len = 200, group = "Treatment")
ggsave("LEFSE_abundbar_subset.cecum.physeq.ten.HF.v.NC.pdf", height=8, width=10, device="pdf")

# subset.cecum.physeq.ten.LF.v.NC

lefse_results <- run_lefse(subset.cecum.physeq.ten.LF.v.NC,   group = 'Treatment',   
                           subgroup = NULL,  taxa_rank = "all", transform = c("log10"), norm = "CPM", 
                           norm_para = list(), kw_cutoff = 0.05, lda_cutoff = 2,  bootstrap_n = 30, 
                           bootstrap_fraction = 2/3, wilcoxon_cutoff = 0.05, multigrp_strat = FALSE, strict = c("2"), 
                           sample_min = 10, only_same_subgrp = FALSE, curv = FALSE)

lef_out <- run_lefse(subset.cecum.physeq.ten.LF.v.NC, group = "Treatment", taxa_rank = "Genus", norm = "CPM",
                     kw_cutoff = 0.05, lda_cutoff = 2)

plot_cladogram(lefse_results, color=c("orange", "blue"),only_marker = FALSE,
               branch_size = 0.2, alpha = 0.2, node_size_scale = 1, node_size_offset = 1, 
               clade_label_level = 3,clade_label_font_size = 4,  annotation_shape = 22,
               annotation_shape_size = 5,  group_legend_param = list(),  
               marker_legend_param = list())
ggsave("LEFSE_clado_subset.cecum.physeq.ten.LF.v.NC.pdf", height=8, width=20, device="pdf")

plot_ef_bar(lef_out)
ggsave("LEFSE_EFbar_subset.cecum.physeq.ten.LF.v.NC.pdf", height=8, width=10, device="pdf")

plot_abundance(lef_out, label_level=6, max_label_len = 200, group = "Treatment")
ggsave("LEFSE_abundbar_subset.cecum.physeq.ten.LF.v.NC.pdf", height=8, width=10, device="pdf")

# subset.cecum.physeq.ten.HF.v.LF

lefse_results <- run_lefse(subset.cecum.physeq.ten.HF.v.LF,   group = 'Treatment',   
                           subgroup = NULL,  taxa_rank = "all", transform = c("log10"), norm = "CPM", 
                           norm_para = list(), kw_cutoff = 0.05, lda_cutoff = 2,  bootstrap_n = 30, 
                           bootstrap_fraction = 2/3, wilcoxon_cutoff = 0.05, multigrp_strat = FALSE, strict = c("2"), 
                           sample_min = 10, only_same_subgrp = FALSE, curv = FALSE)

lef_out <- run_lefse(subset.cecum.physeq.ten.HF.v.LF, group = "Treatment", taxa_rank = "Genus", norm = "CPM",
                     kw_cutoff = 0.05, lda_cutoff = 2)

plot_cladogram(lefse_results, color=c("orange", "blue"),only_marker = FALSE,
               branch_size = 0.2, alpha = 0.2, node_size_scale = 1, node_size_offset = 1, 
               clade_label_level = 3,clade_label_font_size = 4,  annotation_shape = 22,
               annotation_shape_size = 5,  group_legend_param = list(),  
               marker_legend_param = list())
ggsave("LEFSE_clado_subset.cecum.physeq.ten.HF.v.LF.pdf", height=8, width=20, device="pdf")

plot_ef_bar(lef_out)
ggsave("LEFSE_EFbar_subset.cecum.physeq.ten.HF.v.LF.pdf", height=8, width=10, device="pdf")

plot_abundance(lef_out, label_level=6, max_label_len = 200, group = "Treatment")
ggsave("LEFSE_abundbar_subset.cecum.physeq.ten.HF.v.LF.pdf", height=8, width=10, device="pdf")

# subset.cecum.physeq.twe.HF.v.NC

lefse_results <- run_lefse(subset.cecum.physeq.twe.HF.v.NC,   group = 'Treatment',   
                           subgroup = NULL,  taxa_rank = "all", transform = c("log10"), norm = "CPM", 
                           norm_para = list(), kw_cutoff = 0.05, lda_cutoff = 2,  bootstrap_n = 30, 
                           bootstrap_fraction = 2/3, wilcoxon_cutoff = 0.05, multigrp_strat = FALSE, strict = c("2"), 
                           sample_min = 10, only_same_subgrp = FALSE, curv = FALSE)

lef_out <- run_lefse(subset.cecum.physeq.twe.HF.v.NC, group = "Treatment", taxa_rank = "Genus", norm = "CPM",
                     kw_cutoff = 0.05, lda_cutoff = 2)

plot_cladogram(lefse_results, color=c("orange", "blue"),only_marker = FALSE,
               branch_size = 0.2, alpha = 0.2, node_size_scale = 1, node_size_offset = 1, 
               clade_label_level = 3,clade_label_font_size = 4,  annotation_shape = 22,
               annotation_shape_size = 5,  group_legend_param = list(),  
               marker_legend_param = list())
ggsave("LEFSE_clado_subset.cecum.physeq.twe.HF.v.NC.pdf", height=8, width=10, device="pdf")

plot_ef_bar(lef_out)
ggsave("LEFSE_EFbar_subset.cecum.physeq.twe.HF.v.NC.pdf", height=8, width=10, device="pdf")

plot_abundance(lef_out, label_level=6, max_label_len = 200, group = "Treatment")
ggsave("LEFSE_abundbar_subset.cecum.physeq.twe.HF.v.NC.pdf", height=8, width=10, device="pdf")

# subset.cecum.physeq.twe.LF.v.NC

lefse_results <- run_lefse(subset.cecum.physeq.twe.LF.v.NC,   group = 'Treatment',   
                           subgroup = NULL,  taxa_rank = "all", transform = c("log10"), norm = "CPM", 
                           norm_para = list(), kw_cutoff = 0.05, lda_cutoff = 2,  bootstrap_n = 30, 
                           bootstrap_fraction = 2/3, wilcoxon_cutoff = 0.05, multigrp_strat = FALSE, strict = c("2"), 
                           sample_min = 10, only_same_subgrp = FALSE, curv = FALSE)

lef_out <- run_lefse(subset.cecum.physeq.twe.LF.v.NC, group = "Treatment", taxa_rank = "Genus", norm = "CPM",
                     kw_cutoff = 0.05, lda_cutoff = 2)

plot_cladogram(lefse_results, color=c("orange", "blue"),only_marker = FALSE,
               branch_size = 0.2, alpha = 0.2, node_size_scale = 1, node_size_offset = 1, 
               clade_label_level = 3,clade_label_font_size = 4,  annotation_shape = 22,
               annotation_shape_size = 5,  group_legend_param = list(),  
               marker_legend_param = list())
ggsave("LEFSE_clado_subset.cecum.physeq.twe.LF.v.NC.pdf", height=8, width=20, device="pdf")

plot_ef_bar(lef_out)
ggsave("LEFSE_EFbar_subset.cecum.physeq.twe.LF.v.NC.pdf", height=8, width=10, device="pdf")

plot_abundance(lef_out, label_level=6, max_label_len = 200, group = "Treatment")
ggsave("LEFSE_abundbar_subset.cecum.physeq.twe.LF.v.NC.pdf", height=8, width=10, device="pdf")

# subset.cecum.physeq.twe.HF.v.LF

lefse_results <- run_lefse(subset.cecum.physeq.twe.HF.v.LF,   group = 'Treatment',   
                           subgroup = NULL,  taxa_rank = "all", transform = c("log10"), norm = "CPM", 
                           norm_para = list(), kw_cutoff = 0.05, lda_cutoff = 2,  bootstrap_n = 30, 
                           bootstrap_fraction = 2/3, wilcoxon_cutoff = 0.05, multigrp_strat = FALSE, strict = c("2"), 
                           sample_min = 10, only_same_subgrp = FALSE, curv = FALSE)

lef_out <- run_lefse(subset.cecum.physeq.twe.HF.v.LF, group = "Treatment", taxa_rank = "Genus", norm = "CPM",
                     kw_cutoff = 0.05, lda_cutoff = 2)

plot_cladogram(lefse_results, color=c("orange", "blue"),only_marker = FALSE,
               branch_size = 0.2, alpha = 0.2, node_size_scale = 1, node_size_offset = 1, 
               clade_label_level = 3,clade_label_font_size = 4,  annotation_shape = 22,
               annotation_shape_size = 5,  group_legend_param = list(),  
               marker_legend_param = list())
ggsave("LEFSE_clado_subset.cecum.physeq.twe.HF.v.LF.pdf", height=8, width=20, device="pdf")

plot_ef_bar(lef_out)
ggsave("LEFSE_EFbar_subset.cecum.physeq.twe.HF.v.LF.pdf", height=8, width=10, device="pdf")

plot_abundance(lef_out, label_level=6, max_label_len = 200, group = "Treatment")
ggsave("LEFSE_abundbar_subset.cecum.physeq.twe.HF.v.LF.pdf", height=8, width=10, device="pdf")

### subset.fecal.physeq.HF.v.NC
### subset.fecal.physeq.LF.v.NC
### subset.fecal.physeq.HF.v.LF

# subset.fecal.physeq.HF.v.NC

lefse_results <- run_lefse(subset.fecal.physeq.HF.v.NC,   group = 'Treatment',   
                           subgroup = NULL,  taxa_rank = "all", transform = c("log10"), norm = "CPM", 
                           norm_para = list(), kw_cutoff = 0.05, lda_cutoff = 2,  bootstrap_n = 30, 
                           bootstrap_fraction = 2/3, wilcoxon_cutoff = 0.05, multigrp_strat = FALSE, strict = c("2"), 
                           sample_min = 10, only_same_subgrp = FALSE, curv = FALSE)

lef_out <- run_lefse(subset.fecal.physeq.HF.v.NC, group = "Treatment", taxa_rank = "Genus", norm = "CPM",
                     kw_cutoff = 0.05, lda_cutoff = 2)

plot_cladogram(lefse_results, color=c("orange", "blue"),only_marker = FALSE,
               branch_size = 0.2, alpha = 0.2, node_size_scale = 1, node_size_offset = 1, 
               clade_label_level = 3,clade_label_font_size = 4,  annotation_shape = 22,
               annotation_shape_size = 5,  group_legend_param = list(),  
               marker_legend_param = list())
ggsave("LEFSE_clado_subset.fecal.physeq.HF.v.NC.pdf", height=8, width=10, device="pdf")

plot_ef_bar(lef_out)
ggsave("LEFSE_EFbar_subset.fecal.physeq.HF.v.NC.pdf", height=8, width=10, device="pdf")

plot_abundance(lef_out, label_level=6, max_label_len = 200, group = "Treatment")
ggsave("LEFSE_abundbar_subset.fecal.physeq.HF.v.NC.pdf", height=8, width=10, device="pdf")

# subset.fecal.physeq.LF.v.NC

lefse_results <- run_lefse(subset.fecal.physeq.LF.v.NC,   group = 'Treatment',   
                           subgroup = NULL,  taxa_rank = "all", transform = c("log10"), norm = "CPM", 
                           norm_para = list(), kw_cutoff = 0.05, lda_cutoff = 2,  bootstrap_n = 30, 
                           bootstrap_fraction = 2/3, wilcoxon_cutoff = 0.05, multigrp_strat = FALSE, strict = c("2"), 
                           sample_min = 10, only_same_subgrp = FALSE, curv = FALSE)

lef_out <- run_lefse(subset.fecal.physeq.LF.v.NC, group = "Treatment", taxa_rank = "Genus", norm = "CPM",
                     kw_cutoff = 0.05, lda_cutoff = 2)

plot_cladogram(lefse_results, color=c("orange", "blue"),only_marker = FALSE,
               branch_size = 0.2, alpha = 0.2, node_size_scale = 1, node_size_offset = 1, 
               clade_label_level = 3,clade_label_font_size = 4,  annotation_shape = 22,
               annotation_shape_size = 5,  group_legend_param = list(),  
               marker_legend_param = list())
ggsave("LEFSE_clado_subset.fecal.physeq.LF.v.NC.pdf", height=8, width=20, device="pdf")

plot_ef_bar(lef_out)
ggsave("LEFSE_EFbar_subset.fecal.physeq.LF.v.NC.pdf", height=8, width=10, device="pdf")

plot_abundance(lef_out, label_level=6, max_label_len = 200, group = "Treatment")
ggsave("LEFSE_abundbar_subset.fecal.physeq.LF.v.NC.pdf", height=8, width=10, device="pdf")

# subset.fecal.physeq.HF.v.LF

lefse_results <- run_lefse(subset.fecal.physeq.HF.v.LF,   group = 'Treatment',   
                           subgroup = NULL,  taxa_rank = "all", transform = c("log10"), norm = "CPM", 
                           norm_para = list(), kw_cutoff = 0.05, lda_cutoff = 2,  bootstrap_n = 30, 
                           bootstrap_fraction = 2/3, wilcoxon_cutoff = 0.05, multigrp_strat = FALSE, strict = c("2"), 
                           sample_min = 10, only_same_subgrp = FALSE, curv = FALSE)

lef_out <- run_lefse(subset.fecal.physeq.HF.v.LF, group = "Treatment", taxa_rank = "Genus", norm = "CPM",
                     kw_cutoff = 0.05, lda_cutoff = 2)

plot_cladogram(lefse_results, color=c("orange", "blue"),only_marker = FALSE,
               branch_size = 0.2, alpha = 0.2, node_size_scale = 1, node_size_offset = 1, 
               clade_label_level = 3,clade_label_font_size = 4,  annotation_shape = 22,
               annotation_shape_size = 5,  group_legend_param = list(),  
               marker_legend_param = list())
ggsave("LEFSE_clado_subset.fecal.physeq.HF.v.LF.pdf", height=8, width=20, device="pdf")

plot_ef_bar(lef_out)
ggsave("LEFSE_EFbar_subset.fecal.physeq.HF.v.LF.pdf", height=8, width=10, device="pdf")

plot_abundance(lef_out, label_level=6, max_label_len = 200, group = "Treatment")
ggsave("LEFSE_abundbar_subset.fecal.physeq.HF.v.LF.pdf", height=8, width=10, device="pdf")
