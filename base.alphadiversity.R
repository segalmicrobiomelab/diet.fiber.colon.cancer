### Calculate alpha diversity 
### Calculate diversity
### Re-evaluate the alpha diversity 
### Alpha diversity 

sessionInfo()
detach(igraph)
detach(package:igraph)
sessionInfo()

### ALPHA DIVERSITY SHANNON 
### otu.relative.table

### All Samples Alpha diversity 
Shannon_diversity_cecum = vegan::diversity(otu_table(cecum.otu.relative.table), index = "shannon", MARGIN = 2, base = exp(1))

### Plot the diversity 
pdf(file="adiversity_shannon_cecum.treatment.pdf", width=5, height=5)
boxplot(Shannon_diversity_cecum ~ sample_data(cecum.otu.relative.table)$Treatment)
dev.off()

write.table(Shannon_diversity_cecum, file="Shannon_cecum.treatment.txt", sep="\t")

pdf(file="adiversity_shannon_cecum.treatment.pdf", width=3, height=5)
boxplot(Shannon_diversity_cecum ~ sample_data(cecum.otu.relative.table)$Treatment, col=MyColour, outline=FALSE)
stripchart(Shannon_diversity_cecum ~ sample_data(cecum.otu.relative.table)$Treatment, vertical=TRUE, add=TRUE, method="jitter", pch=20, col="black")
dev.off()

### ALPHA DIVERSITY SHANNON 
### otu.relative.table

rarefy_even_depth(fecal.otu.relative.table, sample.size = min(sample_sums(fecal.otu.relative.table)),
                  
                  ### All Samples Alpha diversity 
                  Shannon_diversity_fecal = vegan::diversity(otu_table(fecal.otu.relative.table), index = "shannon", MARGIN = 2, base = exp(1))
                  
                  ### Plot the diversity 
                  pdf(file="adiversity_shannon_fecal.treatment.pdf", width=5, height=5)
                  boxplot(Shannon_diversity_fecal ~ sample_data(fecal.otu.relative.table)$Treatment)
                  dev.off()
                  
                  write.table(Shannon_diversity_cecum, file="Shannon_fecal.treatment.txt", sep="\t")
                  
                  pdf(file="adiversity_shannon_fecal.treatment.pdf", width=3, height=5)
                  boxplot(Shannon_diversity_fecal ~ sample_data(fecal.otu.relative.table)$Treatment, col=MyColour, outline=FALSE)
                  stripchart(Shannon_diversity_fecal ~ sample_data(fecal.otu.relative.table)$Treatment, vertical=TRUE, add=TRUE, method="jitter", pch=20, col="black")
                  dev.off()
                  
                  pdf(file="adiversity_shannon_fecal.treatment.pdf", width=3, height=5)
                  boxplot(Shannon_diversity_fecal ~ sample_data(fecal.otu.relative.table)$Treatment, col=MyColour.new, outline=FALSE)
                  stripchart(Shannon_diversity_fecal ~ sample_data(fecal.otu.relative.table)$Treatment, vertical=TRUE, add=TRUE, method="jitter", pch=20, col="black")
                  dev.off()
                  
                  ### ALPHA DIVERSITY SHANNON 
                  ### otu.relative.table
                  
                  physeq
                  
                  cecum.physeq = subset_samples(physeq, Source %in% c('cecum'))
                  fecal.physeq = subset_samples(physeq, Source %in% c('fecal'))
                  
                  cecum.physeq.rare <- rarefy_even_depth(cecum.physeq, sample.size = min(sample_sums(cecum.physeq)))
                  fecal.physeq.rare <- rarefy_even_depth(fecal.physeq, sample.size = min(sample_sums(fecal.physeq)))
                  
                  
                  ### All Samples Alpha diversity 
                  Shannon_diversity_fecal = vegan::diversity(otu_table(fecal.physeq.rare), index = "shannon", MARGIN = 2, base = exp(1))
                  
                  ### Plot the diversity 
                  pdf(file="adiversity_shannon_fecal.treatment.rare.pdf", width=5, height=5)
                  boxplot(Shannon_diversity_fecal ~ sample_data(fecal.physeq.rare)$Treatment)
                  dev.off()
                  
                  write.table(Shannon_diversity_fecal, file="Shannon_fecal.treatment.rare.txt", sep="\t")
                  
                  pdf(file="adiversity_shannon_fecal.treatment.rare.pdf", width=3, height=5)
                  boxplot(Shannon_diversity_fecal ~ sample_data(fecal.physeq.rare)$Treatment, col=MyColour, outline=FALSE)
                  stripchart(Shannon_diversity_fecal ~ sample_data(fecal.physeq.rare)$Treatment, vertical=TRUE, add=TRUE, method="jitter", pch=20, col="black")
                  dev.off()
                  
                  ### All Samples Alpha diversity 
                  Shannon_diversity_cecum = vegan::diversity(otu_table(cecum.physeq.rare), index = "shannon", MARGIN = 2, base = exp(1))
                  
                  ### Plot the diversity 
                  pdf(file="adiversity_shannon_cecum.treatment.rare.pdf", width=5, height=5)
                  boxplot(Shannon_diversity_cecum ~ sample_data(cecum.physeq.rare)$Treatment)
                  dev.off()
                  
                  write.table(Shannon_diversity_cecum, file="Shannon_cecum.treatment.rare.txt", sep="\t")
                  
                  pdf(file="adiversity_shannon_cecum.treatment.rare.pdf", width=3, height=5)
                  boxplot(Shannon_diversity_cecum ~ sample_data(cecum.physeq.rare)$Treatment, col=MyColour, outline=FALSE)
                  stripchart(Shannon_diversity_cecum ~ sample_data(cecum.physeq.rare)$Treatment, vertical=TRUE, add=TRUE, method="jitter", pch=20, col="black")
                  dev.off()