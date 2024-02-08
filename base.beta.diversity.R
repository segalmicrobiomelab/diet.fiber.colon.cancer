rownames(sample_data(otu.relative.table))
colnames(sample_data(otu.relative.table))

# From metadata columsn 
sample_data(otu.relative.table)$Sample_type
sample_data(otu.relative.table)$Source
sample_data(otu.relative.table)$Treatment

sample_data(otu.relative.table)$Experiment

### Vegan
### Prints out the available distance methods under UniFrac
### In thise case we are using "wunifrac" for weighted UniFrac 
dist_methods <- unlist(distanceMethodList)
print(dist_methods)

sample_data(otu.relative.table)$Source

wUnif.dist = UniFrac(otu.relative.table, weighted=TRUE, normalized=TRUE)

### estimate number of axes
wUniF.pco = dudi.pco(cailliez(wUnif.dist), scannf = FALSE, nf = 3)

### Plot PCoA for all Strep pneumonaie samples
pdf(file="beta.diversity.samples.sampletype.exposure.pdf", width=6, height=6)
s.class(wUniF.pco $li, interaction(sample_data(otu.relative.table)$Source), col=c("forestgreen", "brown", "red3", "dodgerblue", "red", "dodgerblue", "forestgreen", "blue4", "goldenrod", "orange", "green", "goldenrod4", "brown", "dodgerblue4", "brown4", "yellow", "forestgreen", "pink", "blue", "purple","black", "gray"))
dev.off()

vegan::adonis2(wUnif.dist ~ Source, data=data.frame(sample_data(otu.relative.table)))

Permutation test for adonis under reduced model
Terms added sequentially (first to last)
Permutation: free
Number of permutations: 999

vegan::adonis2(formula = wUnif.dist ~ Source, data = data.frame(sample_data(otu.relative.table)))

#GET AXIS
all.ordinate.wUniF <- ordinate(otu.relative.table, method="PCoA", distance="wUniFrac")

# Scatter plot
pdf(file="ordination.pdf")
plot_ordination(otu.relative.table, all.ordinate.wUniF, color = "Source", shape = "Source")
dev.off()

#################################################################
#################################################################
#################################################################

otu.relative.table
save.image(file="savefile.RData")
load(file="savefile.RData")