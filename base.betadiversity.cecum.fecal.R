sample_data(fecal.otu.relative.table)$Treatment
sample_data(fecal.otu.relative.table)$Time.point

### Plot PCoA for all Strep pneumonaie samples
pdf(file="231006.wUniFrac.fecal.treatment.pdf", width=3, height=3)
s.class(wUniF.fecal.pco $li, interaction(sample_data(fecal.otu.relative.table)$Treatment), col=c("purple", "violet","orange"))
dev.off()

vegan::adonis2(wUnif.fecal.dist ~ Treatment, data=data.frame(sample_data(fecal.otu.relative.table)))

Permutation test for adonis under reduced model
Terms added sequentially (first to last)
Permutation: free
Number of permutations: 999

vegan::adonis2(formula = wUnif.fecal.dist ~ Treatment, data = data.frame(sample_data(fecal.otu.relative.table)))
Df SumOfSqs      R2      F Pr(>F)    
Treatment  2  0.16171 0.91109 56.362  0.001 ***
  Residual  11  0.01578 0.08891                  
Total     13  0.17749 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#GET AXIS
all.fecal.ordinate.wUniF <- ordinate(fecal.otu.relative.table, method="PCoA", distance="wUniFrac")

# Scatter plot
pdf(file="w.UniFrac.fecal.treatment.ordinate.pdf")
plot_ordination(fecal.otu.relative.table, all.fecal.ordinate.wUniF, color = "Treatment", shape = "Treatment")
dev.off()

### Plot PCoA for all Strep pneumonaie samples
pdf(file="wUniFrac.fecal.Time.point.pdf", width=3, height=3)
s.class(wUniF.fecal.pco $li, interaction(sample_data(fecal.otu.relative.table)$Time.point), col=c("purple", "violet","orange"))
dev.off()

vegan::adonis2(wUnif.fecal.dist ~ Time.point, data=data.frame(sample_data(fecal.otu.relative.table)))

#GET AXIS
all.fecal.ordinate.wUniF <- ordinate(fecal.otu.relative.table, method="PCoA", distance="wUniFrac")

# Scatter plot
pdf(file="w.UniFrac.fecal.treatment.ordinate.pdf")
plot_ordination(fecal.otu.relative.table, all.fecal.ordinate.wUniF, color = "Treatment", shape = "Treatment")
dev.off()
