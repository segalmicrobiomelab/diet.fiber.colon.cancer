metabolites <- read_csv("Negative alignment_updated.csv")
matrix <- dist(metabolites)
matrix

cmdscale(matrix)



attach(metabolites)
View(metabolites)

metabolites$group

rownames(metabolites) <- metabolites$group
metabolites$group <- NULL
metabolites$Name <- NULL

LPL$SampleID <- NULL
Genus_II$SampleID <- NULL


#In R, the Hellinger-transform is performed using the decostand() function
metabolites_hellinger <- decostand(metabolites, method = "hellinger")

# Predicting transformed species using all variables contained in LPL
my_rda <- rda(metabolites_hellinger ~ ., data = metabolites)

summary(my_rda)

#show the triplot
quartz()
plot(my_rda)
plot(my_rda,display=c("sp","bp"))

#Type 1: scaling 1 shows similarities between objects in the response matrix. 
#. Sites (numbers) that are closer together have more similar communities; 
#. Species that are closer together occupy more sites in common.

# Type 2: scaling 2 shows the effects of explanatory variables.

plot(my_rda, display=c("sp","bp"), scaling = 1, type = "text", frame = FALSE)

plot(my_rda, display=c("sp","bp"), scaling = 1, type = "none", frame = FALSE)
plot(my_rda, display=c("sp","bp"), choices = c(1,3),scaling = 1, type = "text", frame = FALSE)

plot(my_rda, display=c("sp","bp"), scaling = 1, type = "text", frame = FALSE)

points(Genus_II, pch = 21, col = "black", bg = "steelblue", cex = 1.2)
text(Genus_II, pch = 22, col = 'black', bg = "#f2bd33", cex = 1.2)
plot(my_rda, scaling = 2, type = "text", main = "my_rda - Scaling 2")




par(mfrow=c(1,1))
plot(my_rda, scaling=1,display=c("sp","bp"))
var.sc <- scores(my_rda, choices=1:2, scaling=1, display="sp")
arrows(0, 0, var.sc[,1], var.sc[, 2], length = 0, lty = 1, col = "blue")
plot(my_rda, choices = c(1,3),display=c("sp","bp")) 

# Custom triplot code
### extract % explained by the first 2 axes'
quartz()
pdf(file="metabolite_RDA_plot.pdf", width=8, height=8)
plot(my_rda, display=c("sp","bp"), scaling = 1, type = "none", frame = TRUE)

perc <- round(100*(summary(my_rda)$cont$importance[2, 1:2]), 2)
sc_si <- scores(my_rda, display="sites", choices=c(1,2), scaling=1)
sc_sp <- scores(my_rda, display="species", choices=c(1,2), scaling=1)
sc_bp <- scores(my_rda, display="bp", choices=c(1, 2), scaling=1)

# add points for site scores
points(sc_si, pch = 21, col = "black", bg = "steelblue", cex = 1.2)
# add points for species scores
points(sc_sp, pch = 22, col = "black", bg = "#f2bd33", cex = 1.2)
# add text labels for species abbreviations
text(sc_sp + c(0.03, 0.09), labels = rownames(sc_sp), col = "grey40", font = 2, cex = 0.6)
# add arrows for effects of the expanatory variables
arrows(0,0, sc_bp[,1], sc_bp[,2], col = "red", lwd = 2)
# add text labels for arrows
text(x = sc_bp[,1] -0.1,
     y = sc_bp[,2] - 0.03, 
     labels = rownames(sc_bp), 
     col = "red", 
     cex = 1, 
     font = 2)
dev.off()

pdf(file="metabolite_RDA_plot_large.pdf", width=20, height=20)
plot(my_rda, display=c("sp","bp"), scaling = 1, type = "none", frame = TRUE)

perc <- round(100*(summary(my_rda)$cont$importance[2, 1:2]), 2)
sc_si <- scores(my_rda, display="sites", choices=c(1,2), scaling=1)
sc_sp <- scores(my_rda, display="species", choices=c(1,2), scaling=1)
sc_bp <- scores(my_rda, display="bp", choices=c(1, 2), scaling=1)

# add points for site scores
points(sc_si, pch = 21, col = "black", bg = "steelblue", cex = 1.2)
# add points for species scores
points(sc_sp, pch = 22, col = "black", bg = "#f2bd33", cex = 1.2)
# add text labels for species abbreviations
text(sc_sp + c(0.03, 0.09), labels = rownames(sc_sp), col = "grey40", font = 2, cex = 0.6)
# add arrows for effects of the expanatory variables
arrows(0,0, sc_bp[,1], sc_bp[,2], col = "red", lwd = 2)
# add text labels for arrows
text(x = sc_bp[,1] -0.1,
     y = sc_bp[,2] - 0.03, 
     labels = rownames(sc_bp), 
     col = "red", 
     cex = 1, 
     font = 2)
dev.off()

#RDA: explained variance
RsquareAdj(my_rda)

#Test of significance
# Permutation test. For an RDA, you have to test for three different things: 
# 1, Global RDA significance; 2, Axis significance; 3, Terms (explanatory variables) significance
# Global significanec
anova.cca(my_rda, permutations = 999)
# Axis significance
anova.cca(my_rda, by = "axis")
# Terms significance
anova.cca(my_rda, by = "terms")

#Linear dependencies. As a rule of the thumb, if square of VIF >2 multicollinearity is considered high.
sqrt(vif.cca(my_rda))


#presence-absence transformation to calculate species number per sit
Genus_pa <- decostand(Genus_II, "pa")
#calculate sum per species
Genus_sum <- apply(Genus_pa,2,sum)
sort(Genus_sum)

#remove species that occur at less than 3 sites
Genus_fin <- Genus_II[ , ! Genus_sum<3]

sort(apply(Genus_fin, 2, max))
# strong differences in order of magnitude of species abundances
sort(apply(Genus_fin, 2, sd))

ggsave("metabolite_RDA_plot.pdf", height=8, width=10, device="pdf")
dev.off()


ggplot(res.diagdds.pruned.genus.cecum.physeq.twe.HF.v.LF, aes(x = log2FoldChange, y = sig,label=Row.names)) +
  geom_point(color=cols.res.diagdds.pruned.genus.cecum.physeq.twe.HF.v.LF, size=ifelse(res.diagdds.pruned.genus.cecum.physeq.twe.HF.v.LF$log2FoldChange>=1 & res.diagdds.pruned.genus.cecum.physeq.twe.HF.v.LF$padj < alpha, 100 * res.diagdds.pruned.genus.cecum.physeq.twe.HF.v.LF$abundance.experiment, ifelse(res.diagdds.pruned.genus.cecum.physeq.twe.HF.v.LF$log2FoldChange<=-1 & res.diagdds.pruned.genus.cecum.physeq.twe.HF.v.LF$padj < alpha, 100 * res.diagdds.pruned.genus.cecum.physeq.twe.HF.v.LF$abundance.control, 2)),alpha=0.7) +
  geom_text_repel(aes(label=ifelse(res.diagdds.pruned.genus.cecum.physeq.twe.HF.v.LF$padj < alpha & res.diagdds.pruned.genus.cecum.physeq.twe.HF.v.LF$Genus!="g__", as.character(res.diagdds.pruned.genus.cecum.physeq.twe.HF.v.LF$Genus),'')),size=4,force=25, segment.colour="grey",segment.alpha=0.2) +
  geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") +
  xlab("Effect size: log2(fold-change)") +
  ylab("-log10(adjusted p-value)") +
  theme
dev.off()