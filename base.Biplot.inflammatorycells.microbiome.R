library(vegan)
library(permute)
library(lattice)
library(tidyverse)


#use relative table 

Phylum.rel.table = tax_glom(otu.relative.table, taxrank = "Phylum")
Class.rel.table = tax_glom(otu.relative.table, taxrank = "Class")
Order.rel.table = tax_glom(otu.relative.table, taxrank = "Order")
Family.rel.table = tax_glom(otu.relative.table, taxrank = "Family")
Genus.rel.table = tax_glom(otu.relative.table, taxrank = "Genus")
OTU.rel.table = tax_glom(otu.relative.table, taxrank = "Species")

#use relative table 
#BAL.baseline.table.pruned.prim.rel

otu.relative.table
new_data <- sample_data(otu.relative.table)
new_data$SampleID <- rownames(new_data)
as_tibble(new_data)

#read data 
# new_data <- read_csv(file = "230927.NTM_Master.file.with.NO.MRN.csv")

#filter primary BAL samples only 
# new_data <- new_data %>% filter(Description=="BAL", Primary.lobe!="N.A.")

otu.relative.table@sam_data

#get only necessary variables 
new_data <- new_data %>% select(SampleID, gdT17, CTL, Th17, pgTreg, cTreg, Th1)

#get sample data 
old_data <- data.frame(sample_data(otu.relative.table))

#add neccessary variables to sample data 
old_data <- inner_join(old_data, new_data, by="SampleID")

#add variables to count table 
BAL.baseline.table.pruned.prim@sam_data$NET_Assay_Primary_lobes <- as.numeric(old_data$NET_Assay_Primary_lobes)
BAL.baseline.table.pruned.prim@sam_data$Macrophages <- as.numeric(old_data$Macrophages)
BAL.baseline.table.pruned.prim@sam_data$Neutrophils <- as.numeric(old_data$Neutrophils)
BAL.baseline.table.pruned.prim@sam_data$Eosinophils <- as.numeric(old_data$Eosinophils)
BAL.baseline.table.pruned.prim@sam_data$Lymphocytes <- as.numeric(old_data$Lymphocytes)

#add to relative table 
BAL.baseline.table.pruned.prim.rel@sam_data$NET_Assay_Primary_lobes <- as.numeric(old_data$NET_Assay_Primary_lobes)
BAL.baseline.table.pruned.prim.rel@sam_data$Macrophages <- as.numeric(old_data$Macrophages)
BAL.baseline.table.pruned.prim.rel@sam_data$Neutrophils <- as.numeric(old_data$Neutrophils)
BAL.baseline.table.pruned.prim.rel@sam_data$Eosinophils <- as.numeric(old_data$Eosinophils)
BAL.baseline.table.pruned.prim.rel@sam_data$Lymphocytes <- as.numeric(old_data$Lymphocytes)


#create categories of the varaibles (gdT17, CTL, Th17, pgTreg, cTreg, Th1)
otu.relative.table@sam_data$gdT17_catg <- ifelse(otu.relative.table@sam_data$gdT17 < 
                                                   median(otu.relative.table@sam_data$gdT17, na.rm = TRUE), 
                                                 "Low_gdT17", "High_gdT17")

otu.relative.table@sam_data$CTL_catg <- ifelse(otu.relative.table@sam_data$CTL < 
                                                 median(otu.relative.table@sam_data$CTL, na.rm = TRUE), 
                                               "Low_CTL", "High_CTL")

otu.relative.table@sam_data$Th17_catg <- ifelse(otu.relative.table@sam_data$Th17 < 
                                                  median(otu.relative.table@sam_data$Th17, na.rm = TRUE), 
                                                "Low_Th17", "High_Th17")

otu.relative.table@sam_data$pgTreg_catg <- ifelse(otu.relative.table@sam_data$pgTreg < 
                                                    median(otu.relative.table@sam_data$pgTreg, na.rm = TRUE), 
                                                  "Low_pgTreg", "High_pgTreg")

otu.relative.table@sam_data$cTreg_catg <- ifelse(otu.relative.table@sam_data$cTreg < 
                                                   median(otu.relative.table@sam_data$cTreg, na.rm = TRUE), 
                                                 "Low_cTreg", "High_cTreg")

otu.relative.table@sam_data$Th1_catg <- ifelse(otu.relative.table@sam_data$Th1 < 
                                                 median(otu.relative.table@sam_data$Th1, na.rm = TRUE), 
                                               "Low_Th1", "High_Th1")


#plot ordination from relative table 
ps.rel<- otu.relative.table
ord <- ordinate(ps.rel, "NMDS", "bray")

#plot for gdT17_catg only 
p <- plot_ordination(physeq = ps.rel, ordination = ord, type = "biplot", color = "gdT17_catg", shape = "gdT17_catg")
p <- p+scale_color_manual(values = c("Taxa"="purple", "Low_gdT17"="blue", "High_gdT17"="red"))+
  labs(color = "gdT17_catg",
       shape = "gdT17")+
  theme(panel.background = element_blank(),
        panel.border=element_rect(fill=NA),
        panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "darkgrey"),
        panel.grid.minor = element_blank(),
        strip.background=element_blank(),
        axis.title=element_text(size=30,face="bold"),
        axis.text.x=element_text(colour = "darkgrey", size = rel(0.75)),
        axis.text.y=element_text(colour = "darkgrey", size = rel(0.75)),
        axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), 
        legend.text = element_text(size = 16), 
        legend.title = element_text(size = 20))

#save plot 
pdf(file = "NMDS_biplot_gdT17.pdf", height = 22, width = 20)
p
dev.off()

#repeat for CTL_catg 
p <- plot_ordination(physeq = ps.rel, ordination = ord, type = "biplot", color = "CTL_catg", shape = "CTL_catg")
p <- p+scale_color_manual(values = c("Taxa"="purple", "Low_CTL"="blue", "High_CTL"="red"))+
  labs(color = "CTL_catg",
       shape = "CTL")+
  theme(panel.background = element_blank(),
        panel.border=element_rect(fill=NA),
        panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "darkgrey"),
        panel.grid.minor = element_blank(),
        strip.background=element_blank(),
        axis.title=element_text(size=30,face="bold"),
        axis.text.x=element_text(colour = "darkgrey", size = rel(0.75)),
        axis.text.y=element_text(colour = "darkgrey", size = rel(0.75)),
        axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), 
        legend.text = element_text(size = 16), 
        legend.title = element_text(size = 20))

#save plot 
pdf(file = "NMDS_biplot_CTL.pdf", height = 22, width = 20)
p
dev.off()


#repeat for eosinophils 
p <- plot_ordination(physeq = ps.rel, ordination = ord, type = "biplot", color = "Th17_catg", shape = "Th17_catg")
p <- p+scale_color_manual(values = c("Taxa"="purple", "Low_Th17"="blue", "High_Th17"="red"))+
  labs(color = "Th17_catg",
       shape = "Th17")+
  theme(panel.background = element_blank(),
        panel.border=element_rect(fill=NA),
        panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "darkgrey"),
        panel.grid.minor = element_blank(),
        strip.background=element_blank(),
        axis.title=element_text(size=30,face="bold"),
        axis.text.x=element_text(colour = "darkgrey", size = rel(0.75)),
        axis.text.y=element_text(colour = "darkgrey", size = rel(0.75)),
        axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), 
        legend.text = element_text(size = 16), 
        legend.title = element_text(size = 20))

#save plot 
pdf(file = "NMDS_biplot_Th17.pdf", height = 22, width = 20)
p
dev.off()

#repeat for pgTreg_catg 
p <- plot_ordination(physeq = ps.rel, ordination = ord, type = "biplot", color = "pgTreg_catg", shape = "pgTreg_catg")
p <- p+scale_color_manual(values = c("Taxa"="purple", "Low_pgTreg"="blue", "High_pgTreg"="red"))+
  labs(color = "pgTreg_catg",
       shape = "pgTreg")+
  theme(panel.background = element_blank(),
        panel.border=element_rect(fill=NA),
        panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "darkgrey"),
        panel.grid.minor = element_blank(),
        strip.background=element_blank(),
        axis.title=element_text(size=30,face="bold"),
        axis.text.x=element_text(colour = "darkgrey", size = rel(0.75)),
        axis.text.y=element_text(colour = "darkgrey", size = rel(0.75)),
        axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), 
        legend.text = element_text(size = 16), 
        legend.title = element_text(size = 20))

#save plot 
pdf(file = "NMDS_biplot_pgTreg.pdf", height = 22, width = 20)
p
dev.off()

#repeat for cTreg 
p <- plot_ordination(physeq = ps.rel, ordination = ord, type = "biplot", color = "cTreg_catg", shape = "cTreg_catg")
p <- p+scale_color_manual(values = c("Taxa"="purple", "Low_cTreg"="blue", "High_cTreg"="red"))+
  labs(color = "cTreg_catg",
       shape = "cTreg")+
  theme(panel.background = element_blank(),
        panel.border=element_rect(fill=NA),
        panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "darkgrey"),
        panel.grid.minor = element_blank(),
        strip.background=element_blank(),
        axis.title=element_text(size=30,face="bold"),
        axis.text.x=element_text(colour = "darkgrey", size = rel(0.75)),
        axis.text.y=element_text(colour = "darkgrey", size = rel(0.75)),
        axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), 
        legend.text = element_text(size = 16), 
        legend.title = element_text(size = 20))

#save plot 
pdf(file = "NMDS_biplot_cTreg.pdf", height = 22, width = 20)
p
dev.off()

#repeat for cTreg 
p <- plot_ordination(physeq = ps.rel, ordination = ord, type = "biplot", color = "Th1_catg", shape = "Th1_catg")
p <- p+scale_color_manual(values = c("Taxa"="purple", "Low_Th1"="blue", "High_Th1"="red"))+
  labs(color = "Th1_catg",
       shape = "Th1")+
  theme(panel.background = element_blank(),
        panel.border=element_rect(fill=NA),
        panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "darkgrey"),
        panel.grid.minor = element_blank(),
        strip.background=element_blank(),
        axis.title=element_text(size=30,face="bold"),
        axis.text.x=element_text(colour = "darkgrey", size = rel(0.75)),
        axis.text.y=element_text(colour = "darkgrey", size = rel(0.75)),
        axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), 
        legend.text = element_text(size = 16), 
        legend.title = element_text(size = 20))

#save plot 
pdf(file = "NMDS_biplot_Th1.pdf", height = 22, width = 20)
p
dev.off()



#########
#########
#########
#########
#########

#repeat for Th1 
p <- plot_ordination(physeq = ps.rel, ordination = ord, type = "biplot", color = "Th1_catg", shape = "Treatment")
p <- p+scale_color_manual(values = c("Taxa"="purple", "Low_Th1"="blue", "High_Th1"="red"))+
  labs(color = "Th1_catg",
       shape = "Treatment")+
  scale_shape_manual(values = c(20, 15, 17, 18))+
  theme(panel.background = element_blank(),
        panel.border=element_rect(fill=NA),
        panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "darkgrey"),
        panel.grid.minor = element_blank(),
        strip.background=element_blank(),
        axis.title=element_text(size=30,face="bold"),
        axis.text.x=element_text(colour = "darkgrey", size = rel(0.75)),
        axis.text.y=element_text(colour = "darkgrey", size = rel(0.75)),
        axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), 
        legend.text = element_text(size = 16), 
        legend.title = element_text(size = 20))

#save plot 
pdf(file = "NMDS_biplot_Th1_treatment.pdf", height = 22, width = 20)
p
dev.off()

#plot for gdT17_catg only 
p <- plot_ordination(physeq = ps.rel, ordination = ord, type = "biplot", color = "gdT17_catg", shape = "Treatment")
p <- p+scale_color_manual(values = c("Taxa"="purple", "Low_gdT17"="blue", "High_gdT17"="red"))+
  labs(color = "gdT17_catg",
       shape = "Treatment")+
  scale_shape_manual(values = c(20, 15, 17, 18))+
  theme(panel.background = element_blank(),
        panel.border=element_rect(fill=NA),
        panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "darkgrey"),
        panel.grid.minor = element_blank(),
        strip.background=element_blank(),
        axis.title=element_text(size=30,face="bold"),
        axis.text.x=element_text(colour = "darkgrey", size = rel(0.75)),
        axis.text.y=element_text(colour = "darkgrey", size = rel(0.75)),
        axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), 
        legend.text = element_text(size = 16), 
        legend.title = element_text(size = 20))

#save plot 
pdf(file = "NMDS_biplot_gdT17_Treatment.pdf", height = 22, width = 20)
p
dev.off()

#repeat for CTL_catg 
p <- plot_ordination(physeq = ps.rel, ordination = ord, type = "biplot", color = "CTL_catg", shape = "Treatment")
p <- p+scale_color_manual(values = c("Taxa"="purple", "Low_CTL"="blue", "High_CTL"="red"))+
  labs(color = "CTL_catg",
       shape = "Treatment")+
  scale_shape_manual(values = c(20, 15, 17, 18))+
  theme(panel.background = element_blank(),
        panel.border=element_rect(fill=NA),
        panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "darkgrey"),
        panel.grid.minor = element_blank(),
        strip.background=element_blank(),
        axis.title=element_text(size=30,face="bold"),
        axis.text.x=element_text(colour = "darkgrey", size = rel(0.75)),
        axis.text.y=element_text(colour = "darkgrey", size = rel(0.75)),
        axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), 
        legend.text = element_text(size = 16), 
        legend.title = element_text(size = 20))

#save plot 
pdf(file = "NMDS_biplot_CTL_Treatment.pdf", height = 22, width = 20)
p
dev.off()

#repeat for Th17 
p <- plot_ordination(physeq = ps.rel, ordination = ord, type = "biplot", color = "Th17_catg", shape = "Treatment")
p <- p+scale_color_manual(values = c("Taxa"="purple", "Low_Th17"="blue", "High_Th17"="red"))+
  labs(color = "Th17_catg",
       shape = "Treatment")+
  scale_shape_manual(values = c(20, 15, 17, 18))+
  theme(panel.background = element_blank(),
        panel.border=element_rect(fill=NA),
        panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "darkgrey"),
        panel.grid.minor = element_blank(),
        strip.background=element_blank(),
        axis.title=element_text(size=30,face="bold"),
        axis.text.x=element_text(colour = "darkgrey", size = rel(0.75)),
        axis.text.y=element_text(colour = "darkgrey", size = rel(0.75)),
        axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), 
        legend.text = element_text(size = 16), 
        legend.title = element_text(size = 20))

#save plot 
pdf(file = "NMDS_biplot_Th17_Treatment.pdf", height = 22, width = 20)
p
dev.off()

#repeat for pgTreg_catg 
p <- plot_ordination(physeq = ps.rel, ordination = ord, type = "biplot", color = "pgTreg_catg", shape = "Treatment")
p <- p+scale_color_manual(values = c("Taxa"="purple", "Low_pgTreg"="blue", "High_pgTreg"="red"))+
  labs(color = "pgTreg_catg",
       shape = "Treatment")+
  scale_shape_manual(values = c(20, 15, 17, 18))+
  theme(panel.background = element_blank(),
        panel.border=element_rect(fill=NA),
        panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "darkgrey"),
        panel.grid.minor = element_blank(),
        strip.background=element_blank(),
        axis.title=element_text(size=30,face="bold"),
        axis.text.x=element_text(colour = "darkgrey", size = rel(0.75)),
        axis.text.y=element_text(colour = "darkgrey", size = rel(0.75)),
        axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), 
        legend.text = element_text(size = 16), 
        legend.title = element_text(size = 20))

#save plot 
pdf(file = "NMDS_biplot_pgTreg_Treatment.pdf", height = 22, width = 20)
p
dev.off()

#repeat for cTreg 
p <- plot_ordination(physeq = ps.rel, ordination = ord, type = "biplot", color = "cTreg_catg", shape = "Treatment")
p <- p+scale_color_manual(values = c("Taxa"="purple", "Low_cTreg"="blue", "High_cTreg"="red"))+
  labs(color = "cTreg_catg",
       shape = "Treatment")+
  scale_shape_manual(values = c(20, 15, 17, 18))+
  theme(panel.background = element_blank(),
        panel.border=element_rect(fill=NA),
        panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "darkgrey"),
        panel.grid.minor = element_blank(),
        strip.background=element_blank(),
        axis.title=element_text(size=30,face="bold"),
        axis.text.x=element_text(colour = "darkgrey", size = rel(0.75)),
        axis.text.y=element_text(colour = "darkgrey", size = rel(0.75)),
        axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), 
        legend.text = element_text(size = 16), 
        legend.title = element_text(size = 20))

#save plot 
pdf(file = "NMDS_biplot_cTreg_Treatment.pdf", height = 22, width = 20)
p
dev.off()

#########
#########
#########
#########
#########

ps.rel<- otu.relative.table

tax_table(ps.rel) <- cbind(tax_table(ps.rel), asv=taxa_names(ps.rel))
myranks = c( "Family", "Genus")
mylabels = apply(tax_table(ps.rel)[, myranks], 1, paste, sep="", collapse="_")
mylabels <- paste0(mylabels,paste0(".."), paste0(str_sub(rownames(tax_table(ps.rel)), - 3, - 1)))
tax_table(ps.rel) <- cbind(tax_table(ps.rel), catglab=mylabels)

ord <- ordinate(ps.rel, "NMDS", "bray")


pord <- plot_ordination(physeq = ps.rel, ordination = ord, type = "biplot", color = "Species", shape = "Treatment")

#Put Data into a dataframe
ord <- as.data.frame(pord$data)
#Select only the axes data
df <- subset(ord,select=c("NMDS1","NMDS2"))
#remove NaN or NA 
df <- df[complete.cases(df), ]
#calculated Sample variance for each PC
vars <- apply(df, 2, var)
#Create Variable with the Percent Variance
percentVar <- round(100 * (vars/sum(vars)))

#Rename Variables for PC1 and PC2
colnames(ord)[colnames(ord)=="NMDS1"] <- "PC1"
colnames(ord)[colnames(ord)=="NMDS2"] <- "PC2"

#decide what otu to save 
otu.to.save <-as.character(ord$asv)

#from relative table we should get the mean across the row of the otu table
otu.relative.table.df <- data.frame(otu_table(ps.rel))
otu.relative.table.df.meanRA <- rowMeans(otu.relative.table.df)
#need to subset AND reorder just the otus that we have 
otu.relative.table.df.meanRA.save <- otu.relative.table.df.meanRA[otu.to.save]
#add the abundnace data for the res dataframe
ord$abundance <- otu.relative.table.df.meanRA.save


#Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~Treatment,data= ord, mean)
#Merge the Centroid Data into the PCOA Data
ord <- merge(ord,centroids,by="Treatment",suffixes=c("",".centroid"))

# add contamlist - not going to use this
# December 29, 2023
# ord<- left_join(ord, contamlist, by="asv")

ord$taxa <- ifelse(ord$id.type=="Samples", "",ord$catglab)
ord$taxa2 <- ifelse(grepl("NA_NA.*", ord$taxa), "", ord$taxa)
# ord$contam_color <- ifelse(ord$contaminant=="TRUE", "red", "black")
# ord$contam_color <- factor(ord$contam_color)
label.data <- subset(ord,Treatment=="Taxa")
label.data <- label.data %>% filter(abundance>1.000191e-04)
ord$Treatment <- factor(ord$Treatment)

pdf("Treatment_BiPlot.pdf", height = 15, width = 20)
ggplot(ord, aes(PC1, PC2, color=Treatment)) +
  geom_point(size= ifelse(ord$Treatment=="Taxa", 400 * ord$abundance, 2),alpha=0.7) +    
  geom_point(data=subset(ord,Treatment!="Taxa"),size=5,alpha=0.7) +
  #geom_point(data=subset(ord,FINAL_BRONCHIECTASIS_CULTURE_COMB_Simple_Type=="Taxa"),size= ifelse(FINAL_BRONCHIECTASIS_CULTURE_COMB_Simple_Type=="Taxa", 200 * ord$abundance, 2)) +
  xlab(paste0("NMDS1: ",percentVar[1],"% variance")) +
  ylab(paste0("NMDS2: ",percentVar[2],"% variance")) + 
  #scale_color_manual(values=c("#BEBEBE","#FF5DE7", "#EA3323","#00CED1","#296218")) + 
  #plot point and lines from centroid
  geom_point(data=subset(centroids,Treatment!="Taxa"), aes(x=PC1, y=PC2, color=Treatment), size=0) +
  geom_segment(data=subset(ord,Treatment!="Taxa"), aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color=Treatment))+ 
  #geom_text_repel(aes(label=ifelse(newResults$name%in% c("COPD.0002.BALF.L.untouched","COPD.0030.BAL.L.untouched","COPD.0035.BAL.L.untouched") , as.character(newResults$name),'')),size=3,force=25) +
  #labels centroids 
  #geom_text_repel(aes(label=ifelse(newResults$FINAL_BRONCHIECTASIS_CULTURE_COMB=="Pseudomonas" , as.character(newResults$FINAL_BRONCHIECTASIS_CULTURE_COMB),'')),size=3,force=25) +
  geom_label_repel(data=subset(centroids,Treatment!="Taxa"), aes(x=PC1, y=PC2, label=Treatment), size=10) +
  #geom_label_repel(data = subset(ord, NET_Assay_Primary_lobes_catg=="Taxa"), aes(x=PC1, y=PC2, label=Macrophages_catg, size=4, color=Macrophages_catg))+
  #geom_label_repel(data = subset(ord, NET_Assay_Primary_lobes_catg=="Taxa"), aes(x=PC1, y=PC2, label=Neutrophils_catg, size=4, color=Neutrophils_catg))+
  #geom_label_repel(data = subset(ord, NET_Assay_Primary_lobes_catg=="Taxa"), aes(x=PC1, y=PC2, label=Eosinophils_catg, size=4, color=Eosinophils_catg))+
  #geom_label_repel(data = subset(ord, NET_Assay_Primary_lobes_catg=="Taxa"), aes(x=PC1, y=PC2, label=Lymphocytes_catg, size=4, color=Lymphocytes_catg))+
  scale_color_manual(values = c("black", "goldenrod", "blue", "red", "black"))+
  geom_text_repel(data=label.data,aes(label=taxa2, size=4))+
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), 
        legend.position="none")
dev.off()

#################################################################
#################################################################
#################################################################

# Categories: 
#create categories of the varaibles (gdT17, CTL, Th17, pgTreg, cTreg, Th1)
# gdT17_catg gdT17 Low_gdT17 High_gdT17
# CTL_catg CTL Low_CTL "High_CTL
# Th17_catg Th17 Low_Th17 High_Th17
# pgTreg_catg pgTreg Low_pgTreg High_pgTreg
# cTreg_catg cTreg Low_cTreg High_cTreg
# Th1_catg Th1 Low_Th1 High_Th1"

#repeat for gdT17_catg
ps.rel<- otu.relative.table

tax_table(ps.rel) <- cbind(tax_table(ps.rel), asv=taxa_names(ps.rel))
myranks = c( "Family", "Genus")
mylabels = apply(tax_table(ps.rel)[, myranks], 1, paste, sep="", collapse="_")
mylabels <- paste0(mylabels,paste0(".."), paste0(str_sub(rownames(tax_table(ps.rel)), - 3, - 1)))
tax_table(ps.rel) <- cbind(tax_table(ps.rel), catglab=mylabels)

ord <- ordinate(ps.rel, "NMDS", "bray")


pord <- plot_ordination(physeq = ps.rel, ordination = ord, type = "biplot", color = "Species", shape = "gdT17_catg")

#Put Data into a dataframe
ord <- as.data.frame(pord$data)
#Select only the axes data
df <- subset(ord,select=c("NMDS1","NMDS2"))
#remove NaN or NA 
df <- df[complete.cases(df), ]
#calculated Sample variance for each PC
vars <- apply(df, 2, var)
#Create Variable with the Percent Variance
percentVar <- round(100 * (vars/sum(vars)))

#Rename Variables for PC1 and PC2
colnames(ord)[colnames(ord)=="NMDS1"] <- "PC1"
colnames(ord)[colnames(ord)=="NMDS2"] <- "PC2"

#decide what otu to save 
otu.to.save <-as.character(ord$asv)

#from relative table we should get the mean across the row of the otu table
otu.relative.table.df <- data.frame(otu_table(ps.rel))
otu.relative.table.df.meanRA <- rowMeans(otu.relative.table.df)
#need to subset AND reorder just the otus that we have 
otu.relative.table.df.meanRA.save <- otu.relative.table.df.meanRA[otu.to.save]
#add the abundnace data for the res dataframe
ord$abundance <- otu.relative.table.df.meanRA.save

#Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~gdT17_catg,data= ord, mean)
#Merge the Centroid Data into the PCOA Data
ord <- merge(ord,centroids,by="gdT17_catg",suffixes=c("",".centroid"))

#add contamlist 
#ord<- left_join(ord, contamlist, by="asv")

ord$taxa <- ifelse(ord$id.type=="Samples", "",ord$catglab)
ord$taxa2 <- ifelse(grepl("NA_NA.*", ord$taxa), "", ord$taxa)
# ord$contam_color <- ifelse(ord$contaminant=="TRUE", "red", "black")
# ord$contam_color <- factor(ord$contam_color)
label.data <- subset(ord,gdT17_catg=="Taxa")
label.data <- label.data %>% filter(abundance>1.000191e-04)
ord$gdT17_catg <- factor(ord$gdT17_catg)

#check median 
median(ord$gdT17, na.rm = TRUE)

pdf("gdT17_catg_biplot.pdf", height = 15, width = 20)
ggplot(ord, aes(PC1, PC2, color=gdT17_catg)) +
  geom_point(size= ifelse(ord$gdT17_catg=="Taxa", 400 * ord$abundance, 2),alpha=0.7) +    
  geom_point(data=subset(ord,gdT17_catg!="Taxa"),size=5,alpha=0.7) +
  #geom_point(data=subset(ord,FINAL_BRONCHIECTASIS_CULTURE_COMB_Simple_Type=="Taxa"),size= ifelse(FINAL_BRONCHIECTASIS_CULTURE_COMB_Simple_Type=="Taxa", 200 * ord$abundance, 2)) +
  xlab(paste0("NMDS1: ",percentVar[1],"% variance")) +
  ylab(paste0("NMDS2: ",percentVar[2],"% variance")) + 
  #scale_color_manual(values=c("#BEBEBE","#FF5DE7", "#EA3323","#00CED1","#296218")) + 
  #plot point and lines from centroid
  geom_point(data=subset(centroids,gdT17_catg!="Taxa"), aes(x=PC1, y=PC2, color=gdT17_catg), size=0) +
  geom_segment(data=subset(ord,gdT17_catg!="Taxa"), aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color=gdT17_catg))+ 
  #geom_text_repel(aes(label=ifelse(newResults$name%in% c("COPD.0002.BALF.L.untouched","COPD.0030.BAL.L.untouched","COPD.0035.BAL.L.untouched") , as.character(newResults$name),'')),size=3,force=25) +
  #labels centroids 
  #geom_text_repel(aes(label=ifelse(newResults$FINAL_BRONCHIECTASIS_CULTURE_COMB=="Pseudomonas" , as.character(newResults$FINAL_BRONCHIECTASIS_CULTURE_COMB),'')),size=3,force=25) +
  geom_label_repel(data=subset(centroids,gdT17_catg!="Taxa"), aes(x=PC1, y=PC2, label=gdT17_catg), size=10) +
  #geom_label_repel(data = subset(ord, NET_Assay_Primary_lobes_catg=="Taxa"), aes(x=PC1, y=PC2, label=Macrophages_catg, size=4, color=Macrophages_catg))+
  #geom_label_repel(data = subset(ord, NET_Assay_Primary_lobes_catg=="Taxa"), aes(x=PC1, y=PC2, label=Neutrophils_catg, size=4, color=Neutrophils_catg))+
  #geom_label_repel(data = subset(ord, NET_Assay_Primary_lobes_catg=="Taxa"), aes(x=PC1, y=PC2, label=Eosinophils_catg, size=4, color=Eosinophils_catg))+
  #geom_label_repel(data = subset(ord, NET_Assay_Primary_lobes_catg=="Taxa"), aes(x=PC1, y=PC2, label=Lymphocytes_catg, size=4, color=Lymphocytes_catg))+
  scale_color_manual(values = c("black", "goldenrod", "blue", "red", "black"))+
  geom_text_repel(data=label.data,aes(label=taxa2, size=4))+
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), 
        legend.position="none")
dev.off()

#repeat for CTL_catg
ps.rel<- otu.relative.table

tax_table(ps.rel) <- cbind(tax_table(ps.rel), asv=taxa_names(ps.rel))
myranks = c( "Family", "Genus")
mylabels = apply(tax_table(ps.rel)[, myranks], 1, paste, sep="", collapse="_")
mylabels <- paste0(mylabels,paste0(".."), paste0(str_sub(rownames(tax_table(ps.rel)), - 3, - 1)))
tax_table(ps.rel) <- cbind(tax_table(ps.rel), catglab=mylabels)

ord <- ordinate(ps.rel, "NMDS", "bray")


pord <- plot_ordination(physeq = ps.rel, ordination = ord, type = "biplot", color = "Species", shape = "CTL_catg")

#Put Data into a dataframe
ord <- as.data.frame(pord$data)
#Select only the axes data
df <- subset(ord,select=c("NMDS1","NMDS2"))
#remove NaN or NA 
df <- df[complete.cases(df), ]
#calculated Sample variance for each PC
vars <- apply(df, 2, var)
#Create Variable with the Percent Variance
percentVar <- round(100 * (vars/sum(vars)))

#Rename Variables for PC1 and PC2
colnames(ord)[colnames(ord)=="NMDS1"] <- "PC1"
colnames(ord)[colnames(ord)=="NMDS2"] <- "PC2"

#decide what otu to save 
otu.to.save <-as.character(ord$asv)

#from relative table we should get the mean across the row of the otu table
otu.relative.table.df <- data.frame(otu_table(ps.rel))
otu.relative.table.df.meanRA <- rowMeans(otu.relative.table.df)
#need to subset AND reorder just the otus that we have 
otu.relative.table.df.meanRA.save <- otu.relative.table.df.meanRA[otu.to.save]
#add the abundnace data for the res dataframe
ord$abundance <- otu.relative.table.df.meanRA.save

#Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~CTL_catg,data= ord, mean)
#Merge the Centroid Data into the PCOA Data
ord <- merge(ord,centroids,by="CTL_catg",suffixes=c("",".centroid"))

#add contamlist 
#ord<- left_join(ord, contamlist, by="asv")

ord$taxa <- ifelse(ord$id.type=="Samples", "",ord$catglab)
ord$taxa2 <- ifelse(grepl("NA_NA.*", ord$taxa), "", ord$taxa)
# ord$contam_color <- ifelse(ord$contaminant=="TRUE", "red", "black")
# ord$contam_color <- factor(ord$contam_color)
label.data <- subset(ord,CTL_catg=="Taxa")
label.data <- label.data %>% filter(abundance>1.000191e-04)
ord$CTL_catg <- factor(ord$CTL_catg)

#check median 
median(ord$CTL, na.rm = TRUE)

pdf("CTL_catg_biplot.pdf", height = 15, width = 20)
ggplot(ord, aes(PC1, PC2, color=CTL_catg)) +
  geom_point(size= ifelse(ord$CTL_catg=="Taxa", 400 * ord$abundance, 2),alpha=0.7) +    
  geom_point(data=subset(ord,CTL_catg!="Taxa"),size=5,alpha=0.7) +
  #geom_point(data=subset(ord,FINAL_BRONCHIECTASIS_CULTURE_COMB_Simple_Type=="Taxa"),size= ifelse(FINAL_BRONCHIECTASIS_CULTURE_COMB_Simple_Type=="Taxa", 200 * ord$abundance, 2)) +
  xlab(paste0("NMDS1: ",percentVar[1],"% variance")) +
  ylab(paste0("NMDS2: ",percentVar[2],"% variance")) + 
  #scale_color_manual(values=c("#BEBEBE","#FF5DE7", "#EA3323","#00CED1","#296218")) + 
  #plot point and lines from centroid
  geom_point(data=subset(centroids,CTL_catg!="Taxa"), aes(x=PC1, y=PC2, color=CTL_catg), size=0) +
  geom_segment(data=subset(ord,CTL_catg!="Taxa"), aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color=CTL_catg))+ 
  #geom_text_repel(aes(label=ifelse(newResults$name%in% c("COPD.0002.BALF.L.untouched","COPD.0030.BAL.L.untouched","COPD.0035.BAL.L.untouched") , as.character(newResults$name),'')),size=3,force=25) +
  #labels centroids 
  #geom_text_repel(aes(label=ifelse(newResults$FINAL_BRONCHIECTASIS_CULTURE_COMB=="Pseudomonas" , as.character(newResults$FINAL_BRONCHIECTASIS_CULTURE_COMB),'')),size=3,force=25) +
  geom_label_repel(data=subset(centroids,CTL_catg!="Taxa"), aes(x=PC1, y=PC2, label=CTL_catg), size=10) +
  #geom_label_repel(data = subset(ord, NET_Assay_Primary_lobes_catg=="Taxa"), aes(x=PC1, y=PC2, label=Macrophages_catg, size=4, color=Macrophages_catg))+
  #geom_label_repel(data = subset(ord, NET_Assay_Primary_lobes_catg=="Taxa"), aes(x=PC1, y=PC2, label=Neutrophils_catg, size=4, color=Neutrophils_catg))+
  #geom_label_repel(data = subset(ord, NET_Assay_Primary_lobes_catg=="Taxa"), aes(x=PC1, y=PC2, label=Eosinophils_catg, size=4, color=Eosinophils_catg))+
  #geom_label_repel(data = subset(ord, NET_Assay_Primary_lobes_catg=="Taxa"), aes(x=PC1, y=PC2, label=Lymphocytes_catg, size=4, color=Lymphocytes_catg))+
  scale_color_manual(values = c("black", "goldenrod", "blue", "red", "black"))+
  geom_text_repel(data=label.data,aes(label=taxa2, size=4))+
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), 
        legend.position="none")
dev.off()

#repeat for Th17_catg
ps.rel<- otu.relative.table

tax_table(ps.rel) <- cbind(tax_table(ps.rel), asv=taxa_names(ps.rel))
myranks = c( "Family", "Genus")
mylabels = apply(tax_table(ps.rel)[, myranks], 1, paste, sep="", collapse="_")
mylabels <- paste0(mylabels,paste0(".."), paste0(str_sub(rownames(tax_table(ps.rel)), - 3, - 1)))
tax_table(ps.rel) <- cbind(tax_table(ps.rel), catglab=mylabels)

ord <- ordinate(ps.rel, "NMDS", "bray")


pord <- plot_ordination(physeq = ps.rel, ordination = ord, type = "biplot", color = "Species", shape = "Th17_catg")

#Put Data into a dataframe
ord <- as.data.frame(pord$data)
#Select only the axes data
df <- subset(ord,select=c("NMDS1","NMDS2"))
#remove NaN or NA 
df <- df[complete.cases(df), ]
#calculated Sample variance for each PC
vars <- apply(df, 2, var)
#Create Variable with the Percent Variance
percentVar <- round(100 * (vars/sum(vars)))

#Rename Variables for PC1 and PC2
colnames(ord)[colnames(ord)=="NMDS1"] <- "PC1"
colnames(ord)[colnames(ord)=="NMDS2"] <- "PC2"

#decide what otu to save 
otu.to.save <-as.character(ord$asv)

#from relative table we should get the mean across the row of the otu table
otu.relative.table.df <- data.frame(otu_table(ps.rel))
otu.relative.table.df.meanRA <- rowMeans(otu.relative.table.df)
#need to subset AND reorder just the otus that we have 
otu.relative.table.df.meanRA.save <- otu.relative.table.df.meanRA[otu.to.save]
#add the abundnace data for the res dataframe
ord$abundance <- otu.relative.table.df.meanRA.save

#Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~Th17_catg,data= ord, mean)
#Merge the Centroid Data into the PCOA Data
ord <- merge(ord,centroids,by="Th17_catg",suffixes=c("",".centroid"))

#add contamlist 
#ord<- left_join(ord, contamlist, by="asv")

ord$taxa <- ifelse(ord$id.type=="Samples", "",ord$catglab)
ord$taxa2 <- ifelse(grepl("NA_NA.*", ord$taxa), "", ord$taxa)
# ord$contam_color <- ifelse(ord$contaminant=="TRUE", "red", "black")
# ord$contam_color <- factor(ord$contam_color)
label.data <- subset(ord,Th17_catg=="Taxa")
label.data <- label.data %>% filter(abundance>1.000191e-04)
ord$Th17_catg <- factor(ord$Th17_catg)

#check median 
median(ord$Th17, na.rm = TRUE)

pdf("Th17_catg_biplot.pdf", height = 15, width = 20)
ggplot(ord, aes(PC1, PC2, color=Th17_catg)) +
  geom_point(size= ifelse(ord$Th17_catg=="Taxa", 400 * ord$abundance, 2),alpha=0.7) +    
  geom_point(data=subset(ord,Th17_catg!="Taxa"),size=5,alpha=0.7) +
  #geom_point(data=subset(ord,FINAL_BRONCHIECTASIS_CULTURE_COMB_Simple_Type=="Taxa"),size= ifelse(FINAL_BRONCHIECTASIS_CULTURE_COMB_Simple_Type=="Taxa", 200 * ord$abundance, 2)) +
  xlab(paste0("NMDS1: ",percentVar[1],"% variance")) +
  ylab(paste0("NMDS2: ",percentVar[2],"% variance")) + 
  #scale_color_manual(values=c("#BEBEBE","#FF5DE7", "#EA3323","#00CED1","#296218")) + 
  #plot point and lines from centroid
  geom_point(data=subset(centroids,Th17_catg!="Taxa"), aes(x=PC1, y=PC2, color=Th17_catg), size=0) +
  geom_segment(data=subset(ord,Th17_catg!="Taxa"), aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color=Th17_catg))+ 
  #geom_text_repel(aes(label=ifelse(newResults$name%in% c("COPD.0002.BALF.L.untouched","COPD.0030.BAL.L.untouched","COPD.0035.BAL.L.untouched") , as.character(newResults$name),'')),size=3,force=25) +
  #labels centroids 
  #geom_text_repel(aes(label=ifelse(newResults$FINAL_BRONCHIECTASIS_CULTURE_COMB=="Pseudomonas" , as.character(newResults$FINAL_BRONCHIECTASIS_CULTURE_COMB),'')),size=3,force=25) +
  geom_label_repel(data=subset(centroids,Th17_catg!="Taxa"), aes(x=PC1, y=PC2, label=Th17_catg), size=10) +
  #geom_label_repel(data = subset(ord, NET_Assay_Primary_lobes_catg=="Taxa"), aes(x=PC1, y=PC2, label=Macrophages_catg, size=4, color=Macrophages_catg))+
  #geom_label_repel(data = subset(ord, NET_Assay_Primary_lobes_catg=="Taxa"), aes(x=PC1, y=PC2, label=Neutrophils_catg, size=4, color=Neutrophils_catg))+
  #geom_label_repel(data = subset(ord, NET_Assay_Primary_lobes_catg=="Taxa"), aes(x=PC1, y=PC2, label=Eosinophils_catg, size=4, color=Eosinophils_catg))+
  #geom_label_repel(data = subset(ord, NET_Assay_Primary_lobes_catg=="Taxa"), aes(x=PC1, y=PC2, label=Lymphocytes_catg, size=4, color=Lymphocytes_catg))+
  scale_color_manual(values = c("black", "goldenrod", "blue", "red", "black"))+
  geom_text_repel(data=label.data,aes(label=taxa2, size=4))+
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), 
        legend.position="none")
dev.off()

#repeat for pgTreg_catg
ps.rel<- otu.relative.table

tax_table(ps.rel) <- cbind(tax_table(ps.rel), asv=taxa_names(ps.rel))
myranks = c( "Family", "Genus")
mylabels = apply(tax_table(ps.rel)[, myranks], 1, paste, sep="", collapse="_")
mylabels <- paste0(mylabels,paste0(".."), paste0(str_sub(rownames(tax_table(ps.rel)), - 3, - 1)))
tax_table(ps.rel) <- cbind(tax_table(ps.rel), catglab=mylabels)

ord <- ordinate(ps.rel, "NMDS", "bray")


pord <- plot_ordination(physeq = ps.rel, ordination = ord, type = "biplot", color = "Species", shape = "pgTreg_catg")

#Put Data into a dataframe
ord <- as.data.frame(pord$data)
#Select only the axes data
df <- subset(ord,select=c("NMDS1","NMDS2"))
#remove NaN or NA 
df <- df[complete.cases(df), ]
#calculated Sample variance for each PC
vars <- apply(df, 2, var)
#Create Variable with the Percent Variance
percentVar <- round(100 * (vars/sum(vars)))

#Rename Variables for PC1 and PC2
colnames(ord)[colnames(ord)=="NMDS1"] <- "PC1"
colnames(ord)[colnames(ord)=="NMDS2"] <- "PC2"

#decide what otu to save 
otu.to.save <-as.character(ord$asv)

#from relative table we should get the mean across the row of the otu table
otu.relative.table.df <- data.frame(otu_table(ps.rel))
otu.relative.table.df.meanRA <- rowMeans(otu.relative.table.df)
#need to subset AND reorder just the otus that we have 
otu.relative.table.df.meanRA.save <- otu.relative.table.df.meanRA[otu.to.save]
#add the abundnace data for the res dataframe
ord$abundance <- otu.relative.table.df.meanRA.save

#Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~pgTreg_catg,data= ord, mean)
#Merge the Centroid Data into the PCOA Data
ord <- merge(ord,centroids,by="pgTreg_catg",suffixes=c("",".centroid"))

#add contamlist 
#ord<- left_join(ord, contamlist, by="asv")

ord$taxa <- ifelse(ord$id.type=="Samples", "",ord$catglab)
ord$taxa2 <- ifelse(grepl("NA_NA.*", ord$taxa), "", ord$taxa)
# ord$contam_color <- ifelse(ord$contaminant=="TRUE", "red", "black")
# ord$contam_color <- factor(ord$contam_color)
label.data <- subset(ord,pgTreg_catg=="Taxa")
label.data <- label.data %>% filter(abundance>1.000191e-04)
ord$pgTreg_catg <- factor(ord$pgTreg_catg)

#check median 
median(ord$pgTreg, na.rm = TRUE)

pdf("pgTreg_catg_biplot.pdf", height = 15, width = 20)
ggplot(ord, aes(PC1, PC2, color=pgTreg_catg)) +
  geom_point(size= ifelse(ord$pgTreg_catg=="Taxa", 400 * ord$abundance, 2),alpha=0.7) +    
  geom_point(data=subset(ord,pgTreg_catg!="Taxa"),size=5,alpha=0.7) +
  #geom_point(data=subset(ord,FINAL_BRONCHIECTASIS_CULTURE_COMB_Simple_Type=="Taxa"),size= ifelse(FINAL_BRONCHIECTASIS_CULTURE_COMB_Simple_Type=="Taxa", 200 * ord$abundance, 2)) +
  xlab(paste0("NMDS1: ",percentVar[1],"% variance")) +
  ylab(paste0("NMDS2: ",percentVar[2],"% variance")) + 
  #scale_color_manual(values=c("#BEBEBE","#FF5DE7", "#EA3323","#00CED1","#296218")) + 
  #plot point and lines from centroid
  geom_point(data=subset(centroids,pgTreg_catg!="Taxa"), aes(x=PC1, y=PC2, color=pgTreg_catg), size=0) +
  geom_segment(data=subset(ord,pgTreg_catg!="Taxa"), aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color=pgTreg_catg))+ 
  #geom_text_repel(aes(label=ifelse(newResults$name%in% c("COPD.0002.BALF.L.untouched","COPD.0030.BAL.L.untouched","COPD.0035.BAL.L.untouched") , as.character(newResults$name),'')),size=3,force=25) +
  #labels centroids 
  #geom_text_repel(aes(label=ifelse(newResults$FINAL_BRONCHIECTASIS_CULTURE_COMB=="Pseudomonas" , as.character(newResults$FINAL_BRONCHIECTASIS_CULTURE_COMB),'')),size=3,force=25) +
  geom_label_repel(data=subset(centroids,pgTreg_catg!="Taxa"), aes(x=PC1, y=PC2, label=pgTreg_catg), size=10) +
  #geom_label_repel(data = subset(ord, NET_Assay_Primary_lobes_catg=="Taxa"), aes(x=PC1, y=PC2, label=Macrophages_catg, size=4, color=Macrophages_catg))+
  #geom_label_repel(data = subset(ord, NET_Assay_Primary_lobes_catg=="Taxa"), aes(x=PC1, y=PC2, label=Neutrophils_catg, size=4, color=Neutrophils_catg))+
  #geom_label_repel(data = subset(ord, NET_Assay_Primary_lobes_catg=="Taxa"), aes(x=PC1, y=PC2, label=Eosinophils_catg, size=4, color=Eosinophils_catg))+
  #geom_label_repel(data = subset(ord, NET_Assay_Primary_lobes_catg=="Taxa"), aes(x=PC1, y=PC2, label=Lymphocytes_catg, size=4, color=Lymphocytes_catg))+
  scale_color_manual(values = c("black", "goldenrod", "blue", "red", "black"))+
  geom_text_repel(data=label.data,aes(label=taxa2, size=4))+
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), 
        legend.position="none")
dev.off()

#repeat for cTreg_catg
ps.rel<- otu.relative.table

tax_table(ps.rel) <- cbind(tax_table(ps.rel), asv=taxa_names(ps.rel))
myranks = c( "Family", "Genus")
mylabels = apply(tax_table(ps.rel)[, myranks], 1, paste, sep="", collapse="_")
mylabels <- paste0(mylabels,paste0(".."), paste0(str_sub(rownames(tax_table(ps.rel)), - 3, - 1)))
tax_table(ps.rel) <- cbind(tax_table(ps.rel), catglab=mylabels)

ord <- ordinate(ps.rel, "NMDS", "bray")


pord <- plot_ordination(physeq = ps.rel, ordination = ord, type = "biplot", color = "Species", shape = "cTreg_catg")

#Put Data into a dataframe
ord <- as.data.frame(pord$data)
#Select only the axes data
df <- subset(ord,select=c("NMDS1","NMDS2"))
#remove NaN or NA 
df <- df[complete.cases(df), ]
#calculated Sample variance for each PC
vars <- apply(df, 2, var)
#Create Variable with the Percent Variance
percentVar <- round(100 * (vars/sum(vars)))

#Rename Variables for PC1 and PC2
colnames(ord)[colnames(ord)=="NMDS1"] <- "PC1"
colnames(ord)[colnames(ord)=="NMDS2"] <- "PC2"

#decide what otu to save 
otu.to.save <-as.character(ord$asv)

#from relative table we should get the mean across the row of the otu table
otu.relative.table.df <- data.frame(otu_table(ps.rel))
otu.relative.table.df.meanRA <- rowMeans(otu.relative.table.df)
#need to subset AND reorder just the otus that we have 
otu.relative.table.df.meanRA.save <- otu.relative.table.df.meanRA[otu.to.save]
#add the abundnace data for the res dataframe
ord$abundance <- otu.relative.table.df.meanRA.save

#Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~cTreg_catg,data= ord, mean)
#Merge the Centroid Data into the PCOA Data
ord <- merge(ord,centroids,by="cTreg_catg",suffixes=c("",".centroid"))

#add contamlist 
#ord<- left_join(ord, contamlist, by="asv")

ord$taxa <- ifelse(ord$id.type=="Samples", "",ord$catglab)
ord$taxa2 <- ifelse(grepl("NA_NA.*", ord$taxa), "", ord$taxa)
# ord$contam_color <- ifelse(ord$contaminant=="TRUE", "red", "black")
# ord$contam_color <- factor(ord$contam_color)
label.data <- subset(ord,cTreg_catg=="Taxa")
label.data <- label.data %>% filter(abundance>1.000191e-04)
ord$cTreg_catg <- factor(ord$cTreg_catg)

#check median 
median(ord$cTreg, na.rm = TRUE)

pdf("cTreg_catg_biplot.pdf", height = 15, width = 20)
ggplot(ord, aes(PC1, PC2, color=cTreg_catg)) +
  geom_point(size= ifelse(ord$cTreg_catg=="Taxa", 400 * ord$abundance, 2),alpha=0.7) +    
  geom_point(data=subset(ord,cTreg_catg!="Taxa"),size=5,alpha=0.7) +
  #geom_point(data=subset(ord,FINAL_BRONCHIECTASIS_CULTURE_COMB_Simple_Type=="Taxa"),size= ifelse(FINAL_BRONCHIECTASIS_CULTURE_COMB_Simple_Type=="Taxa", 200 * ord$abundance, 2)) +
  xlab(paste0("NMDS1: ",percentVar[1],"% variance")) +
  ylab(paste0("NMDS2: ",percentVar[2],"% variance")) + 
  #scale_color_manual(values=c("#BEBEBE","#FF5DE7", "#EA3323","#00CED1","#296218")) + 
  #plot point and lines from centroid
  geom_point(data=subset(centroids,cTreg_catg!="Taxa"), aes(x=PC1, y=PC2, color=cTreg_catg), size=0) +
  geom_segment(data=subset(ord,cTreg_catg!="Taxa"), aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color=cTreg_catg))+ 
  #geom_text_repel(aes(label=ifelse(newResults$name%in% c("COPD.0002.BALF.L.untouched","COPD.0030.BAL.L.untouched","COPD.0035.BAL.L.untouched") , as.character(newResults$name),'')),size=3,force=25) +
  #labels centroids 
  #geom_text_repel(aes(label=ifelse(newResults$FINAL_BRONCHIECTASIS_CULTURE_COMB=="Pseudomonas" , as.character(newResults$FINAL_BRONCHIECTASIS_CULTURE_COMB),'')),size=3,force=25) +
  geom_label_repel(data=subset(centroids,cTreg_catg!="Taxa"), aes(x=PC1, y=PC2, label=cTreg_catg), size=10) +
  #geom_label_repel(data = subset(ord, NET_Assay_Primary_lobes_catg=="Taxa"), aes(x=PC1, y=PC2, label=Macrophages_catg, size=4, color=Macrophages_catg))+
  #geom_label_repel(data = subset(ord, NET_Assay_Primary_lobes_catg=="Taxa"), aes(x=PC1, y=PC2, label=Neutrophils_catg, size=4, color=Neutrophils_catg))+
  #geom_label_repel(data = subset(ord, NET_Assay_Primary_lobes_catg=="Taxa"), aes(x=PC1, y=PC2, label=Eosinophils_catg, size=4, color=Eosinophils_catg))+
  #geom_label_repel(data = subset(ord, NET_Assay_Primary_lobes_catg=="Taxa"), aes(x=PC1, y=PC2, label=Lymphocytes_catg, size=4, color=Lymphocytes_catg))+
  scale_color_manual(values = c("black", "goldenrod", "blue", "red", "black"))+
  geom_text_repel(data=label.data,aes(label=taxa2, size=4))+
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), 
        legend.position="none")
dev.off()

#repeat for Th1_catg
ps.rel<- otu.relative.table

tax_table(ps.rel) <- cbind(tax_table(ps.rel), asv=taxa_names(ps.rel))
myranks = c( "Family", "Genus")
mylabels = apply(tax_table(ps.rel)[, myranks], 1, paste, sep="", collapse="_")
mylabels <- paste0(mylabels,paste0(".."), paste0(str_sub(rownames(tax_table(ps.rel)), - 3, - 1)))
tax_table(ps.rel) <- cbind(tax_table(ps.rel), catglab=mylabels)

ord <- ordinate(ps.rel, "NMDS", "bray")


pord <- plot_ordination(physeq = ps.rel, ordination = ord, type = "biplot", color = "Species", shape = "Th1_catg")

#Put Data into a dataframe
ord <- as.data.frame(pord$data)
#Select only the axes data
df <- subset(ord,select=c("NMDS1","NMDS2"))
#remove NaN or NA 
df <- df[complete.cases(df), ]
#calculated Sample variance for each PC
vars <- apply(df, 2, var)
#Create Variable with the Percent Variance
percentVar <- round(100 * (vars/sum(vars)))

#Rename Variables for PC1 and PC2
colnames(ord)[colnames(ord)=="NMDS1"] <- "PC1"
colnames(ord)[colnames(ord)=="NMDS2"] <- "PC2"

#decide what otu to save 
otu.to.save <-as.character(ord$asv)

#from relative table we should get the mean across the row of the otu table
otu.relative.table.df <- data.frame(otu_table(ps.rel))
otu.relative.table.df.meanRA <- rowMeans(otu.relative.table.df)
#need to subset AND reorder just the otus that we have 
otu.relative.table.df.meanRA.save <- otu.relative.table.df.meanRA[otu.to.save]
#add the abundnace data for the res dataframe
ord$abundance <- otu.relative.table.df.meanRA.save

#Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~Th1_catg,data= ord, mean)
#Merge the Centroid Data into the PCOA Data
ord <- merge(ord,centroids,by="Th1_catg",suffixes=c("",".centroid"))

#add contamlist 
#ord<- left_join(ord, contamlist, by="asv")

ord$taxa <- ifelse(ord$id.type=="Samples", "",ord$catglab)
ord$taxa2 <- ifelse(grepl("NA_NA.*", ord$taxa), "", ord$taxa)
# ord$contam_color <- ifelse(ord$contaminant=="TRUE", "red", "black")
# ord$contam_color <- factor(ord$contam_color)
label.data <- subset(ord,Th1_catg=="Taxa")
label.data <- label.data %>% filter(abundance>1.000191e-04)
ord$Th1_catg <- factor(ord$Th1_catg)

#check median 
median(ord$Th1, na.rm = TRUE)

pdf("Th1_catg_biplot.pdf", height = 15, width = 20)
ggplot(ord, aes(PC1, PC2, color=Th1_catg)) +
  geom_point(size= ifelse(ord$Th1_catg=="Taxa", 400 * ord$abundance, 2),alpha=0.7) +    
  geom_point(data=subset(ord,Th1_catg!="Taxa"),size=5,alpha=0.7) +
  #geom_point(data=subset(ord,FINAL_BRONCHIECTASIS_CULTURE_COMB_Simple_Type=="Taxa"),size= ifelse(FINAL_BRONCHIECTASIS_CULTURE_COMB_Simple_Type=="Taxa", 200 * ord$abundance, 2)) +
  xlab(paste0("NMDS1: ",percentVar[1],"% variance")) +
  ylab(paste0("NMDS2: ",percentVar[2],"% variance")) + 
  #scale_color_manual(values=c("#BEBEBE","#FF5DE7", "#EA3323","#00CED1","#296218")) + 
  #plot point and lines from centroid
  geom_point(data=subset(centroids,Th1_catg!="Taxa"), aes(x=PC1, y=PC2, color=Th1_catg), size=0) +
  geom_segment(data=subset(ord,Th1_catg!="Taxa"), aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color=Th1_catg))+ 
  #geom_text_repel(aes(label=ifelse(newResults$name%in% c("COPD.0002.BALF.L.untouched","COPD.0030.BAL.L.untouched","COPD.0035.BAL.L.untouched") , as.character(newResults$name),'')),size=3,force=25) +
  #labels centroids 
  #geom_text_repel(aes(label=ifelse(newResults$FINAL_BRONCHIECTASIS_CULTURE_COMB=="Pseudomonas" , as.character(newResults$FINAL_BRONCHIECTASIS_CULTURE_COMB),'')),size=3,force=25) +
  geom_label_repel(data=subset(centroids,Th1_catg!="Taxa"), aes(x=PC1, y=PC2, label=Th1_catg), size=10) +
  #geom_label_repel(data = subset(ord, NET_Assay_Primary_lobes_catg=="Taxa"), aes(x=PC1, y=PC2, label=Macrophages_catg, size=4, color=Macrophages_catg))+
  #geom_label_repel(data = subset(ord, NET_Assay_Primary_lobes_catg=="Taxa"), aes(x=PC1, y=PC2, label=Neutrophils_catg, size=4, color=Neutrophils_catg))+
  #geom_label_repel(data = subset(ord, NET_Assay_Primary_lobes_catg=="Taxa"), aes(x=PC1, y=PC2, label=Eosinophils_catg, size=4, color=Eosinophils_catg))+
  #geom_label_repel(data = subset(ord, NET_Assay_Primary_lobes_catg=="Taxa"), aes(x=PC1, y=PC2, label=Lymphocytes_catg, size=4, color=Lymphocytes_catg))+
  scale_color_manual(values = c("black", "goldenrod", "blue", "red", "black"))+
  geom_text_repel(data=label.data,aes(label=taxa2, size=4))+
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), 
        legend.position="none")
dev.off()
