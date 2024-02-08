### Code with all samples: 

#Create Distance Matrix with Bray (or wUniFrac depending what you are using)
vegdist = vegdist(t(otu_table(cecum.otu.relative.table)), method = "bray")

#Formulate principal component co-ordinates for PCOA plot, k as the choice of PCs
CmdScale <- cmdscale(vegdist, k =10)
#calculated Sample variance for each PC
vars <- apply(CmdScale, 2, var)
#Create Variable with the Percent Variance
percentVar <- round(100 * (vars/sum(vars)))

#Merge PC Data with MetaData
require(data.table)
newResults <- merge(x = CmdScale, y = sample_data(cecum.otu.relative.table), by = "row.names", all.x = TRUE)
#Rename Variables for PC1 and PC2
colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"

#Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~ Treatment,data= newResults, mean) #Here you would use your grouping variable (e.g., days)
#Merge the Centroid Data into the PCOA Data
newResults <- merge(newResults,centroids,by="Treatment",suffixes=c("",".centroid"))#Here you would use your grouping variable (e.g., days)

MyColour <- c("#BFB5D8", "#AED9A3", "#A3B8E4") 
names(MyColour) <- c("chow", "low fiber", "high fiber")

pdf("231010_Murine_cecum_treatment.pdf", height = 8, width = 8)
ggplot(newResults, aes(PC1, PC2, color=Treatment)) + # Graph PC1 and PC2
  geom_point(size=5) + # Set the size of the points
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + #Label PC1
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + #Label PC2
  #Set colors for each category, should be the same number of categories
  scale_color_manual(values=c("#A3B8E4", "#AED9A3", "#BFB5D8" )) +
  #plot point and lines from centroid
  geom_point(data=centroids, aes(x=PC1, y=PC2, color=Treatment), size=0) +#Here you would use your grouping variable (e.g., days)
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= Treatment)) + #Here you would use your grouping variable (e.g., days)
  #If you want to identify specific samples use the code bellow
  #geom_text_repel(aes(label=ifelse(newResults$name%in% c("COPD.0002.BAL.L.untouched","COPD.0030.BAL.L.untouched","COPD.0035.BAL.L.untouched") , as.character(newResults$name),'')),size=3,force=25) +
  #labels centroids should be same number of categories
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2), label=c('high fiber', 'low fiber', 'chow'), size=8) + #Here you can label the way you want
  #Use the code bellow if you want to switch the X axis around
  #scale_x_reverse() +
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "grey80", size = rel(0.75)),axis.text.y=element_text(colour = "grey80", size = rel(0.75)),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
dev.off()

ggplot(newResults, aes(PC1, PC2, color=Treatment)) + # Graph PC1 and PC2
  geom_point(size=5) + # Set the size of the points
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + #Label PC1
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + #Label PC2
  #Set colors for each category, should be the same number of categories
  scale_color_manual(values=MyColour) +
  #plot point and lines from centroid
  geom_point(data=centroids, aes(x=PC1, y=PC2, color=Treatment), size=0) +#Here you would use your grouping variable (e.g., days)
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= Treatment)) + #Here you would use your grouping variable (e.g., days)
  #If you want to identify specific samples use the code bellow
  #geom_text_repel(aes(label=ifelse(newResults$name%in% c("COPD.0002.BAL.L.untouched","COPD.0030.BAL.L.untouched","COPD.0035.BAL.L.untouched") , as.character(newResults$name),'')),size=3,force=25) +
  #labels centroids should be same number of categories
  # geom_label_repel(data = centroids, aes(x=PC1, y=PC2), labels=Exposure), size=10) + #Here you can label the way you want
  #Use the code bellow if you want to switch the X axis around
  #scale_x_reverse() +
  theme()

sample_data(cecum.otu.relative.table)
write.table(sample_data(cecum.otu.relative.table), file="export.map.txt", sep="\t")

vegan::adonis2(vegdist ~ Treatment, data=data.frame(sample_data(cecum.otu.relative.table)))

vegan::adonis2(formula = vegdist ~ Treatment, data = data.frame(sample_data(cecum.otu.relative.table)))
