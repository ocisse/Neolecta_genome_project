library(ggplot2)

palette <- c('#ffffff','#0033BB')

data <- read.table("data.txt", header = T, sep = "\t")
d <- qplot(GenomeSizeMb, ProteomeSize, data = data, colour = factor(Phenotype))

#d + layer(geom = "point")
#d + geom_point()
#d + geom_jitter(size = 5)
#d + geom_smooth(size = 0.1)
d + theme_set(theme_bw())


#d + stat_bin2d(bins = 30)
#d + stat_binhex()
#d + stat_density2d()

# new from Greg
dataCM <- read.table("Multicell_Only.txt", header = T, sep = "\t")
d2 <- qplot (Species, data = dataCM, geom = "bar", weight = ProteomeSize, 
             colour = factor(LINEAGE)) 
d2 + theme_set(theme_bw())
d2 + ggtitle("COMPLEX MULTICELLULAR FUNGI") + 
  labs(x="SPECIES", y="PROTEOME SIZE")

 
dataCM2 <- read.table("Multicell_many.txt", header = T, sep = "\t")
d3 <- qplot (Species, data = dataCM2, geom = "bar", weight = ProteomeSize, 
             colour = factor(LINEAGE)) 
d3 + theme_set(theme_bw())
d3 + ggtitle("COMPLEX MULTICELLULAR FUNGI") + 
  labs(x="SPECIES", y="PROTEOME SIZE")


# with more species
dataCM3 <- read.table("Multicell_Only.tab", header = T, sep = "\t")
d4<- qplot (Acronyme, data = dataCM3, geom = "bar", weight = ProteomeSize,
           colour = factor(dataCM3$LINEAGE)) 
dataCM3$NamesOrdered <-  factor(dataCM3$Names,as.character(dataCM3$Names))

d4 <- qplot (NamesOrdered, data = dataCM3, geom = "bar", weight = ProteomeSize,
           colour = factor(LINEAGE), xlab="Species", ylab="Proteome size", title="OK") 
d4 + theme_set(theme_bw())
d4 + coord_flip() 
