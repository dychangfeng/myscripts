
GC <- read.table(file = "~/ncbi_genomes/ncbi-genomes-representive1200/GC_length_two.txt")
all_G4 <- read.table(file = "~/ncbi_genomes/ncbi-genomes-representive1200/all_G4.txt")
spare_G4 <- read.table(file = "~/ncbi_genomes/ncbi-genomes-representive1200/spare_G4.txt")
library(readr)
##--------change from factor to character----------
all_G4$V2 <- as.character(all_G4$V2)
spare_G4$V2 <- as.character(spare_G4$V2)

## n <- sapply(all_G4$V2, function(x) unlist(strsplit(x, '[_]'))[1])
for (i in seq(1:1464)){
  all_G4[i,2] <- unlist(strsplit(all_G4[i,2], '[_]'))[2]
}
for (i in seq(1:1464)){
  spare_G4[i,2] <- unlist(strsplit(spare_G4[i,2], '[_]'))[2]
}
spare_G4$V2 == all_G4$V2 ## now make sure the order are the same in both file

colnames(GC) <- c('GC_percentage', 'genome_length')
colnames(all_G4) <- c('all.g4', 'file')
GC$all.g4 <- all_G4[1:1464, 1]
GC$spare.g4 <- spare_G4[1:1464, 1]
GC$file.name <- all_G4[1:1464, 2]
length(unique(GC$file.name)) ## make sure file.name is unique identifier
library('tidyr')
library('dplyr')
library('ggplot2')

bacteria.1400 <- GC%>%mutate(genome.mbp = genome_length/1000000, GC.density = all.g4/genome.mbp, spare.percentage = 100*spare.g4/all.g4)
write.table(bacteria.1400, file = "~/ncbi_genomes/bacteria_1464.txt")
b1 <- ggplot(data = bacteria_1464_with_info, 
             aes(x=GC_percentage, y=GC_density, size = spare_percentage, fill = proteins_density)) +
  geom_point(shape=21, colour="#000000") +
  scale_fill_continuous(low = 'blue', high = 'red')+
  facet_wrap(~Group, nrow = 3)+
  ggtitle('G4 density over GC percentage for bacteria') +
  labs(x = 'GC Percentage', y = 'GC Density')+
  
  facet_wrap(~Group, nrow = 5)
b2 <- b1 + scale_y_log10() 
bacteria_GC.density_over500 <- filter(bacteria.1400, GC.density>=500)
ggplot(filter(bacteria.1400, GC.density>=200), 
       aes(x=GC_percentage, y=GC.density, size = genome.mbp, fill = spare.percentage)) +
  geom_point(shape=21, colour="#000000") +
  scale_fill_discrete(low = 'blue', high = 'red')+
  ggtitle('G4 density over GC percentage for bacteria') +
  labs(x = 'GC Percentage', y = 'G4 Density')+
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20,face = "bold"),
        plot.title = element_text(size = 20, face = 'bold'))

bacteria_117_with_info <- read.csv("~/ncbi_genomes/bacteria_117_with_info.csv")
#bacteria_1464_with_info <- read.csv("~/ncbi_genomes/bacteria_1464_with_info.csv")
#----------------------------------------------
