GC.2 <- read.table(file = "~/ncbi_genomes/ncbi-genomes-2017-08-10/GC_length2.txt")
all_G4.2 <- read.table(file = "~/ncbi_genomes/ncbi-genomes-2017-08-10/all_G4.txt")
spare_G4.2 <- read.table(file = "~/ncbi_genomes/ncbi-genomes-2017-08-10/spare_G4.txt")

##--------change from factor to character----------
all_G4.2$V2 <- as.character(all_G4.2$V2)
spare_G4.2$V2 <- as.character(spare_G4.2$V2)

n <- sapply(all_G4$V2, function(x) unlist(strsplit(x, '[.]'))[1])
for (i in seq(1:1464)){
  all_G4[i,2] <- unlist(strsplit(all_G4[i,2], '[.]'))[1]
}
for (i in seq(1:1464)){
  spare_G4[i,2] <- unlist(strsplit(spare_G4[i,2], '[.]'))[1]
}
spare_G4$V2 == all_G4$V2 ## now make sure the order are the same in both file

colnames(GC.2) <- c('GC_percentage', 'genome_length')
colnames(all_G4.2) <- c('all.g4', 'file')
GC.2$all.g4 <- all_G4.2[1:117, 1]
GC.2$spare.g4 <- spare_G4.2[1:117, 1]
GC.2$file.name <- all_G4.2[1:117, 2]
length(unique(GC.2$file.name)) ## make sure file.name is unique identifier
library('tidyr')
library('dplyr')
library('ggplot2')

bacteria.117.ref <- GC.2%>%mutate(genome.mbp = genome_length/1000000, GC.density = all.g4/genome.mbp, spare.percentage = 100*spare.g4/all.g4)
write.table(bacteria.117.ref, file = "~/ncbi_genomes/bacteria_117_ref.txt")
b1 <- ggplot(data = bacteria.117.ref, 
             aes(x=GC_percentage, y=GC.density, size = genome.mbp, fill = spare.percentage)) +
  geom_point(shape=21, colour="#000000") +
  scale_fill_continuous(low = 'blue', high = 'red')+
  ggtitle('G4 density over GC percentage for bacteria') +
  labs(x = 'GC Percentage', y = 'GC Density')
b2 <- b1 + scale_y_log10() 
