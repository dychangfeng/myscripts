library('ggplot2')
library('readr')
devtools::install_github("davidgohel/gdtools")
install.packages('gdtools')
install.packages('svglite')
library(svglite)
test_genomes <- read_tsv("/Users/Yun/Documents/bacteria_G4/test_genome/test_genomes.txt", col_names = TRUE)
colnames(test_genomes)[c(4,8)] <- c("spare_percentage","genome_length_log10")

#test_genomes['Species'] <- c('H_pylori', 'Actinomyces', 'Clostridium','D_radiodurans', 'Ecoli', 'Archea', 'yeast', 'araTha', 'Chlamy', 'Rice', 'Chicken', 'Maize', 'Trout', 'Alligator', 'Dog', 'Cat', 'Mouse', 'Human')

test_plot=ggplot(test_genomes)+ geom_point(shape=21, color='#000001',
                                       aes(x=genome_length_log10,y=G4_density,
                                           size=spare_percentage,
                                           fill=GC))+
  labs(x='Genome length (logMbp)', y='PQS density')+
  geom_text(aes(x=genome_length_log10, y=G4_density, label=Species), size=6,nudge_x = -0.15, vjust=0, nudge_y = 0.5)+
  scale_fill_continuous(high = "red", low = "blue")+
  ggtitle('PQS density over genomes length for test_genomes')+
  theme(axis.text = element_text(size = 18),axis.title = element_text(size = 20,face='bold'),plot.title = element_text(size = 20, face = 'bold'))

ggsave('~/Desktop/test_genome_plog2.svg', plot = test_plot, width = 12, height = 8,
       units = 'in')
