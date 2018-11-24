
library(dplyr)
library(tidyr)
library(ggplot2)
dt_GC <- read.delim('~/ncbi_genomes/Deinococcus_thermus/Deinococcus_thermus.txt', stringsAsFactors =FALSE) ## use the file with removed white spaces
dt_G4 <- read.delim('~/ncbi_genomes/Deinococcus_thermus/all_G4.3.txt', header = FALSE)
dt_spare_G4 <- read.delim('~/ncbi_genomes/Deinococcus_thermus/spare_G4.3.txt', header = FALSE)
dt_G4 <- dt_G4[1:83,]%>%separate('V2', into=c('GCA','GCA_number','ASM','n1','n2'), 
                                 sep="_", remove=TRUE)%>%unite('GCA_access',c('GCA','GCA_number'))
dt_spare_G4 <- dt_spare_G4[1:83,]%>%separate('V2', into=c('n1','n2','n3','n4','n5'),
                                             sep='_', remove=TRUE)%>%unite('GCA_access',c('n1','n2'))
dt_GC['GCA_access'] <- as.character(dt_GC$Assembly.Accession)
colnames(dt_G4) <- c('all_G4','GCA_access', 'n1')
colnames(dt_spare_G4) <- c('spare_G4', 'GCA_access','n1')
G4 <- merge(dt_G4,dt_spare_G4,by='GCA_access')
G4 <- G4[c('GCA_access','all_G4','spare_G4')]
denicoccus_thermus_all <- merge(G4,dt_GC, by='GCA_access')
denicoccus_thermus_simple <- denicoccus_thermus_all[c('GCA_access', 'all_G4',
                                                      'spare_G4', 'X.Organism.Name',
                                                      'Size..Mb.','GC.', 'Genes',
                                                      'Proteins')]
colnames(denicoccus_thermus_simple) <- c('GCA_access','all_G4','spare_G4','Name',
                                         'genome_size_mb','GC','Genes','Proteins')
denicoccus_thermus_simple <- denicoccus_thermus_simple%>%mutate(G4_density=all_G4/genome_size_mb, spare_percentage=
                                     spare_G4*100/all_G4)
denicoccus_thermus_simple$Proteins <- as.numeric(denicoccus_thermus_simple$Proteins)
denicoccus_thermus_simple$GC <- as.numeric(denicoccus_thermus_simple$GC)
denicoccus_thermus_simple$G4_density <- as.numeric(denicoccus_thermus_simple$G4_density)
denicoccus_thermus_simple <- denicoccus_thermus_simple%>%mutate(protein_density=Proteins/genome_size_mb)
write.table(denicoccus_thermus_simple, file = "~/ncbi_genomes/Deinococcus_thermus/D_T_with_genus_G4.txt",
            sep='\t', row.names = FALSE)
denicoccus_thermus_simple = read.table("~/Documents/bacteria_G4/D_thermus/D_T_with_genus_G4.txt", sep='\t', header = TRUE)

dt_plot=ggplot(denicoccus_thermus_simple)+ geom_point(shape=21, color='#000001',
                                                 aes(x=GC,y=G4_density,
                                                     size=spare_percentage,
                                                     fill=protein_density))+
  scale_fill_gradient(high = "red", low = "blue")+
  ggtitle('G4 density over GC percentage \n for Denicoccus Thermus')+
  theme(axis.text = element_text(size = 18),axis.title = element_text(size = 20,face='bold'),plot.title = element_text(size = 20, face = 'bold'))
  
ggsave('~/Desktop/denicoccus_thermus_all.png', plot = dt_plot, width = 10, height = 6,
       units = 'in')
denicoccus_thermus_simple%>%filter(grepl('Thermus', Name))
denicoccus_thermus_simple=denicoccus_thermus_simple%>%separate('Name', into=c('genus','species','n1','n2'), 
                                     sep=" ", remove=FALSE)

dt_genus=ggplot(denicoccus_thermus_simple)+ geom_point(shape=21, color='#000001',
                                                      aes(x=GC,y=G4_density,
                                                          size=spare_percentage,
                                                          fill=genus))+
  scale_fill_hue()+
  ggtitle('PQS density over GC percentage for Denicoccus Thermus')+
  labs(x='%GC', y='PQS density/Mbp')+
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 20,face='bold'),
        plot.title = element_text(size = 20, face = 'bold'),
        legend.text = element_text(size = 18, face = 'bold'),
        legend.title = element_text(size=18, face='bold'))+
  guides(color=guide_legend(override.aes = list(size=10)))

ggsave('~/Desktop/denicoccus_thermus_genus.png', plot = dt_genus, width = 10, height = 5
      ,units = 'in')
random_G4 <- read.csv("~/Desktop/G4_count_random_2m_20_80_2.csv")
colnames(random_G4) <- seq(18,78,2)
random_G4$`18` <- NULL ## remove first column
random_G4 <- random_G4%>%gather() ## gather all the column names as key , values in one columns
random_G4$key <- as.numeric(random_G4$key)
##-------------------fit the random curve with a function--------------
fit <- lm(random_G4$key~random_G4$value) ## should only put one y for one x to fit
plot(random_G4$key, fit$fitted.values)

##-------------------------------------------
dt_genum_with_control = dt_genus + geom_boxplot(data = random_G4, aes(x=key, y=value, group=cut_width(key,2)), color ='black')
ggsave('~/Desktop/denicoccus_thermus_genus_with_control.png', plot = dt_genum_with_control, width = 10, height = 5
       ,units = 'in')
library(svglite)
ggsave('~/Desktop/denicoccus_thermus.svg', plot = dt_genum_with_control, width = 12, height = 8,
       units = 'in')
filter(denicoccus_thermus_simple, genus=='Deinococcus')%>%select(c('Name','spare_percentage', "G4_density"))
ggplot(denicoccus_thermus_simple%>%filter(genus=='Deinococcus'))+ geom_point(shape=21, color='#000001',
                                              aes(x=GC,y=G4_density,
                                                  size=spare_percentage,
                                                  fill=protein_density))+
  scale_fill_gradient(high = "red", low = "blue")+
  ggtitle('G4 density over GC percentage \n for Denicoccus Thermus')+
  theme(axis.text = element_text(size = 18),axis.title = element_text(size = 20,face='bold'),plot.title = element_text(size = 20, face = 'bold'))
