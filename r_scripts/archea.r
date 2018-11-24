library(readr)
library(tidyr)
library(dplyr)
library(devtools)
library(ggbiplot) ## to make a pretty ggbiplot
library(ggplot2)
library(tidyr)
library(mosaic)
library(stringi)
library(stringr)
library(pvclust) ## calculate p values for clustering
library(reshape2)
archea_allg4 <- read_tsv('~/Documents/archaea/all_G4_2.txt', col_names = c('all_G4_counts', 'file_name'))
archea_spareg4 <- read_tsv("~/Documents/archaea/spare_G4_2.txt", col_names = c('spare_G4_counts', 'file_name'))
archea_GC <- read_tsv("~/Documents/archaea/GC_genome_length_2.txt", col_names = c('GC_percentage', 'genome_L'))
archea_meta <- read_tsv("~/Documents/archaea/archea_complete_genome.txt", col_names = TRUE)

## scree plot
PVEplot <- function(x){qplot(c(1:10), x[1:10]) + 
    geom_line() + 
    xlab("Principal Component") + 
    ylab("PVE") +
    ggtitle("Scree Plot") +
    ylim(0, 1)}
## biplot
my_theme = theme(legend.direction = 'horizontal', legend.position = 'top',
                 axis.text = element_text(size = 14),
                 axis.title = element_text(size = 20,face = "bold"),
                 plot.title = element_text(size = 20, face = 'bold'),
                 strip.text.x = element_text(size=10, face="bold"),
                 strip.background = element_rect(colour="burlywood", fill="burlywood"),
                 legend.text = element_text(size = 16, face = 'bold' ),
                 legend.title = element_text(size=16, face='bold'))
gbiplot <- function(pca.a, groups){ggbiplot(pca.a, obs.scale = 1, var.scale = 1, 
                                            ellipse = TRUE, groups = groups,
                                            circle = TRUE, inherit.aes = FALSE) + scale_color_discrete(name = '') 
}

## function to calculate pve
pve = function(pca.b) {pca.b$sdev^2/sum(pca.b$sdev^2)}
save_path = '/Users/Yun/Documents/archaea/figures/'
norm = function(x){
  x/sum(x)
}
##-------------get accession number out---------
archea_allg4 <- archea_allg4[1:184,]%>%separate(file_name, into=c('GCA','GCA_number','ASM','n1','n2'), sep="_", remove=TRUE)%>%unite('GCA_access',c('GCA','GCA_number'))%>%select(all_G4_counts, GCA_access)
archea_spareg4 <- archea_spareg4[1:184,]%>%separate(file_name, into=c('GCA','GCA_number','ASM','n1','n2'), sep="_", remove=TRUE)%>%unite('GCA_access',c('GCA','GCA_number'))%>%select(spare_G4_counts, GCA_access)
archea_g4 <- G4 <- merge(archea_allg4,archea_spareg4,by='GCA_access')

archea_meta['GCA_access']=archea_meta$`Assembly Accession`
archea_meta <- archea_meta%>%select("Organism/Name", "Group","SubGroup", "Size (Mb)","GC%","GCA_access","Reference", "Genes",
                                    "Proteins" )
archea_meta['genome_size']=archea_meta$`Size (Mb)`
archea_meta['GC'] =archea_meta$`GC%`
## get the genus from the name of the archea
archea_meta$genus = gsub(' .+', '', archea_meta$`Organism/Name`)
archea <- merge(archea_meta, archea_g4, by="GCA_access")%>%mutate(G4_density=all_G4_counts/genome_size, 
                                                                  protein_density = Proteins/genome_size,
                                                                spare_percentage=spare_G4_counts*100/all_G4_counts)
write_delim(archea, '~/Documents/archaea/archaea_info.txt', delim = '\t')

g4 = read_delim('~/Documents/archaea/G4_bins_tss.txt', delim = '\t')
# remove the rows with G4 number of 0
g4 = g4[!apply(g4[,2:21], 1, function(x){sum(x) == 0}),]
colnames(g4) = c('GCA_access',seq(-285,285, 30))
g4[,2:21] = as.data.frame(t(apply(g4[2:21], 1, norm)))
g4_Subgroup = merge(g4, archea[c('GCA_access', 'SubGroup', 'genus')],by.x = 'GCA_access', by.y = 'GCA_access')

melt_archaea =melt(g4_Subgroup[,2:23], id.vars = c('SubGroup', 'genus'),
     variable.name = 'distance_to_tss')
ggplot(melt_archaea, aes(x=as.factor(distance_to_tss), y = value)) + geom_boxplot(aes(fill = SubGroup)) + 
  my_theme + scale_y_continuous(limits = c(0,0.3)) +
  facet_wrap(~SubGroup, ncol = 2) + 
  labs(x= 'Distance to TSS', y = 'G4 fraction around TSS')
major_groups = melt_archaea$SubGroup %in% c('Methanomicrobia', 'Thermococci','Halobacteria', 'Archaeoglobi')
ggplot(melt_archaea[major_groups,], 
       aes(x=as.factor(distance_to_tss), y = value)) + geom_violin(aes(fill = SubGroup)) + 
  my_theme + 
  facet_wrap(~SubGroup, ncol = 2) + 
  labs(x= 'Distance to TSS', y = 'G4 fraction around TSS')
ggsave('violin plot of archaea major groups.png', path = save_path, scale = 1, width = 10,
       height = 8, units='in')
##-----------PCA analysis------------------
rownames(g4) = g4$X1
pca.ar = prcomp(g4[,2:21],scale = TRUE)
pca.ar$sdev
scores.ar = data.frame(pca.ar$x) ## all PCs
pve.ar = pve(pca.ar)
PVEplot(pve.ar)
ar.bi = gbiplot(pca.ar, g4_Subgroup$SubGroup) + my_theme
ar.bi


random_G4 <- read.csv("~/Documents/bacteria_G4/G4_count_random_2m_20_80_2.csv")

colnames(random_G4) <- seq(18,78,2)
random_G4$`18` <- NULL ## remove first column
random_G4 <- random_G4%>%gather() ## gather all the column names as key , values in one columns
random_G4$key <- as.numeric(random_G4$key)

archea_plot=ggplot(archea)+ geom_point(shape=21, color='#000001',
                                       aes(x=GC,y=G4_density))+
  scale_fill_hue()+
  facet_wrap(~SubGroup, nrow = 3)+xlab('GC%') + ylab('PQS density') + 
  ggtitle('PQS density over GC percentage for Archea')+my_theme

archea_control=archea_plot+geom_boxplot(data=random_G4[random_G4$key<70,],aes(x=key,y=value, group=cut_width(key,2)), color='red', size = 0.2, outlier.shape = 19, outlier.size = 0.5)
ggsave('~/Documents/archaea/archea_control.png', plot = archea_control, width = 10, height = 8,
       units = 'in')
archea_thermococci = archea%>%filter(SubGroup=='Thermococci')
archea_thermococci['Name']=archea_thermococci$`Organism/Name`
archea_thermococci=archea_thermococci%>%separate('Name', into=c('genus','species','n1','n2'), sep=' ', remove=FALSE)
themococci=ggplot(archea_thermococci)+ geom_point(shape=21, color='#000001',
                                       aes(x=GC,y=G4_density,
                                           size=spare_percentage,
                                           fill=genus))+
  scale_fill_hue()+
  facet_wrap(~SubGroup, nrow = 3)+
  ggtitle('G4 density over GC percentage for Archea')+
  theme(axis.text = element_text(size = 18),axis.title = element_text(size = 20,face='bold'),plot.title = element_text(size = 20, face = 'bold'))
archea['Name']=archea$`Organism/Name`
archea=archea%>%separate('Name', c('genus','species', 'n1', 'n2'), sep=' ', remove=FALSE)
ggplot(archea%>%filter(SubGroup=='Methanomicrobia'))+ geom_point(shape=21, color='#000001',
                                       aes(x=GC,y=G4_density,
                                           size=spare_percentage,
                                           fill=genus))+
  scale_fill_hue()+
  ggtitle('G4 density over GC percentage for Archea')+
  theme(axis.text = element_text(size = 18),axis.title = element_text(size = 20,face='bold'),plot.title = element_text(size = 20, face = 'bold'))
##------------------thermococci--------
thermococci_meta = archea_meta[archea_meta$SubGroup=='Thermococci', ]

ggplot(melt_archaea[melt_archaea$SubGroup=='Thermococci',], 
       aes(x=as.factor(distance_to_tss), y = value)) + geom_boxplot(aes(fill = genus)) + 
  my_theme + 
  facet_wrap(~SubGroup, ncol = 2) + 
  labs(x= 'Distance to TSS', y = 'G4 fraction around TSS')
ggsave('archea_thermococci_PQS_distribution.png', path = save_path,width = 10, height = 8,
       units = 'in')
ggplot(archea[archea$SubGroup=='Thermococci', ])+ geom_point(shape=21, color='#000001',
                           aes(x=GC,y=G4_density, fill=genus), size = 5, alpha = 0.7)+
  scale_fill_hue()+
  xlab('GC%') + ylab('PQS density') + geom_boxplot(data=random_G4[random_G4$key<70,],aes(x=key,y=value, group=cut_width(key,2)), color='red', size = 0.2, outlier.shape = 19, outlier.size = 0.5) + 
  ggtitle('PQS density over GC percentage for Thermococci')+my_theme
ggsave('archea_thermococci.png', path = save_path,width = 10, height = 8,
       units = 'in')
