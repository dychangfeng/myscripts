library(tidyr)
library(stringr)
library(dplyr)
library(ggplot2)
bacteria_1464_with_info <- read.csv("~/Documents/bacteria_G4/bacteria_representive/bacteria_1464_with_info.csv")
bacteria_1464_with_info$SubGroup <- as.factor(bacteria_1464_with_info$SubGroup)
random_G4 <- read.csv("~/Documents/bacteria_G4/G4_count_random_2m_20_80_2.csv")
colnames(random_G4) <- seq(18,78,2)
random_G4$`18` <- NULL ## remove first column
random_G4 <- random_G4%>%gather() ## gather all the column names as key , values in one columns
random_G4$key <- as.numeric(random_G4$key)

##-----------------------------------process file---------

bacteria_1464_with_info <-  bacteria_1464_with_info%>%mutate(proteins_density = Proteins/genome_mbp)%>%
  separate('X.Organism.Name', into=c('genus','n0','n1','n2'),sep=" ", remove=FALSE)%>%select(-n0, -n1, -n2)
my_theme <-   theme(axis.text = element_text(size = 16),
                    axis.title = element_text(size = 18,face = "bold"),
                    plot.title = element_text(size = 16, face = 'bold'),
                    strip.text.x = element_text(size=14, face="bold"),
                    strip.background = element_rect(colour="burlywood", fill="burlywood"),
                    legend.text = element_text(size = 16, face = 'bold' ),
                    legend.title = element_text(size=16, face='bold'))
##----check out each group---------
group_no_care = c('Caldiserica','Calditrichaeota', 'Chrysiogenetes','Dictyoglomi','Elusimicrobia'
                  ,'Thermodesulfobacteria')

bb <- ggplot(bacteria_1464_with_info[!bacteria_1464_with_info$Group %in% group_no_care,]) + geom_point(aes(x = GC_percentage, 
                                                       y = GC_density),shape=21, colour="#000000", fill = 'grey') +my_theme+
  ggtitle('')+
  labs(x = 'GC %', y = 'PQS Density (N/Mbp)')+
  facet_wrap(~Group, nrow = 4)

bacteria_control=bb + geom_boxplot(data = random_G4, aes(x=key,y=value, group=cut_width(key,2)), color='red', alpha=0.2)
ggsave('G4_bacteria_with_control2.svg', bacteria_control, path = '~/Documents/bacteria_G4/bacteria_representive/figures/', scale = 1, width = 7,
       height = 7, units='in')
terra_panel <-  ggplot(Terrabacteria) + geom_point(shape = 21, colour = '#000000',aes(x = GC_percentage, y = GC_density, 
                                                                                      size = spare_percentage,fill = proteins_density)) +
  my_theme+facet_wrap(~SubGroup, nrow = 3) + scale_fill_continuous(low = 'blue', high = 'red')
##----------------Terrabacteria----------------------
Terrabacteria <- bacteria_1464_with_info%>%filter(Group=='Terrabacteria group')
library(ggplot2)
t1 <- ggplot(data = Terrabacteria) +
  geom_point(shape=21, aes(x=GC_percentage, y=GC_density, fill = SubGroup),colour="#000000", size =3) +
  facet_wrap(~SubGroup, ncol = 2)+
  ggtitle('G4 density over GC percentage for Terrabacteria') +
  labs(x = 'GC Percentage', y = 'G4 Density') +
  guides(fill =guide_legend(title.theme = element_text(
    size = 18,
    face = "italic",
    colour = "red",
    angle = 0
  )
  )
  ) +my_theme
ggsave('G4_terrabacteria_protein_density.png', t1, path = '~/Desktop/', scale = 1, width = 10,
       height = 8, units='in')

t1_control= t1+geom_boxplot(data= random_G4, aes(x=key, y=value, group=cut_width(key,2)), color='black')

t2 <- t1 + scale_y_log10() 
#-----------filter out the subgroup with less genomes-----
terra_filter <- Terrabacteria%>%filter(SubGroup=='Actinobacteria'|SubGroup=='Deinococcus-Thermus'|SubGroup=='Firmicutes'|SubGroup=='Cyanobacteria/Melainabacteria group')
tt1 <- ggplot(data = terra_filter) +
  geom_point(shape=21, colour="#000000",aes(x=GC_percentage, y=GC_density, fill = SubGroup), size =5) +
  facet_wrap(~SubGroup, ncol = 2)+
  ggtitle('G4 density over GC percentage for Terrabacteria with filter') +
  labs(x = 'GC Percentage', y = 'G4 Density')+ my_theme
ggsave('G4_terrabacteria_control.png', terra_control, path = '~/Desktop/', scale = 1, width = 10,
       height = 8, units='in')
terra_control=tt1+geom_boxplot(data= random_G4, aes(x=key, y=value, group=cut_width(key,2)), color='black')


##-----Proteobacteria-------------------
Proteobacteria <- bacteria_1464_with_info%>%filter(Group=='Proteobacteria')
p1 <- ggplot(data = Proteobacteria) +
  geom_point(shape=21, colour="#000000",aes(x=GC_percentage, y=GC_density, size = spare_percentage, fill = proteins_density)) +
  scale_fill_continuous(low='blue', high ='red')+
  facet_wrap(~SubGroup, nrow = 2)+
  ggtitle('G4 density over GC percentage for Proteobacteria') +
  my_theme+
  labs(x = 'GC Percentage', y = 'G4 Density')
ggsave('G4_proteobacteria_control.png', plot = proteo_control, path = '~/Desktop/', scale = 1, width = 10,
       height = 8, units='in')
proteo_control=p1+geom_boxplot(data=random_G4,aes(x=key, y=value,group=cut_width(key,2)), color='black')

p2 <- p1 + scale_y_log10()
##------------FCB group----------------
FCB_group <- bacteria_1464_with_info%>%filter(Group=='FCB group')
f1 <- ggplot(data = FCB_group) +
  geom_point(aes(x=GC_percentage, y=GC_density, size = spare_percentage, fill = SubGroup),shape=21, colour="#000000") +
  scale_fill_hue()+
  facet_wrap(~SubGroup, nrow = 2)+
  ggtitle('G4 density over GC percentage for FCB group') +
  my_theme+
  labs(x = 'GC Percentage', y = 'G4 Density')
ggsave('G4_FCB_group.png', plot = f1, path = '~/Desktop/', scale = 1, width = 10,
       height = 8, units='in')
f1_control = f1 + geom_boxplot(data = random_G4, aes(x=key, y=value,group=cut_width(key,2)), color ='black')
##-------------check out each Group-----------
pro_panel <- ggplot(Proteobacteria) + geom_point(shape = 21, colour = '#000000',aes(x = GC_percentage, y = GC_density, 
                                                                                    size = spare_percentage,fill = proteins_density)) +
  facet_wrap(~SubGroup, nrow = 4) + scale_fill_continuous(low = 'blue', high = 'red')
##----------Actinobacteria-----
Actinobacteria <- bacteria_1464_with_info%>%filter(SubGroup=='Actinobacteria')
Actinobacteria <- Actinobacteria%>%separate('X.Organism.Name', into=c('genus','species','n1','n2'), 
                                            sep=" ", remove=FALSE)
a1 <- ggplot(data = Actinobacteria, 
             aes(x=GC_percentage, y=GC_density, size = spare_percentage, fill=genus)) +
  geom_point(shape=21, colour="#000000") +
  scale_fill_discrete()+
  ggtitle('G4 density over GC percentage for Actinobacteria') +
  labs(x = 'GC Percentage', y = 'G4 Density')+my_theme

ggsave('G4_actinobacteria_genus.png', plot = a1, path = '~/Desktop/', scale = 1, width = 10,
       height = 8, units='in')
a2 <- a1 + scale_y_log10()
actino_750 <- Actinobacteria%>%filter(GC_density>=750)
actino_750$X.Organism.Name
##-------------filter genus of actinobacteria--------
actino_genus = ggplot(data = Actinobacteria[duplicated(Actinobacteria$genus),]) +
  geom_point(aes(x=GC_percentage, y=GC_density,  size = spare_percentage, fill=genus),shape=21, colour="#000000") +
  scale_fill_discrete()+
  ggtitle('G4 density over GC percentage for Actinobacteria') +
  labs(x = 'GC Percentage', y = 'G4 Density')+ my_theme

actino_with_control=actino_genus+geom_boxplot(data = random_G4, aes(x=key, y=value,group=cut_width(key,2)), color ='black')
ggsave('G4_actinobacteria_genus_with_control.png', plot = actino_with_control, path = '~/Desktop/', scale = 1, width = 10,
       height = 8, units='in')
##----------------Deinococcus_thermus-------------
Deinococcus_thermus<- bacteria_1464_with_info%>%filter(SubGroup=='Deinococcus-Thermus')
d1 <- ggplot(data = Deinococcus_thermus, 
             aes(x=GC_percentage, y=GC_density, size = spare_percentage, fill=proteins_density)) +
  geom_point(shape=21, colour="#000000") +
  scale_fill_continuous(low='blue', high='red')+
  ggtitle('G4 density over GC percentage for Deinococcus-Thermus') +
  labs(x = 'GC Percentage', y = 'G4 Density')+
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20,face = "bold"),
        plot.title = element_text(size = 20, face = 'bold'),
        strip.text = element_text(size = 16))
ggsave('G4_D_thermus.png', plot = d1, path = '~/Desktop/', scale = 1, width = 10,
       height = 8, units='in')
d2 <- a1 + scale_y_log10()
##---------check out each group------------
Actino_panel <- ggplot(Actinobacteria) + geom_point(aes(x = GC_percentage, y = GC_density)) +
  facet_wrap(~SubGroup, nrow = 6)
##---------Firmicutes-------------------
bacteria_1464_with_info['G4_density'] <- bacteria_1464_with_info$GC_density
Firmicutes <- bacteria_1464_with_info%>%filter(SubGroup=='Firmicutes')
which(Firmicutes$G4_density==max(Firmicutes$G4_density))
f1 <- ggplot(data = Firmicutes) +
  geom_point(aes(x=GC_percentage, y=GC_density, size = spare_percentage, fill=proteins_density),shape=21, colour="#000000") +
  scale_fill_continuous(low='blue', high='red')+
  ggtitle('G4 density over GC percentage for Firmicutes') +
  labs(x = 'GC Percentage', y = 'G4 Density')+
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20,face = "bold"),
        plot.title = element_text(size = 16, face = 'bold'),
        strip.text = element_text(size = 16))
ggsave('G4_Firmicutes_with_control.png', plot = Frimicutes_with_control, path = '~/Desktop/', scale = 1, width = 10,
       height = 8, units='in')
Frimicutes_with_control=f1 + geom_boxplot(data = random_G4, aes(x=key, y=value,group=cut_width(key,2)), color ='black') 
f2 <- f1 + scale_y_log10()
##----------Alphaproteobacteria---
Alphaproteobacteria <- bacteria_1464_with_info%>%filter(SubGroup=='Alphaproteobacteria')
al1 <- ggplot(data = Alphaproteobacteria, 
              aes(x=GC_percentage, y=GC_density, size = spare_percentage, fill=Proteins)) +
  geom_point(shape=21, colour="#000000") +
  scale_fill_continuous(low='blue', high='red')+
  ggtitle('G4 density over GC percentage for Alphaproteobacteria') +
  labs(x = 'GC Percentage', y = 'G4 Density')
al2<- p1 + scale_y_log10()

