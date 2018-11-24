library(tidyr)
library(stringr)
library(dplyr)
library(ggplot2)
##-----------plotting setup----------------
my_theme <-   theme(legend.direction = 'horizontal', legend.position = 'top',
                    axis.text = element_text(size = 18),
                    axis.title = element_text(size = 20,face = "bold"),
                    plot.title = element_text(size = 20, face = 'bold'),
                    strip.text.x = element_text(size=20, face="bold"),
                    strip.background = element_rect(colour="burlywood", fill="burlywood"),
                    legend.text = element_text(size = 16, face = 'bold' ),
                    legend.title = element_text(size=16, face='bold'))
myplot2 <- function(mydata, color_by){ggplot(data = mydata) +
    geom_point(shape=21, colour="#000000",aes(x=G4_density, y=fold_enrichment, fill = color_by),size = 5) + 
     my_theme}
mybar_family <- function(mydata){
  ggplot(mydata) + geom_bar(aes(family, fold_enrichment, fill=family ),
        position = "dodge", stat = "summary", fun.y = "mean")+theme(axis.text.x = element_text(angle = 60, hjust = 1))+my_theme
}
mybar_order <- function(mydata){
  ggplot(mydata) + geom_bar(aes(order, fold_enrichment, fill=order ),
                            position = "dodge", stat = "summary", fun.y = "mean")+
    theme(axis.text.x = element_text(angle = 60, hjust = 1))+my_theme
}
save_path = '/Users/Yun/Documents/bacteria_G4/D_thermus/stats_fig/'
##--------------------------------get data--------------------------------------------
bacterial_final <- read.delim('~/Documents/bacteria_G4/bacteria_with_order.txt', header = T)
bacterial_final <- bacterial_final%>%filter(G4_density>10)
random_G4 <- read.csv("~/Documents/bacteria_G4/G4_count_random_2m_20_80_2.csv")
colnames(random_G4) <- seq(18,78,2)
random_G4$`18` <- NULL ## remove first column
random_G4 <- random_G4%>%gather() ## gather all the column names as key , values in one columns
random_G4$key <- as.numeric(random_G4$key)
random_p = geom_boxplot(data= random_G4, aes(x=key, y=value, group=cut_width(key,2)), color='red', alpha = 0.2)
mygroup <- c('Acidobacteria', 'Aquificae', 'FCB group', 'Nitrospirae', 'Proteobacteria', 
             'PVC group', 'Spirochaetes',  'Terrabacteria group', 'Thermotogae')
ggplot(bacterial_final%>%filter(Group %in% mygroup)) +
  geom_point(shape=21, colour="#000000",aes(x=G4_density, y=fold_enrichment, fill=Group),size = 2) +
  facet_wrap(~Group, nrow = 3)+
  my_theme + xlab('G4_density') + ylab('Fold Enrichment')
##-----------CDS vs genome length----------------
ggplot(bacterial_final) + geom_point(aes(x=genome_mbp,y=cds_Mb/genome_mbp), size = 2,
                                     fill='cyan4')+my_theme+
  xlab('Genome length (mbp)')+
  ylab('CDS fraction')
##---------------proteobacteira----------------
pro_subgroup = bacterial_final%>%filter(Group=='Proteobacteria')
myplot <- function(mydata){ggplot(data = mydata) +
  geom_point(shape=21, colour="#000000",aes(x=G4_density, y=fold_enrichment, fill = order),size = 4) +
  facet_wrap(~SubGroup, nrow = 3)+my_theme}
acetobacteraceae = pro_subgroup%>%filter(family=='Acetobacteraceae')
myplot(pro_subgroup)
myplot2(acetobacteraceae, color_by = acetobacteraceae$genus)
myplot2(pro_subgroup%>%filter(SubGroup=='Alphaproteobacteria'))
myplot2(pro_subgroup%>%filter(SubGroup=='Betaproteobacteria'))
myplot2(pro_subgroup%>%filter(SubGroup=='Gammaproteobacteria'))
myplot2(pro_subgroup%>%filter(SubGroup=='delta/epsilon subdivisions'))

mybar_family(pro_subgroup%>%filter(SubGroup=='Alphaproteobacteria'))
mybar_family(pro_subgroup%>%filter(SubGroup=='Betaproteobacteria'))
mybar_order(pro_subgroup%>%filter(SubGroup=='Gammaproteobacteria'))
mybar_order(pro_subgroup%>%filter(SubGroup=='delta/epsilon subdivisions'))
## ------------------Terra--------------
Terra_subgroup = bacterial_final%>%filter(Group=='Terrabacteria group')
myplot(Terra_subgroup)
actino =Terra_subgroup%>%filter(SubGroup=='Actinobacteria')
myplot2(Terra_subgroup%>%filter(SubGroup=='Cyanobacteria/Melainabacteria group')) 
myplot2(actino, actino$family)
myplot2(Terra_subgroup%>%filter(SubGroup=='Chloroflexi'))
myplot2(Terra_subgroup%>%filter(SubGroup=='Firmicutes', G4_counts_at_tss>10))
myplot2(Terra_subgroup%>%filter(SubGroup=='Firmicutes', G4_counts_at_tss>10))

ggplot(actino%>%filter(order=='Actinomycetales')) + geom_bar(aes(family, fold_enrichment, fill=family ),
                  position = "dodge", stat = "summary", fun.y = "mean")
mybar(actino%>%filter(order=='Actinomycetales'))
mybar_order(actino)
##-------------cyanobacteria-------------------------------
cyano = read.delim('~/Documents/bacteria_G4/cyno_bacteria/cyano_with_order.txt', sep = '\t')
cyano = filter(cyano, G4_counts_at_tss>10)
myplot2(cyano, color_by = cyano$family)
ggplot(cyano) +
  geom_point(shape=21, colour="#000000",aes(x=G4_density, y=fold_enrichment),size = 3, fill = 'brown') +
  my_theme
mybar_family(cyano)
ggplot(cyano) +
  geom_point(shape=21, colour="#000000",aes(x=GC_percent, y=G4_density, fill = order),size = 5) +
  my_theme + random_p + xlab('GC %') + ylab('PQSs Density (N/Mbp)' )
ggsave('pqss density of cyanobacteria.png', path = save_path, scale = 1, width = 10,
       height = 8, units='in')
##--------------------D_Thermus------
D.thermus <- read.delim('~/Documents/bacteria_G4/D_thermus/G4_in_D_thermus_final_summary.txt', header = T)
mybar_order(D.thermus)
ggplot(D.thermus) +
  geom_point(shape=21, colour="#000000",aes(x=G4_density, y=fold_enrichment, fill = genus),size = 5) +
my_theme

ggplot(D.thermus) +
  geom_point(shape=21, colour="#000000",aes(x=GC_percentage*100, y=G4_density, fill = order),size = 5) +
  my_theme + random_p + xlab('GC %') + ylab('PQSs Density (N/Mbp)' )
ggsave('pqss density of dt.svg', path = save_path, scale = 1, width = 5,
       height = 4, units='in')
ggplot(D.thermus) +
  geom_point(shape=21, colour="#000000",aes(x=G4_in_CDS, y=cds_Mb/Genome_length, fill = order),size = 5) +
  my_theme + xlim(0,1) + ylim(0,1)
##------------------------FCB -----------
FCB_group <- bacterial_final%>%filter(Group=='FCB group')
sub_fcb = FCB_group%>%filter(SubGroup=='Bacteroidetes/Chlorobi group')
myplot(FCB_group%>%filter(G4_density>10))
myplot2(sub_fcb, fill_by = sub_fcb$family)





