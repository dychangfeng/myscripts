library(tidyr)
library(stringr)
library(dplyr)
library(ggplot2)
my_theme <-   theme(axis.text = element_text(size = 16),
                    axis.title = element_text(size = 18,face = "bold"),
                    plot.title = element_text(size = 16, face = 'bold'),
                    strip.text.x = element_text(size=14, face="bold"),
                    strip.background = element_rect(colour="burlywood", fill="burlywood"),
                    legend.text = element_text(size = 16, face = 'bold' ),
                    legend.title = element_text(size=16, face='bold'))
G4_spare =read.csv("~/Documents/bacteria_G4/G4_count_random_2m_20_80_spare_percent.csv")
G4_spare = G4_spare[,2:7]
colnames(G4_spare) <- seq(20,70,10)
random_G4 <- G4_spare%>%gather() ## gather all the column names as key , values in one columns
random_G4$key <- as.numeric(random_G4$key)
ggplot(random_G4)+geom_boxplot(aes(x=as.factor(key), y=value)) + my_theme + xlab('GC%') + ylab('>4 track %')
ggsave('~/Documents/bacteria_G4/figures/random_spare_fraction.png', scale = 1, width = 7,
       height = 7, units='in')
