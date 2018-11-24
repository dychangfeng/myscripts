library(dplyr)
library(readr)
## readr package won't convert character to factor
taxa.info = read_delim('~/Desktop/archive-kingdom-bacteria-bl3/taxa.txt',delim = '\t')
bacterial_final <- read_delim('~/Documents/bacteria_G4/bacteria_representive/G4_in_bacteria_rep_final.txt', delim ='\t')

mean(unique(bacterial_final$genus)%in%taxa_order$genus) ## 84 percent can map back to have an order
taxa_order = taxa.info%>%select(class, order, family, genus)
taxa_order = drop_na(taxa_order)%>%distinct(.keep_all = TRUE) # drop duplicated, and NA
bacterial_all = left_join(bacterial_final, taxa_order, by='genus')
write_delim(bacterial_all, '~/Documents/bacteria_G4/bacteria_with_order.txt', delim = '\t', col_names = TRUE)
##betaproteo_all = betaproteo_all[-which(duplicated(betaproteo_all)),]
myplot2 <- function(mydata){ggplot(data = mydata) +
    geom_point(shape=21, colour="#000000",aes(x=G4_density, y=fold_enrichment, fill = G4_in_CDS),size = 2) +
    scale_fill_continuous(low='blue', high ='red')+
    facet_wrap(~order, nrow = 3)
    ggtitle('G4 density over GC percentage for Proteobacteria') +my_theme}
ggplot(betaproteo_all) + geom_point(shape=21, colour="#000000",aes(x=G4_density, y=fold_enrichment, fill = order),size = 5) + my_theme
##--------------------add to cyano group---------------------------
cyano = read_delim('~/Documents/bacteria_G4/cyno_bacteria/cyano_all_final.txt', delim = '\t')
#cyano = cyano%>%rename(genus = Genus) ## all with lower first letter
library(tools)
#cyano['genus']=lapply(cyano$genus, toTitleCase)%>%unlist() ## convert to title case
taxa_order%>%filter(class == 'Cyanophyceae') ## choose only cyano bacteria
cyano_order = left_join(cyano, taxa_order, by='genus')
write_delim(cyano_order, '~/Documents/bacteria_G4/cyno_bacteria/cyano_with_order.txt', delim = '\t')
