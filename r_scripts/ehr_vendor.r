library(readr)
library(tidyr)
library(dplyr)
library(ggplot2)
ehr <- read_csv("~/Desktop/EHR-vendors-count-dataset.csv")
ehr_total <- ehr%>%select('developer','provider_type', 'tot_provs_report_developer')
ehr_agg <- aggregate(. ~ developer+provider_type, data = ehr_total, sum)
ehr_agg <- ehr_agg[with(ehr_agg, order(-tot_provs_report_developer)),] ## sort by total number
ehr_agg <- ehr_agg%>%mutate(percent=100*tot_provs_report_developer/sum(tot_provs_report_developer))%>%
  mutate(agg_percent=cumsum(percent))
write_csv(ehr_agg, '~/Desktop/ehr_agg.csv' )
