library(tidyverse)
hopped_summary <- as.data.frame(read.table("/Users/ryan/Courses/Bi624/demultiplexing-rylanci/script_data/map.out"))

demulti_summary <- as.data.frame(read.table("/Users/ryan/Courses/Bi624/demultiplexing-rylanci/demulti_output/demulti_summary.out"))

demulti_summary <- demulti_summary$V1 /4 

plot(hopped_summary$V1, hopped_summary$V2 /4)

plot(hopped_summary$V3, hopped_summary$V4 /4)


ggplot(hopped_summary) + 
  geom_point(aes(V1, V2/4)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  xlab("Hopped Index") + ylab("Occurence") +
  ggtitle("Occurence of Index Hopping in R1 Barcodes")


ggplot(hopped_summary) + 
  geom_point(aes(V3, V4/4)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  xlab("Hopped Index") + ylab("Occurence") +
  ggtitle("Occurence of Index Hopping in R2 Barcodes")



ggplot(demulti_summary) + 
  geom_point(aes(V2, V1/4)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  scale_y_log10() + xlab("Barcode Bin") + ylab("log10Record Count") +
  ggtitle("Records per Barcode Bin")

