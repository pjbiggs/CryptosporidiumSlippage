### Plotting the data

## Getting set up

# Load the required packages:

library(dada2)
packageVersion("dada2")
#library(ShortRead)
#packageVersion("ShortRead")
#library(ggplot2)
#packageVersion("ggplot2")
library(tidyverse)
packageVersion("tidyverse")
#library(plyr)
#packageVersion("plyr")
#library("viridis")
#packageVersion("viridis")


## load data back in ##

sampleList <- read_tsv("../Samples.txt")
species <- read_tsv("species.txt") %>% as_tibble()

seqData <- read_tsv("dada_seqsMod.txt", col_names = FALSE)
colnames(seqData) <- c("seqID", "sequence")

# join in the sequences
combinedA <- full_join(species, seqData)

# add seqID\t for the first sequence 
tableCounts <- read_tsv("dada_table.txt") %>% as_tibble()
head(tableCounts)

# proportionate the counts on a per column basis
df <- tableCounts[2:37]
propCounts <- as.data.frame(apply(df, 2, function(i) i/sum(i)))
propCounts$seqID <- tableCounts$seqID
rm(df)

combinedAll <- full_join(combinedA, propCounts)

# as per the manuscript, take the 20 most prevalent sequences only
combinedAll <- combinedAll[1:20, ]

# by sample
rangesFullTrim <- select(combinedAll, -c(seqID, sequence))
write_tsv(rangesFullTrim, "proportionTop20.txt")

# data is then curated to sum by sequence length for each sample #
propTotal <- select(rangesFullTrim, -c(clusterID, indicative)) %>% group_by(length, Species) %>% 
  summarise(across(everything(), sum)) 
head(propTotal)

rangesFullTrimLong <- gather(propTotal, "AltName", "sampleVal", 3:38)
rangesFullTrimLongS <- inner_join(rangesFullTrimLong, sampleList)

ggplot(rangesFullTrimLongS, aes(fill=Species, y=sampleVal, x=length)) + 
  geom_bar(position="stack", stat="identity") +
  facet_grid(HPMix ~ Group) + theme(legend.position="bottom") + 
  scale_fill_manual(values=c("#A1B654", "#2C3B75")) +
  ggtitle("Barplot of the control combinations grouped by\nproportion of C. hominis:C. parvum for BANDSIZE=2") + 
  xlab("amplicon length (bp)") + ylab("proportion of reads across 3 replicates")
ggsave("controlsInOrderAll.jpg", width = 16, height = 20, units = "cm")
