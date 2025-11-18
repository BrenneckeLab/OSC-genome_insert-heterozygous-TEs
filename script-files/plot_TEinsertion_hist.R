library(tidyverse)
library(cowplot)
library(plotly)
theme_set(theme_cowplot())


###################################################################################################
args  =  commandArgs(TRUE);
argmat  =  sapply(strsplit(args, "="), identity)

for (i in seq.int(length=ncol(argmat))) {
  assign(argmat[1, i], argmat[2, i])
}

# available variables
print(ls())

###################################################################################################

TABLE=read_tsv(INPUT , col_names = TRUE)

TABLE = TABLE %>%
  mutate(TEfraction = TElength*100/InsertionLength)

p=ggplot(TABLE, aes(x=InsertionLength, y=TEfraction))+
  geom_point(alpha=0.3, size=1, shape=18)+
  scale_x_log10()+
  labs(y="fraction TE",
       title="insertion length vs TE fraction")
NAME=paste0(OPENdir,"detectedInsertions.length_vs_TE.",EXT,".png")
ggsave(NAME,p)

p=ggplot(TABLE, aes(x=InsertionLength))+
  geom_histogram(bins=50)+
  scale_y_log10()+
  labs(x="insertion length",
     y="n insertions in bin",
     title="binned insertion length")
NAME=paste0(OPENdir,"detectedInsertions.insertionLength_binned.",EXT,".png")
ggsave(NAME,p)

p=ggplot(TABLE, aes(x=TEfraction))+
  geom_histogram(binwidth=5)+
  labs(x="fraction TE per insertion",
       y="n insertions in bin",
       title="binned fraction TE per insertion")
NAME=paste0(OPENdir,"detectedInsertions.TEfraction_binned.",EXT,".png")
ggsave(NAME,p)
