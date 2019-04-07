library(grid)
library(gridExtra)
library(tidyverse)

## TABLE 1: confusion table definitions
classifications <- c("True Positives", "False Positives", "True Negatives", "False Negatives", "MCC")
values <- c(3929, 2619, 3221774, 29422, 0.234)
dat <- cbind(classifications, values)
colnames(dat) <- c("Classifications", "Scores")

t1 <- ttheme_default(core=list(
  fg_params=list(fontface=c(rep("plain", 4), "bold"), cex = 0.75),
                   alpha = rep(c(1,0.5), each=5),
  colhead = list(fg_params=list(cex = 0.75)),
  rowhead = list(fg_params=list(cex = 0.75)))
)

g <- tableGrob(dat, theme = t1)
h <- convertHeight(sum(g$heights), "in", TRUE)
w <- convertWidth(sum(g$widths), "in", TRUE)
ggsave(g, filename = "MCC.png", width=w, height=h, dpi=400)

## TABLE 2: metadata table for each dataset used in analyses
metadata_table <- read_tsv("~/Documents/UBC/YEAR4/MICB449/datasets/raw_metadata.txt")

t2 <- metadata_table %>%
  select(run_accession, read_count, scientific_name, class, recA, rpoB, rps8, pyrG, ribosomal_L10P) %>%
  rename("RPL10" = ribosomal_L10P) %>%
  rename("Dataset ID" = run_accession) %>%
  rename("Read Count" = read_count) %>%
  rename("Taxonomy" = scientific_name) %>%
  rename("Class" = class)

g2 <- tableGrob(t2)
h <- convertHeight(sum(g2$heights), "in", TRUE)
w <- convertWidth(sum(g2$widths), "in", TRUE)
ggsave(g2, filename = "dataset_metadata.png", width=w, height=h)

## TABLE 3: comparison of sanger vs illumina vs nanopore reads
sanger <- c(">800 bp", "Template using pool of molecules", "$0.00125 / bp", "Highest quality, relatively low throughput")
illumina <- c("2x150 bp", "Template clonally derived from single molecule", "$6.67 x 10^-7 / bp", "High quality, extremely high throughput")
nanopore <- c(">6-8 kbp", "Template is single nucleic acid sequence", "$4.75 x 10^-9 / bp", "Poor quality, comparitively low throughput")

t3 <- tibble(sanger, illumina, nanopore)
colnames(t3) <- c("Sanger reads", "Illumina reads", "Oxford Nanopore reads")

g3 <- tableGrob(t3)
h <- convertHeight(sum(g3$heights), "in", TRUE)
w <- convertWidth(sum(g3$widths), "in", TRUE)
ggsave(g3, filename = "comparison_sequence_tech.png", width=w, height=h, dpi=400)

