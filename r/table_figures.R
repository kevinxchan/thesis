library(grid)
library(gridExtra)
library(tidyverse)

## TABLE 1: confusion table definitions
classifications <- c("True Positives", "False Positives", "True Negatives", "False Negatives", "MCC")
values <- c(3929, 2619, 3221774, 29422, 0.234)
dat <- cbind(classifications, values)
colnames(dat) <- c("Classifications", "Scores")

t1 <- ttheme_default(core=list(
  fg_params=list(fontface=c(rep("plain", 4), "bold")),
                   alpha = rep(c(1,0.5), each=5))
)

g <- tableGrob(dat, theme = t1)
png("100kHighRes300dpi.png", units="px", width=1600, height=1600, res=300)
plot(g, main="100,000 points", col=adjustcolor("black", alpha=0.2))
dev.off()
g

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
