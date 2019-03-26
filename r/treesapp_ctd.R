library(tidyverse)
library(RColorBrewer)

raw_dat <- read_tsv("data/ctd_scores.txt")
head(raw_dat)

dat <- mutate_if(raw_dat, is.character, str_replace_all, pattern = "rL10P", replacement = "RPL10")

fig <- ggplot(dat, aes(x = Marker, y = CTD_score)) +
  geom_jitter(aes(fill = Dataset), shape = 21, size = 3, alpha = 2/3, position = position_jitter(0.15)) +
  scale_fill_brewer(palette = "PuOr") +
  xlab("Marker Gene") +
  ylab("Cumulative Taxonomic Distance") +
  ggtitle("CTD Per Marker Gene") +
  theme(panel.background = element_blank(),
        panel.grid.major = element_line(colour = "#bdbdbd", linetype = "dotted"),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5))

ggsave(fig, filename = "treesapp_ctd.png", width = 8, height = 5, dpi = 400)
 