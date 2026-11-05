# plot_modification.R
  height = 665 / 300,
library(dplyr)
library(stringr)
library(IRanges)
library(ggplot2)
library(scales)
library(colorspace)
library(broom)
library(forcats)
library(data.table)

valid_chroms <- c(as.character(1:22), "X", "Y")
chrom_dfs <- list()

for (chrom in valid_chroms) {
  awk_cmd <- sprintf("awk '$1 == \"%s\"' basic_mod_sample2.bedMethyl | head -n 100000", chrom)
  dt <- fread(cmd = awk_cmd, sep = "\t", header = FALSE)
  if (nrow(dt) > 0) {
    chrom_dfs[[chrom]] <- dt
  }
}

bedMethyl <- rbindlist(chrom_dfs, use.names = FALSE, fill = FALSE)

colnames(bedMethyl) <- c(
  "chrom", "start_position", "end_position", "modified_base", "score",
  "strand", "start_position2", "end_position2", "color", "Nvalid_cov",
  "fraction_modified", "Nmod", "Ncanonical", "Nother_mod", "Ndelete",
  "Nfail", "Ndiff", "Nnocall"
)

#100kb windows of chromq
set.seed(123)
sampled <- bedMethyl[, .SD[sample(.N, min(.N, 100000))], by = chrom]
filtered_sampled <- sampled[
  !is.na(fraction_modified) & fraction_modified < 100,
]
#filter methylation status
filtered_sampled$methylation_status <- ifelse(
  filtered_sampled$fraction_modified == 0, "Unmethylated",
  ifelse(filtered_sampled$modified_base == "h", "5hmC", "5mC")
)

#plot
p <- ggplot(filtered_sampled, aes(x = chrom, y = fraction_modified)) +
  geom_jitter(aes(color = methylation_status), width = 0.2, alpha = 0.7, size = 0.8) +
  scale_color_manual(
    values = c("5hmC" = "deeppink2", "5mC" = "cadetblue3", "Unmethylated" = "grey"),
    name = "Methylation Status"
  ) +
  labs(
    x = "Chromosome",
    y = "Probability of 5hmC/5mC Modification (%)",
    title = "Genomic Distribution of Modified Cytosine of Normal Human Astrocytes"
  ) +
  stat_summary(
    fun = mean,
    geom = "crossbar",
    width = 0.5,
    color = "blue",
    position = position_dodge(width = 0.9)
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 5, angle = 0, hjust = 1, color = "black"),
    axis.text.y = element_text(size = 5, color = "black"),
    axis.title.x = element_text(size = 5, color = "black"),
    axis.title.y = element_text(size = 5, color = "black"),
    legend.text = element_text(size = 5, color = "black"),
    legend.title = element_text(size = 5, color = "black"),
    legend.position = "top",
    legend.justification = "center",
    legend.direction = "horizontal",
    plot.title = element_text(hjust = 0.5, color = "black", size = 8)
  )
#save
ggsave(
  filename = "mod_fraction_horizontal_plot.png",
  plot = p,
  width = 2338 / 300,
  height = 665 / 300,
  dpi = 300,
  units = "in")
