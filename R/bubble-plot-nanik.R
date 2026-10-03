
#########################Thorsten_Nanik###
# Install packages if needed
if (!requireNamespace("readxl", quietly = TRUE)) install.packages("readxl")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("grid", quietly = TRUE)) install.packages("grid")

# Load libraries
library(readxl)
library(ggplot2)
library(dplyr)
library(grid)

# Load dataset
data <- read_excel("xxxx.xlsx")

# Check column names
str(data)

# Convert GeneRatio from "61/429" to numeric ratio
data <- data %>%
  mutate(
    GeneRatio_numeric = sapply(GeneRatio, function(x) {
      parts <- strsplit(as.character(x), "/")[[1]]
      as.numeric(parts[1]) / as.numeric(parts[2])
    }),
    p.adjust = as.numeric(p.adjust),
    Count = as.numeric(Count),
    logP = -log10(p.adjust)
  ) %>%
  filter(
    !is.na(Description),
    !is.na(GeneRatio_numeric),
    !is.na(p.adjust),
    !is.na(Count),
    p.adjust > 0
  )

# Order pathways by significance
data$Description <- factor(
  data$Description,
  levels = data$Description[order(data$logP)]
)

# Create bubble chart
bubble_chart <- ggplot(
  data,
  aes(
    x = GeneRatio_numeric,
    y = Description
  )
) +
  geom_point(
    aes(
      size = Count,
      color = logP
    ),
    alpha = 0.9
  ) +
  scale_size_continuous(
    range = c(4, 12)
  ) +
  scale_color_gradient(
    low = "darkblue",
    high = "#FF5722"
  ) +
  scale_x_continuous(
    limits = c(0, max(data$GeneRatio_numeric, na.rm = TRUE) * 1.15),
    breaks = pretty(
      c(0, max(data$GeneRatio_numeric, na.rm = TRUE)),
      n = 3  ##### To control x-axis ticks-number
    ),
    expand = expansion(mult = c(0.02, 0.08))
  ) +
  labs(
    title = "KEGG",
    x = "Gene ratio",
    y = "",
    size = "Protein count",
    color = "-log10(adjusted P-value)"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 16, family = "sans", face = "bold", color = "black"),
    axis.text.x = element_text(size = 14, family = "sans", face = "bold", color = "black"),
    axis.title = element_text(size = 12, family = "sans", face = "bold", color = "black"),
    plot.title = element_text(size = 14, family = "sans", face = "bold", color = "black", hjust = 0.5),
    legend.title = element_text(size = 14, family = "sans", face = "bold", color = "black"),
    legend.text = element_text(size = 12, family = "sans", color = "black"),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.ticks.length = unit(0.1, "cm"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  ) +
  guides(
    size = guide_legend(order = 1),
    color = guide_colorbar(order = 2)
  )

# Display plot
print(bubble_chart)

# Save as high-resolution TIFF
tiff(
  filename = "xxxx.tif",
  width = 9.2,
  height = 6,
  units = "in",
  res = 600,
  compression = "lzw"
)

print(bubble_chart)

dev.off()

