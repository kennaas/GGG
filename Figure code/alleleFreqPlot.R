library(knitr)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(grid)

loter_run = 1

load(file = paste0("Runs/AlleleFreqs/alleleFreqs", loter_run, ".RData"))
#head(comp)
#kable(comp)

my_colors = brewer.pal(3, "Dark2")[1:3]


ggplot(data = data.frame(pInner, pOuter, pOther)) +
  stat_density(aes(x = pInner, color = "inner"),
               geom = "line", size = 1) + 
  stat_density(aes(x = pOuter, color = "outer"),
               geom = "line", size = 1) + 
  stat_density(aes(x = pOther, color = "other"),
               geom = "line", size = 1) +
  theme_light() +
  scale_color_manual(values = my_colors, 
                     breaks = c("inner", "outer", "other")) +
  ylab("Density") +
  xlab(bquote(paste("Group-specific allele frequencies ", italic(hat(p)[mr])))) +
  theme(panel.border = element_rect(color = "black", size = 1),
        panel.grid = element_blank(),
        plot.title = element_text(size = 10, face = "bold"),
        legend.position = c(0.85, 0.7),
        legend.spacing.y = unit(0, "cm"),
        legend.key.height = unit(-0.02, "cm"),
        legend.key.width = unit(.88, "cm"),
        legend.box.background = element_rect(color = "black",
                                             linetype = "solid",
                                             size = 0.5),
        legend.background = element_blank(),
        legend.margin = margin(0.1,0.1,0.13,0.1, unit="cm")
  ) +
  guides(color = guide_legend(title = NULL))

dev.off()
ggsave(path = "Figures",
       filename = paste0("alleleFreqPlot", loter_run, ".pdf"),
       device = "pdf", units = "in",
       width = 4, height = 3)
rm(list = ls())
