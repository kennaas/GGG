library(knitr)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(grid)

loter_run = 1

load(file = "Runs/morphData_pedigree_version.RData")
load(file = paste0("Runs/q&pi_comp", loter_run, ".RData"))
load(file = paste0("Data/Founders", loter_run, "/founders.RData"))
#head(comp)
#kable(comp)

my_colors = brewer.pal(3, "Dark2")[1:3]

scatterplot_single = function(admOnly, comp, group) {
  if (admOnly) {
    pop = admixedInds[admixedInds %in% morphData$ringnr]
  } else {
    pop = morphData$ringnr
  }
  groupNum = switch(group,
                    inner = 1,
                    outer = 2,
                    other = 3)
  col = switch(group, 
               inner = my_colors[1],
               outer = my_colors[2],
               other = my_colors[3])
  qVsPiPlot = ggplot(data = comp[comp$ringnr %in% pop, ],
                     aes(x = get(group), 
                         y = get(paste0("pi_", group)))) +
    geom_point(color = col, alpha = 0.2) +
    xlim(0, 1) +
    ylim(0, 1) +
    theme_light() +
    ggtitle(group) +
    # ggtitle(paste0(group, ": ", ifelse(admOnly, 
    #                                    "admixed population only",
    #                                    "phenotyped population"))) + #"Group membership proportions: ", group), 
    #       subtitle = "Genome-based vs. pedigree-based") +
    ylab(bquote(hat(italic(pi))[italic(i)*.(groupNum)])) + #ylab(expression(pi[iparse(groupNum)])) +
    xlab(bquote(italic(q)[italic(i)*.(groupNum)])) +
    guides(color = NULL) +
    #geom_smooth(method = "auto", color = col, 
    #            se = TRUE, fullrange = FALSE, level = 0.95) +
    theme(panel.border = element_rect(color = "black", size = 1),
          panel.grid = element_blank(),
          plot.title = element_text(size = 10, face = "plain"))
  
  
  return(qVsPiPlot)
}

pop_corr = function(pop, comp) {
  comp_temp = comp[comp$ringnr %in% pop, ]

  
  print(dim(comp_temp)[1])
  
  c(cor(comp_temp$inner, comp_temp$pi_inner), 
    cor(comp_temp$outer, comp_temp$pi_outer), 
    cor(comp_temp$other, comp_temp$pi_other)) 
}

pop = morphData$ringnr
popAdm = admixedInds[admixedInds %in% morphData$ringnr]

pop_corr(pop, comp)
pop_corr(popAdm, comp)

inner = scatterplot_single(admOnly = FALSE, comp, "inner")
outer = scatterplot_single(admOnly = FALSE, comp, "outer")
other = scatterplot_single(admOnly = FALSE, comp, "other")

innerAdm = scatterplot_single(admOnly = TRUE, comp, "inner")
outerAdm = scatterplot_single(admOnly = TRUE, comp, "outer")
otherAdm = scatterplot_single(admOnly = TRUE, comp, "other")

title = textGrob("Group membership proportions (phenotyped admixed population)", 
                 gp = gpar(fontsize = 15, fontface = 2L))
subtitle = textGrob("Genome-based vs. pedigree-based", gp = gpar(fontsize = 10, fontface = 3L))
margin = unit(0.1, "line")
# scatterPlotFull = 
#   grid.arrange(title,
#                subtitle,
#                inner, innerAdm,
#                outer, outerAdm,
#                other, otherAdm,
#                layout_matrix = matrix(c(1,1,1,
#                                         2,2,2,
#                                         3,5,7,
#                                         4,6,8), byrow = TRUE,
#                                       nrow = 4, ncol = 3),
#                heights = c(margin, margin, 1, 1))

scatterPlotFull = 
  grid.arrange(title,
               subtitle,
               innerAdm, outerAdm, otherAdm,
               layout_matrix = matrix(c(1,1,1,
                                        2,2,2,
                                        3,4,5), byrow = TRUE,
                                      nrow = 3, ncol = 3),
               heights = c(margin, margin, 1))

dev.off()
ggsave(plot = scatterPlotFull, path = "Figures", 
       filename = paste0("groupScatterplot", loter_run, ".pdf"),
       device = "pdf", units = "in", 
       width = 8.27, height = 3.3)

rm(list = ls())
