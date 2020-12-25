library(knitr)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(grid)
load(file = "Runs/morphData_pedigree_version.RData")
load(file = "Runs/q&pi_comp3.RData")
load(file = "Data/Founders3/founders.RData")
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
    ggtitle(paste0(group, ": ", ifelse(admOnly, 
                                       "admixed population only",
                                       "phenotyped population"))) + #"Group membership proportions: ", group), 
    #       subtitle = "Genome-based vs. pedigree-based") +
    ylab(bquote(pi[i*.(groupNum)])) + #ylab(expression(pi[iparse(groupNum)])) +
    xlab(bquote(q[i*.(groupNum)])) +
    guides(color = NULL) +
    #geom_smooth(method = "auto", color = col, 
    #            se = TRUE, fullrange = FALSE, level = 0.95) +
    theme(panel.border = element_rect(color = "black", size = 1),
          panel.grid = element_blank(),
          plot.title = element_text(size = 10, face = "plain"))
  
  
  return(qVsPiPlot)
}


scatterplot_big = function(pop, comp) {
  qVsPiPlot = ggplot(data = comp[comp$ringnr %in% pop, ]) +
    geom_point(aes(x = inner, y = pi_inner, color = "inner")) +
    geom_point(aes(x = outer, y = pi_outer, color = "outer")) +
    geom_point(aes(x = other, y = pi_other, color = "other")) +
    theme_light() +
    ggtitle("Group membership proportions", 
            subtitle = "Genome-based vs. pedigree-based") +
    scale_color_manual(values = my_colors, 
                       breaks = c("inner", "outer", "other")) +
    ylab(expression(pi[ir])) +
    xlab(expression(q[ir]))
  
  return(qVsPiPlot)
}

pop_corr = function(pop, comp) {
  comp_temp = comp[comp$ringnr %in% pop, ]
  scatterplot_big(pop, comp)
  
  print(dim(comp_temp)[1])
  
  c(cor(comp_temp$inner, comp_temp$pi_inner), 
    cor(comp_temp$outer, comp_temp$pi_outer), 
    cor(comp_temp$other, comp_temp$pi_other)) 
}

pop = morphData$ringnr
popAdm = admixedInds[admixedInds %in% morphData$ringnr]

pop_corr(pop, comp)
#scatterplot_big(pop, comp)
inner = scatterplot_single(admOnly = FALSE, comp, "inner")
outer = scatterplot_single(admOnly = FALSE, comp, "outer")
other = scatterplot_single(admOnly = FALSE, comp, "other")

innerAdm = scatterplot_single(admOnly = TRUE, comp, "inner")
outerAdm = scatterplot_single(admOnly = TRUE, comp, "outer")
otherAdm = scatterplot_single(admOnly = TRUE, comp, "other")

title = textGrob("Group membership proportions", 
                 gp = gpar(fontsize = 15, fontface = 2L))
subtitle = textGrob("Genome-based vs. pedigree-based", gp = gpar(fontsize = 10, fontface = 3L))
margin = unit(0.1, "line")
scatterPlotFull = 
  grid.arrange(title,
               subtitle,
               inner, innerAdm,
               outer, outerAdm,
               other, otherAdm,
               layout_matrix = matrix(c(1,1,1,
                                        2,2,2,
                                        3,5,7,
                                        4,6,8), byrow = TRUE,
                                      nrow = 4, ncol = 3),
               heights = c(margin, margin, 1, 1))

dev.off()
ggsave(plot = scatterPlotFull, path = "Figures", 
       filename = "groupScatterplot.pdf",
       device = "pdf", units = "in", 
       width = 8.27, height = 6)

rm(list = ls())
