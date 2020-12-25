library(ggplot2)
library(gridExtra)
library(INLA)
library(RColorBrewer)
save = TRUE
run_num = 3


load(paste0("Runs/wing/inla_result_hetGG_rio", run_num, "_VR_wing.RData"))
wingRioVR = model
load("Runs/wing/inla_result_hetGG_pedigree_wing.RData")
wingGGPed = model

load(paste0("Runs/mass/inla_result_hetGG_rio", run_num, "_VR_mass.RData"))
massRioVR = model
load("Runs/mass/inla_result_hetGG_pedigree_mass.RData")
massGGPed = model

load(paste0("Runs/tarsus/inla_result_hetGG_rio", run_num, "_VR_tarsus.RData"))
tarsusRioVR = model
load("Runs/tarsus/inla_result_hetGG_pedigree_tarsus.RData")
tarsusGGPed = model

my_colors = brewer.pal(3, "Dark2")[1:3]

if (save) { 
  path = paste0("Figures/VA_rio", run_num, ".pdf")
  #pdf(file = path, width = 8.27, height = 9)
}

get_posterior = function(model, precision) {
  return(inla.smarginal(inla.tmarginal(
    function(x) 1 / x,
    model$marginals.hyperpar[[precision]])))
}

make_post_plot = function(response, legend) {
  innerVarGen = get_posterior(get(paste0(response, "RioVR")), 5)
  outerVarGen = get_posterior(get(paste0(response, "RioVR")), 6)
  otherVarGen = get_posterior(get(paste0(response, "RioVR")), 7)
  innerVarPed = get_posterior(get(paste0(response, "GGPed")), 5)
  outerVarPed = get_posterior(get(paste0(response, "GGPed")), 6)
  otherVarPed = get_posterior(get(paste0(response, "GGPed")), 7)
  
  plot = ggplot() +
    geom_line(data = subset(data.frame(innerVarGen), y > 0.02),
              aes(x, y, color = "inner", linetype = "Genomic"),
              size = 1) +
    geom_line(data = subset(data.frame(outerVarGen), y > 0.02), 
              aes(x, y, color = "outer", linetype = "Genomic"),
              size = 1) +
    geom_line(data = subset(data.frame(otherVarGen), y > 0.02),
              aes(x, y, color = "other", linetype = "Genomic"),
              size = 1) +
    geom_line(data = subset(data.frame(innerVarPed), y > 0.02),
              aes(x, y, color = "inner", linetype = "Pedigree"),
              size = 1) +
    geom_line(data = subset(data.frame(outerVarPed), y > 0.02),
              aes(x, y, color = "outer", linetype = "Pedigree"),
              size = 1) +
    geom_line(data = subset(data.frame(otherVarPed), y > 0.02),
              aes(x, y, color = "other", linetype = "Pedigree"),
              size = 1) +
    #  xlim(0, 4) +
    theme_light() +
    ggtitle(switch(response, 
                   wing = "Wing length",
                   mass = "Body mass",
                   tarsus = "Tarsus length"))+ 
    scale_color_manual(values = my_colors, 
                       breaks = c("inner", "outer", "other")) +
    scale_linetype_manual(values = c("solid", "dotted")) +
    ylab("Density") +
    xlab(expression(sigma[G[r]]^{2})) +
    scale_y_continuous(labels = function(x) {sprintf("%.1f", x)}) +
    theme(panel.border = element_rect(color = "black", size = 1),
          panel.grid = element_blank(),
          legend.position = c(0.85, 0.7),
          legend.spacing.y = unit(0, "cm"),
          #legend.box = "vertical",
          legend.key.height = unit(-0.06, "cm"),
          legend.key.width = unit(.88, "cm"),
          legend.box.background = element_rect(color = "black",
                                               linetype = "solid",
                                               size = 0.5),
          legend.background = element_blank(),
          legend.margin = margin(0.1,0.1,0.13,0.1, unit="cm"))
  
  
  if (legend) {
    plot = plot + guides(color = guide_legend(title = NULL), 
                         linetype = guide_legend(title = NULL))
  } else {
    plot = plot + guides(color = FALSE, 
                         linetype = FALSE)
  }
  
  return(plot)
}

wingPlot = make_post_plot("wing", legend = TRUE)
massPlot = make_post_plot("mass", legend = FALSE)
tarsusPlot = make_post_plot("tarsus", legend = FALSE)

allPlot = grid.arrange(wingPlot, massPlot, tarsusPlot, nrow = 3)

if (save) {
  dev.off()
  ggsave(plot = allPlot, path = "Figures", filename = "VA_rio3.pdf",
         device = "pdf", units = "in", 
         width = 8.27, height = 9)
  #system(paste0('start "', path, '"'))
} else {
  rm(list = ls())
}