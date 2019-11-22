# plot variation in flow variables used in analyses

# load required data
flow_compiled <- readRDS("data/compiled-flow-to-plot.rds")
sys_year <- readRDS("data/system-year-lookup.rds")

# set system names for plots
systems_to_keep <- c("broken", "goulburn", "king", "murray", "ovens")
system_order <- c(4, 5, 2, 3, 1)
colour_palette_order <- c(5, 3, 4, 1, 2)
system_names <- systems_to_keep
for (i in seq_along(system_names)) {
  init <- toupper(substr(system_names[i], 1, 1))
  remainder <- substr(system_names[i], 2, nchar(system_names[i]))
  system_names[i] <- paste0(init, remainder)
}

# convert from rebased years to actual years for plotting
actual_years <- sys_year$year + 1992

# which systems do we want to plot?
sys_to_plot <- c(1:5)

# which variables?
vars_to_plot <- colnames(flow_compiled)

# match variable codes with full names
var_names <- c("spawning_variability" = "Proportional variability in spawning flows",
               "prop_spring_lt" = "Spring flow relative to long-term flow",
               "prop_summer_lt" = "Summer flow relative to long-term flow",
               "prop_winter_lt" = "Winter flow relative to long-term flow",
               "prop_max_antecedent_lt" = "Antecedent maximum flow relative to long-term flow",
               "spawning_temp" = "Temperature during spawning months",
               "adult_cpue" = "Adult CPUE")
var_names <- var_names[vars_to_plot]

# plot flow variation for each system
jpeg(file = "outputs/plots/flow_variation.jpg",
     width = 8.2, height = 10, units = "in", res = 300)
layout(matrix(1:7, ncol = 1), heights = c(rep(1, 6), 0.45))
par(mar = c(3.3, 4.0, 1.5, 1.1))
col_pal <- viridis::inferno(256)[seq(1, 250, length = max(sys_to_plot))]
for (i in seq_along(vars_to_plot)) {
  
  test_var <- vars_to_plot[i]
  
  for (j in seq_along(unique(sys_to_plot))) {
    
    idx <- sys_year$system == j
    
    yplot <- flow_compiled[idx, test_var]
    xplot <- actual_years[idx]
    
    if (j == 1) {
      plot(yplot ~ xplot,
           pch = 16, bty = "l", las = 1, type = "n",
           xlab = "", ylab = "",
           xaxt = "n",
           xlim = range(actual_years),
           ylim = range(flow_compiled[, test_var], na.rm = TRUE))
      year_labels <- seq(min(actual_years), max(actual_years), by = 2)
      axis(1, at = year_labels, labels = year_labels)
      mtext("Year", side = 1, line = 2.1, adj = 0.5, cex = 0.8)
      mtext("Value", side = 2, line = 2.9, adj = 0.5, cex = 0.8)
      mtext(var_names[i], side = 3, line = -0.05, adj = 0.02, cex = 0.9)
    }
    
    lines(yplot ~ xplot,
          col = col_pal[colour_palette_order[j]], lwd = 3)
    
  } 
  
}

# add a legend in the final panel
plot(1, 1, type = "n", bty = "n", xaxt = "n", yaxt = "n",
     xlab = "", ylab = "")
legend(legend = system_names[system_order],
       x = "center",
       cex = 1.25,
       horiz = TRUE,
       lty = 1, lwd = 2, bty = "n",
       col = col_pal,
       xpd = TRUE)

dev.off()
