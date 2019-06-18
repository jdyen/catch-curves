# "observed" is observed abundances (years in rows, age classes in columns)
# "fitted" is modelled abundances (years in rows, age classes in columns)
# "residuals" are the observed minus modelled abundances (years in rows, age classes in columns)

par(mfrow = c(3, 1), mar = c(4.5, 4.5, 3.1, 1.1))

# plot observed values
abs_range <- max(abs(range(observed, na.rm = TRUE)))
image_breaks <- c(seq(0, abs_range, length = 100), abs_range + 0.0001)
col_pal_image <- colorRampPalette(c("#f7f7f7",
                                    "#2166ac", "#053061"))(99)
col_pal_image <- c(col_pal_image, ggplot2::alpha("black", 0.7))
plot_values <- ifelse(is.na(observed), abs_range + 0.00005, observed)
fields::image.plot(plot_values,
                   col = col_pal_image,
                   breaks = image_breaks,
                   xaxt = "n", yaxt = "n",
                   xlab = "Year", ylab = "Age")
axis(1, at = seq(0.05, 1, length = 10), labels = seq(2000, 2018, by = 2))
axis(2, at = seq(0, 1, length = n_obs_tmp), labels = seq(0, n_obs_tmp - 1, by = 1),
     las = 1)
mtext("Observed abundances", side = 3, line = 1, adj = 1, cex = 1.1)

# plot fitted values
abs_range <- max(abs(range(fitted, na.rm = TRUE)))
image_breaks <- c(seq(0, abs_range, length = 100), abs_range + 0.0001)
col_pal_image <- colorRampPalette(c("#f7f7f7",
                                    "#2166ac", "#053061"))(99)
col_pal_image <- c(col_pal_image, ggplot2::alpha("black", 0.7))
plot_values <- ifelse(is.na(fitted), abs_range + 0.00005, fitted)
fields::image.plot(plot_values,
                   col = col_pal_image,
                   breaks = image_breaks,
                   xaxt = "n", yaxt = "n",
                   xlab = "Year", ylab = "Age")
axis(1, at = seq(0.05, 1, length = 10), labels = seq(2000, 2018, by = 2))
axis(2, at = seq(0, 1, length = n_obs_tmp), labels = seq(0, n_obs_tmp - 1, by = 1),
     las = 1)
mtext("Modelled abundances", side = 3, line = 1, adj = 1, cex = 1.1)

# plot residuals
abs_range <- max(abs(range(residuals, na.rm = TRUE)))
image_breaks <- c(seq(-abs_range, abs_range, length = 100), abs_range + 0.0001)
col_pal_image <- colorRampPalette(c("#67001f", "#b2182b",
                                    "#f7f7f7",
                                    "#2166ac", "#053061"))(99)
col_pal_image <- c(col_pal_image, ggplot2::alpha("black", 0.7))
plot_values <- ifelse(is.na(residuals), abs_range + 0.00005, residuals)
fields::image.plot(plot_values,
                   col = col_pal_image,
                   breaks = image_breaks,
                   xaxt = "n", yaxt = "n",
                   xlab = "Year", ylab = "Age")
axis(1, at = seq(0.05, 1, length = 10), labels = seq(2000, 2018, by = 2))
axis(2, at = seq(0, 1, length = n_obs_tmp), labels = seq(0, n_obs_tmp - 1, by = 1),
     las = 1)
mtext("Residual (age class strength)", side = 3, line = 1, adj = 1, cex = 1.1)
