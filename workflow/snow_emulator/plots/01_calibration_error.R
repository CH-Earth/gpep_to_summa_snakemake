
################################################################################
#
# Script to plot the errors of the emulators for snow
#
# nicolas.vasquez@ucalgary.ca
#
################################################################################

library("ggplot2")
library("dplyr")
library("sf")
library("patchwork")
library("cowplot")

################################################################################
# Variables required by the script to run
################################################################################

dirSave   <- file.path( "/Users/nvasquez/Library/CloudStorage/Dropbox-Personal/UofC/03_manuscripts/02_Dave" )

load( file.path( dirSave, "my_environment.RData"))

dirSave   <- file.path( "/Users/nvasquez/Library/CloudStorage/Dropbox-Personal/UofC/03_manuscripts/02_Dave" )


# Load Shapefiles
na_shp <- st_read("/Users/nvasquez/Library/CloudStorage/Dropbox-Personal/UofC/01_Databases/NorthAmerica/boundaries_p_2021_v3.shp")
bd_shp <- st_read("/Users/nvasquez/Library/CloudStorage/Dropbox-Personal/UofC/01_Databases/WB_Adm0_boundary/WB_Adm0_boundary_lines_10m.gpkg")
ci_shp <- st_read("/Users/nvasquez/Library/CloudStorage/Dropbox-Personal/UofC/01_Databases/Caribe_Islands.gpkg")

na_shp <- st_simplify(x = na_shp, preserveTopology = FALSE, dTolerance = 10000)
na_shp <- st_transform(x = na_shp, crs = "epsg:4326")

################################################################################
# Maps
################################################################################

fig.map.a <- ggplot() +
                  geom_sf(data = na_shp, fill = "gray85", colour = "gray88") +
                  geom_sf(data = bd_shp, colour = "gray30", linetype = "dashed", linewidth = 0.3) +
                  geom_sf(data = ci_shp, fill = "gray85", colour = "gray88") +
                  geom_point(data   = stationsReg, aes(x = lon, y = lat), stroke = 0.1, colour = "black", fill = "gray50", shape = 21, size = 1) +
                  geom_point(data   = stations,    aes(x = lon, y = lat), stroke = 0.1, colour = "black", fill = "red", shape = 21, size = 1.5) +
                  scale_x_continuous(name = "Longitude") +
                  scale_y_continuous(name = "Latitude") +
                  coord_sf(datum       = sf::st_crs(4326),
                           default_crs = sf::st_crs(4326),
                           crs         = sf::st_crs("EPSG:9311"),
                           xlim        = c(-140, -53),
                           ylim        = c(26.5, 75)) +
                  labs(title    = "a) Location of stations used for training",
                       subtitle = "(in red)") +
                  theme(plot.background        = element_blank(),
                        panel.background       = element_blank(),
                        panel.border           = element_rect(fill = NA, colour = "black",  linewidth = 0.6),
                        panel.grid.major       = element_line(           colour = "gray80", linewidth = 0.1),
                        axis.ticks             = element_line(           colour = "black",  linewidth = 0.2),
                        axis.ticks.length      = unit(0.15, "cm"),
                        legend.position        = "bottom",
                        strip.background       = element_blank(),
                        text                   = element_text(colour = "black")) 

fig.map.b <- ggplot() +
              geom_sf(data = na_shp, fill = "gray85", colour = "gray88") +
              geom_sf(data = bd_shp, colour = "gray30", linetype = "dashed", linewidth = 0.3) +
              geom_sf(data = ci_shp, fill = "gray85", colour = "gray88") +
              geom_point(data   = stations,    aes(x = lon, y = lat), stroke = 0.1, colour = "black", fill = "gray50", shape = 21, size = 1) +
              geom_point(data   = stationsReg, aes(x = lon, y = lat), stroke = 0.1, colour = "black", fill = "red",    shape = 21, size = 1.5) +
              labs(title    = "c) Location of stations used for regionalization",
                   subtitle = "(in red)") +
              scale_x_continuous(name = "Longitude") +
              scale_y_continuous(name = "Latitude") +
              coord_sf(datum       = sf::st_crs(4326),
                       default_crs = sf::st_crs(4326),
                       crs         = sf::st_crs("EPSG:9311"),
                       xlim        = c(-140, -53),
                       ylim        = c(26.5, 75)) +
              theme(plot.background        = element_blank(),
                    panel.background       = element_blank(),
                    panel.border           = element_rect(fill = NA, colour = "black",  linewidth = 0.6),
                    panel.grid.major       = element_line(           colour = "gray80", linewidth = 0.1),
                    axis.ticks             = element_line(           colour = "black",  linewidth = 0.2),
                    axis.ticks.length      = unit(0.15, "cm"),
                    legend.position        = "bottom",
                    strip.background       = element_blank(),
                    text                   = element_text(colour = "black")) 

################################################################################
#                         Compute calibration error
################################################################################

fig.1 <- ggplot( data = plotData[rows, ], aes( x = TSSv1, y = targetMetricRF ) ) +
          geom_hex( bins = 200 ) +
          geom_abline( slope = 1, intercept = 0, colour = "black", linetype = "dashed", linewidth = 0.1 ) +
          scale_x_continuous( name = "TSS from SUMMA",
                              breaks = seq(from = 0, to = 1, by = 0.2)) +
          scale_y_continuous( name = "TSS from Random Forest Emulator" ,
                              breaks = seq(from = 0, to = 1, by = 0.2)) +
          scale_fill_viridis_c( name = "Count", option = "plasma", trans = "log10" ) +
          labs(title    = "b) Large Sample Emulator",
               subtitle = "Training error for i = 5") +
          annotate( "text",
                    x     = -0.01, y    = 1.05,
                    hjust = 0,    vjust = 1,
                    label    = paste0( "Training error:\nRMSE = ", 
                                       round(x = errorsCalib$rmse[  rowError], digits = 3), " (TSS units)\nMAE  = ", 
                                       round(x = errorsCalib$mae[   rowError], digits = 3), " (TSS units)\nEvaluation error:\nRMSE = ",
                                       round(x = errorsSearch$rmse[ 5],        digits = 3), " (TSS units)\nMAE  = ", 
                                       round(x = errorsSearch$mae[  5],        digits = 3), " (TSS units)"),
                    size     = 3.5,
                    colour   = "black") +
          theme( panel.background = element_blank(),
                 panel.border     = element_rect(fill = NA, colour = "black", linewidth = 0.5),
                 panel.grid       = element_line(colour = "gray80", linewidth = 0.2, linetype = "dotted"),
                 axis.ticks       = element_line( colour = "black", linewidth = 0.4),
                 axis.ticks.length = unit(0.15, "cm"),
                 legend.position   = "inside",
                 legend.position.inside = c(0.85, 0.3),
                 axis.text              = element_text(size = 10))


# Regionalization
errorsReg <- group_by( .data = plotDataReg, step)
errorsReg <- summarize( .data = errorsReg, 
                           rmse  = sqrt( mean( (TSSv1 - emulatedMetric)^2 , na.rm = TRUE)),
                           mae   = mean( abs(   TSSv1 - emulatedMetric)   , na.rm = TRUE))

aux.1 <- data.frame(step = errorsReg$step, value = errorsReg$rmse, metric = "RMSE")
aux.2 <- data.frame(step = errorsReg$step, value = errorsReg$mae,  metric = "MAE")

aux <- rbind(aux.1, aux.2)

fig.2 <- ggplot( data = aux) +
          geom_line(  aes(x = step, y = value, colour = metric), linetype = "dashed") +
          geom_point( aes(x = step, y = value, colour = metric), size = 2) +
          scale_x_continuous(name = "Iteration step") +
          scale_y_continuous(name = "Error (TSS units)") +
          scale_colour_manual(name = "Metric",
                              values = c("purple", "red")) +
          labs( title = "d) Error for unseen stations") +
          theme( panel.background = element_blank(),
                 panel.border     = element_rect(fill = NA, colour = "black", linewidth = 0.5),
                 panel.grid       = element_line(colour = "gray80", linewidth = 0.2, linetype = "dotted"),
                 axis.ticks       = element_line( colour = "black", linewidth = 0.4),
                 axis.ticks.length = unit(0.15, "cm"),
                 legend.position   = "inside",
                 legend.position.inside = c(0.85, 0.8),
                 axis.text              = element_text(size = 10))



fig <- (fig.map.a | fig.1) / (fig.map.b | fig.2) 

fig

fig <- plot_grid(fig.map.a, fig.1, 
                 fig.map.b, fig.2,
                 ncol  = 2, 
                 align = "hv",
                 axis  = "tblr")

ggsave( plot     = fig,
        filename = file.path( dirSave, "errors.png"),
        width    = 16,
        height   = 15,
        units    = "cm",
        scale    = 1.45,
        dpi      = 400)

ggsave( plot     = fig,
        filename = file.path( dirSave, "errors.pdf"),
        width    = 16,
        height   = 15,
        units    = "cm",
        scale    = 1.45)


################################################################################
#                                                            Fin. ¡Buena suerte!
################################################################################
