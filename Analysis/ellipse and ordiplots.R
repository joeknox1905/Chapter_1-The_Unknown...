#This script creates plots of migrants vs residents in trait space using ellipse volumes and 
performs a PERMANOVA BY migratory status, ecoregion and season as a preliminary analysis akin to distance between centroids

load libraries------------
library(vegan)
Library(ggplot)


#loading traitspaces
avospace<- read.csv("Data/Processed/eco_region_traitspace__PCA_BI.csv")

#@@@ JK n.b. will also want to join pcs to assemb for seasonal splits and analysis e.g. BI seasonal traitspace from avospace

#######################################################################
#visualising overlap of groups in traitspace with ggplot------------------
#######################################################################

ellipse_plot1 <- ggplot(avospace, aes(PC1, PC2, colour = status)) +
  stat_ellipse(aes(fill = status), geom = "polygon", alpha = 0.2) +
  geom_point() +
  facet_wrap(~ season + ecoregion) +
  theme_light()

ggsave("ecoregion_seasonal_status_pca_ellipses_plot1.png", plot = ellipse_plot1, dpi=300)

ellipse_plot2 <- ggplot(avospace, aes(PC3, PC4, colour = status)) +
  stat_ellipse(aes(fill = status), geom = "polygon", alpha = 0.2) +
  geom_point() +
  facet_wrap(~ season + ecoregion) +
  theme_light() 

  ggsave("ecoregion_seasonal_status_pca_ellipses_plot2.png", plot = ellipse_plot2, dpi=300)

ellipse_plot3 <- ggplot(avospace, aes(PC5, PC6, colour = status)) +
  stat_ellipse(aes(fill = status), geom = "polygon", alpha = 0.2) +
  geom_point() +
  facet_wrap(~ season + ecoregion) +
  theme_light() 

  ggsave("ecoregion_seasonal_status_pca_ellipses_plot3.png", plot = ellipse_plot3, dpi=300)

####################################################################
#visualise distance to centroids of group with vegan---------------------
######################################################################


#@@@JK will need to save plots if want to keep

avospace <- avospace %>% mutate(status = factor(status))
pal <- setNames(hcl.colors(nlevels(avospace$status), "Dark 3"), levels(avospace$status))

panels <- expand.grid(ecoregion = sort(unique(avospace$ecoregion)),
                      season = sort(unique(avospace$season)),
                      stringsAsFactors = FALSE)

par(mfrow = c(n_distinct(avospace$season), n_distinct(avospace$ecoregion)),
    mar = c(3, 3, 2, 1))

for (i in seq_len(nrow(panels))) {
  d <- avospace %>% filter(ecoregion == panels$ecoregion[i], season == panels$season[i])
  plot(range(avospace$PC1), range(avospace$PC2), type = "n", las = 1,
       main = paste(panels$ecoregion[i], panels$season[i], sep = " | "), cex.main = 0.8)
  if (nrow(d) == 0) next

  xy <- as.matrix(d[, c("PC1", "PC2")])
  g  <- droplevels(d$status)

  points(xy, pch = 16, col = adjustcolor(pal[as.character(g)], 0.6))
  ordiellipse(xy, g, display = "sites", kind = "sd", conf = 0.95,
              draw = "polygon", col = pal[levels(g)], alpha = 50, label = FALSE)
  ordispider(xy, g, display = "sites", col = pal[levels(g)], lwd = 1)

  if (i == 1) legend("topright", legend = names(pal), fill = pal, bty = "n", cex = 0.7)
}

#PC 3V4:

avospace <- avospace %>% mutate(status = factor(status))
pal <- setNames(hcl.colors(nlevels(avospace$status), "Dark 3"), levels(avospace$status))

panels <- expand.grid(ecoregion = sort(unique(avospace$ecoregion)),
                      season = sort(unique(avospace$season)),
                      stringsAsFactors = FALSE)

par(mfrow = c(n_distinct(avospace$season), n_distinct(avospace$ecoregion)),
    mar = c(3, 3, 2, 1))

for (i in seq_len(nrow(panels))) {
  d <- avospace %>% filter(ecoregion == panels$ecoregion[i], season == panels$season[i])
  plot(range(avospace$PC1), range(avospace$PC2), type = "n", las = 1,
       main = paste(panels$ecoregion[i], panels$season[i], sep = " | "), cex.main = 0.8)
  if (nrow(d) == 0) next

  xy <- as.matrix(d[, c("PC3", "PC4")])
  g  <- droplevels(d$status)

  points(xy, pch = 16, col = adjustcolor(pal[as.character(g)], 0.6))
  ordiellipse(xy, g, display = "sites", kind = "sd", conf = 0.95,
              draw = "polygon", col = pal[levels(g)], alpha = 50, label = FALSE)
  ordispider(xy, g, display = "sites", col = pal[levels(g)], lwd = 1)

  if (i == 1) legend("topright", legend = names(pal), fill = pal, bty = "n", cex = 0.7)
}

# pc 5V6

avospace <- avospace %>% mutate(status = factor(status))
pal <- setNames(hcl.colors(nlevels(avospace$status), "Dark 3"), levels(avospace$status))

panels <- expand.grid(ecoregion = sort(unique(avospace$ecoregion)),
                      season = sort(unique(avospace$season)),
                      stringsAsFactors = FALSE)

par(mfrow = c(n_distinct(avospace$season), n_distinct(avospace$ecoregion)),
    mar = c(3, 3, 2, 1))

for (i in seq_len(nrow(panels))) {
  d <- avospace %>% filter(ecoregion == panels$ecoregion[i], season == panels$season[i])
  plot(range(avospace$PC1), range(avospace$PC2), type = "n", las = 1,
       main = paste(panels$ecoregion[i], panels$season[i], sep = " | "), cex.main = 0.8)
  if (nrow(d) == 0) next

  xy <- as.matrix(d[, c("PC5", "PC6")])
  g  <- droplevels(d$status)

  points(xy, pch = 16, col = adjustcolor(pal[as.character(g)], 0.6))
  ordiellipse(xy, g, display = "sites", kind = "sd", conf = 0.95,
              draw = "polygon", col = pal[levels(g)], alpha = 50, label = FALSE)
  ordispider(xy, g, display = "sites", col = pal[levels(g)], lwd = 1)

  if (i == 1) legend("topright", legend = names(pal), fill = pal, bty = "n", cex = 0.7)
}


