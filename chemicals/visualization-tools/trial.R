setwd("~/Documents/chemicals/visualization-tools/")

packages <- c("palmerpenguins", "ochRe", "GGally", "tourr", "liminal")
package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)

ggscatmat(penguins, 
          columns = 3:6, 
          col="species") +
  scale_colour_ochre(
    palette="nolan_ned") +
  theme(aspect.ratio=1,
        legend.position="bottom")

clrs <- ochre_pal(
  palette="nolan_ned")(3)
col <- clrs[
  as.numeric(
    penguins$species)]
animate_xy(penguins[,3:6], 
           col=col, 
           axes="off", 
           fps=15)

animate_xy(flea[,-7], col=flea$species)

limn_tour(penguins, cols = 3:6)
