library(ggplot2)
library(raster)
library(tidyr)
library(dplyr)
library(sp)
library(rgdal)
library(Hmisc)

# cargar datos
data <- read.table("file_list_proportion.csv", sep = ",", header = TRUE, 
  stringsAsFactors = FALSE)

# data structure
head(data)


# conglomerate variable
cgl <- unlist(strsplit(data$archivo,split="\\-HEM"))[2*(1:nrow(data))-1]
# date variable
mes <- as.numeric(substring(data$archivo, nchar(data$archivo) - 7, 
  nchar(data$archivo) - 6))

# sites variable
site <- sub("^[0-9]+-HEM-(S[0-9]).*", "\\1", data$archivo) 

# add new variables to data.frame
data$Cgl <- as.character(cgl)
data$site <- as.character(site)
data$season <- ifelse(mes %in% c(12, 1:4), "Dry", "Wet")

# eliminate duplicates
data_sub <- data[!duplicated(data[,c("Cgl","site")]),
  c("Cgl","site","proportion", "season", "comparacion")]

# create variables by plot (cgl)
prop_cgl <- data_sub %>%
  group_by(Cgl) %>%
  summarise(
    file_comp = first(comparacion),
    season = first(season),
    mean_prop = mean(proportion),
    range_prop = max(proportion) - min(proportion)
  )
  
# add coordinates to plots (cgl)
data_coords <- read.table("../data/congs_coords.csv",sep=",",header=TRUE)
coords <- data_coords
head(coords)
coordinates(coords) <- ~ x + y

# we will get the Holdrige class for each coordinate
raster_zvh <- raster("../data/zvh_31_lcc.tif")
projection(coords) <- projection(raster_zvh)
plot(raster_zvh)
points(coords)
data_coords$zvh <- raster::extract(raster_zvh, coords)

# we repeat for a coarser classification of zvh
raster_zvh_c <- raster("../data/zvh_p_hgw.tif")
# the raster is in a different projection
projection(raster_zvh_c)
# reproject
raster_zvh_rep <- projectRaster(raster_zvh_c, crs = projection(raster_zvh), 
  method = "ngb")
plot(raster_zvh_rep)
points(coords)

data_coords$zvh_c <- raster::extract(raster_zvh_rep, coords)
table(data_coords$zvh_c)

# we do the same for ecological integrity
raster_ei <- raster("/Volumes/ARCHIVOS_MAC/Dropbox/Datos Redes Bayesianas/EI_maps/Stage_3/Final_net_Scores.tif")
plot(raster_ei)
points(coords)
data_coords$ei <- raster::extract(raster_ei, coords)

data_coords$Cgl <- as.character(coords$Cgl)

prop_coords <- prop_cgl %>% 
  left_join(data_coords) %>%
  mutate(zvh_c = factor(zvh_c)) %>%
  group_by(zvh_c) %>%
  mutate(
    n = n() 
    ) %>%
  filter(n > 30, !is.na(zvh_c)) %>%
  mutate(
    prop_class = cut(mean_prop, 
      breaks = quantile(mean_prop, probs = c(0, 0.333, 0.666, 1), na.rm = TRUE),
      include.lowest = TRUE, labels = c("low", "medium", "high")),
    ei_class = cut2(ei, cuts = c(0, 0.5, 0.75, 1))
  )

head(prop_coords)

ggplot(prop_coords, aes(x = mean_prop, y = ei)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~ zvh_c, nrow = 2)

head(prop_coords)

table(prop_coords$prop_class, prop_coords$ei_class, prop_coords$zvh_c)

round(prop.table(table(prop_coords$prop_class, prop_coords$ei_class), 1) * 100, 1)

fit <- lm(ei ~ mean_prop + range_prop + zvh_c + mean_prop * zvh_c + 
    range_prop * zvh_c, data = prop_coords)
summary(fit)

# almost all observations come from Wet season
table(prop_coords$season)


##### Shape 
file_coords <- prop_cgl %>%
  left_join(data_coords) %>%
  filter(!is.na(x), !is.na(y)) %>%
  as.data.frame()

coordinates(file_coords) <- ~ x + y
projection(file_coords) <- projection(raster_zvh)

writeOGR(file_coords, "./hemispheric_points", "hemispheric_pics", 
  driver = "ESRI Shapefile", verbose = FALSE, overwrite_layer = TRUE)

## selección de conglomerados Pedro
sel_shp <- readOGR("fotosHemisf_p" , "selP") %>%
  as.data.frame() 
sel_files <- sel_shp %>%
  mutate(file_comp = as.character(file_comp)) %>%
  select(comparacion = file_comp) %>%
  left_join(data)

# write filenames to a text file named selected files
write.table(sel_files, row.names = FALSE, sep = ",")

# Copy images from ameyalli (access in Julián's computer)
selected_files <- read.csv("~/Downloads/imagenes_zvh/selected_files.txt")
setwd("Q:/FOTOGRAFIAS HEMISFERICAS/2012")
exito <- file.copy(from = paste("Q:/FOTOGRAFIAS HEMISFERICAS/2012/", 
  as.character(selected_files$archivo), sep = ""), 
  to = "Q:/FOTOGRAFIAS HEMISFERICAS/2012/seleccion_zvh")
