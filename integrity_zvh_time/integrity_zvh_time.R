# En este script buscamos relacionar la presencia de especies domésticas
# (vacas, borregos, perros, caballos) con integridad ecológica (ROBIN)

library(Hmisc)
library(rgdal)
library(sp)
library(raster)
library(tidyr)
library(plyr)
library(dplyr)
library(lubridate)

# y copiamos la proyección del raster de integridad para reproyectar a lcc
raster_anio <- raster("../data/integrity_time/out_bn_ie_tabfinal_2004_1.tif")
raster_zvh <- raster("../data/zvh_31_lcc.tif")
raster_zvh_crop <- crop(raster_zvh, raster_anio)
raster_stack <- stack(raster_anio, raster_zvh_crop)
zvh_anio <- as.data.frame(raster_stack)

ieReadSum <- function(file_path){
  raster_anio <- raster(file_path)
  raster_stack <- stack(raster_anio, raster_zvh_crop)
  zvh_anio <- as.data.frame(raster_stack)
  colnames(zvh_anio) <- c("ie", "zvh")
  zvh_summary <- zvh_anio %>%
      filter(!is.na(zvh), !is.na(ie)) %>%
      mutate(ie_class = cut2(ie, cuts = c(0, 0.25, 0.5, 0.75, 1))) %>%
    group_by(zvh) %>%
    summarise(
      n = n(),
      ie_mean = mean(ie), 
      ie_median = median(ie),
      ie_min = min(ie), 
      ie_max = max(ie), 
      ie_high = sum(ie_class == "[0.00,0.25)") / n,
      ie_medhigh = sum(ie_class == "[0.25,0.50)") / n,
      ie_medlow = sum(ie_class == "[0.50,0.75)") / n, 
      ie_low = sum(ie_class == "[0.75,1.00]") / n
    ) %>%
    mutate(
      file_path = basename(file_path)
    )
  zvh_summary
}

paths <- dir("../data/integrity_time", pattern = "\\.tif", full.names = TRUE)
zvh_ie <- ldply(paths, ieReadSum)

# year month date
zvh_ie$year <- as.numeric(substring(zvh_ie$file_path, 20, 23))
zvh_ie$month <- extract_numeric(substring(zvh_ie$file_path, 25, 26))
zvh_ie$date <- ymd(paste(zvh_ie$year, zvh_ie$month, "01", sep = "-"))

zvh_ie_f <- zvh_ie %>%
  filter(n > 1000) %>%
  group_by(zvh) %>%
  mutate(max_dif = as.numeric((max(ie_median) - min(ie_median)) > 0.06))

zvh_ie_f %>%
  group_by(zvh) %>%
  summarise(max_dif = max(ie_median) - min(ie_median)) %>%
  arrange(-max_dif)

zvh_ie_f$max_dif[zvh_ie_f$zvh %in% c(26, 27)] <- 2

zvh_ie_f$zvh_text <- "nada"
zvh_ie_f$zvh_text[zvh_ie_f$zvh == 26] <- "Subtropical thorn woodland"
zvh_ie_f$zvh_text[zvh_ie_f$zvh == 27] <- "Subtropical dry forest"
zvh_ie_f$zvh_text[zvh_ie_f$zvh == 21] <- "Warm temperate moist forest"
zvh_ie_f$zvh_text[zvh_ie_f$zvh == 8] <- "Boreal moist forest"

zvh_labels <- filter(zvh_ie_f, zvh %in% c(26, 21, 8, 27), 
    year == 2004, month == 2)
zvh_labels$ie_median[zvh_labels$zvh == 26] <- 0.45
zvh_labels$ie_median[zvh_labels$zvh == 27] <- 0.58
zvh_labels$ie_median[zvh_labels$zvh == 21] <- 0.91
zvh_labels$ie_median[zvh_labels$zvh == 8] <- 0.78


ggplot(zvh_ie_f, aes(x = date, y = ie_median, group = factor(zvh))) +
  geom_line(show_guide = FALSE, aes(color = factor(max_dif), 
    size = factor(max_dif))) +
  scale_color_manual(values = c("darkgray", "red", "blue")) +
  scale_size_manual(values = c(0.5, 1, 1)) +
  geom_text(data = zvh_labels, 
    aes(x = date, y = ie_median, label = zvh_text), vjust = 0, hjust = 0, 
    fontface = "bold") +
  labs(title = "Median of Integrity per Holdrige life zone", x = "time", y = "")


ggplot(filter(zvh_ie, n > 30), aes(x = date, y = ie_median, color = factor(zvh))) +
  geom_line()
  
ggplot(filter(zvh_ie, n > 30), aes(x = reorder(factor(zvh), ie_mean), y = ie_mean)) +
  geom_boxplot()
  
ggplot(filter(zvh_ie, n > 30), aes(x = year, y = ie_median)) +
  geom_boxplot()
  
library(googleVis)

zvh_ie$date_D <- as.Date(zvh_ie$date)
zvh_ie$zvh_cat<- factor(zvh_ie$zvh)
zvh_ie_f <- filter(zvh_ie, n > 30, month == 1) %>%
  mutate(
    zvh_cat = as.character(zvh)
  ) %>%
  select(zvh_cat, date_D, ie_min, year, ie_mean, n)

prueba <- gvisMotionChart(data = zvh_ie_f, idvar = "zvh_cat", timevar = "date_D", 
  xvar = "ie_min", yvar = "ie_mean", sizevar = "n")
  #date.format = "%Y-%m-%d")
plot(prueba)




#### ahora intentamos las zonas con cambios (cherry picking)
positive <- readOGR("../data/trends", "positive_ie_trend")
positive@proj4string
negative <- readOGR("../data/trends", "negative_ie_trend")
negative@proj4string

ieReadSumExtract <- function(file_path, trend){
  raster_anio <- raster(file_path)
  raster_stack <- stack(raster_anio, raster_zvh_crop)
  raster_shape <- raster::extract(raster_stack, trend)
  zvh_anio <- as.data.frame(raster_stack)
  colnames(zvh_anio) <- c("ie", "zvh")
  zvh_summary <- zvh_anio %>%
      filter(!is.na(zvh), !is.na(ie)) %>%
      mutate(ie_class = cut2(ie, cuts = c(0, 0.25, 0.5, 0.75, 1))) %>%
    group_by(zvh) %>%
    summarise(
      n = n(),
      ie_mean = mean(ie), 
      ie_median = median(ie),
      ie_min = min(ie), 
      ie_max = max(ie), 
      ie_high = sum(ie_class == "[0.00,0.25)") / n,
      ie_medhigh = sum(ie_class == "[0.25,0.50)") / n,
      ie_medlow = sum(ie_class == "[0.50,0.75)") / n, 
      ie_low = sum(ie_class == "[0.75,1.00]") / n
    ) %>%
    mutate(
      file_path = basename(file_path)
    )
  zvh_summary
}

zvh_ie_positive <- ldply(paths, ieReadSumExtract)
