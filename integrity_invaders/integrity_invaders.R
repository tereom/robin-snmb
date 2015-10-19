# En este script buscamos relacionar la presencia de especies domésticas
# (vacas, borregos, perros, caballos) con integridad ecológica (ROBIN)

library(Hmisc)
library(rgdal)
library(sp)
library(raster)
library(dplyr)

# conexión a la base de datos (snmb)
PASS_SNMB = Sys.getenv("PASS_SNMB")
base_input <- src_postgres(dbname = "snmb", host = "dbms", user = "snmb", 
  password = PASS_SNMB)

# Para obtener el valor de integridad que corresponde a cada conglomerado 
# INFyS usamos el mapa de integridad y extraemos el valor correspondiente a
# las coordenadas del "Centro" del conglomerado (ingresadas en campo)

# encontrar coordenadas de sitios
conglomerado <- tbl(base_input, "conglomerado_muestra") %>%
  collect() %>%
  filter(id != 638 & id != 632) %>% # quitamos cgls repetidos
  # filter(nombre == "110296" | nombre == "158970" )
  select(conglomerado_muestra_id = id, nombre) 

sitio <- collect(tbl(base_input, "sitio_muestra")) %>%
  filter(sitio_numero == "Centro") %>%
  inner_join(conglomerado, by = "conglomerado_muestra_id") %>%
  mutate(
    lat = lat_grado + lat_min/60 + lat_seg/3600, 
    lon = ifelse(lon_grado > 0,  lon_grado + lon_min/60 + lon_seg/3600, 
      -(lon_grado - lon_min/60 - lon_seg/3600)), 
    lon = -lon, 
    cgl = as.numeric(nombre)
  ) %>%
  dplyr::select(cgl, conglomerado_muestra_id, lon, lat) 


# copiamos la proyección de este shape
malla_real <- readOGR("/Volumes/ARCHIVOS_C/Mac_Pro/SNMB/datos/malla_real", 
  "infys_2cWGS")
malla_real@proj4string

# y copiamos la proyección del raster de integridad para reproyectar a lcc
raster_ei <- raster("/Volumes/ARCHIVOS_MAC/Dropbox/Datos Redes Bayesianas/EI_maps/Stage_3/Final_net_Scores.tif")

sitio_shape <- as.data.frame(sitio)
coordinates(sitio_shape) <- ~ lon + lat
sitio_shape@proj4string <- malla_real@proj4string
sitio_lcc <- spTransform(sitio_shape, projection(raster_ei))

plot(raster_ei)
points(sitio_lcc)

sitio$ei <- raster::extract(raster_ei, sitio_lcc)
sitio$id <- as.character(sitio$cgl)
sitio$x_s <- sitio_lcc@coords[, 1]
sitio$y_s <- sitio_lcc@coords[, 2]

# Ahora usamos la siguiente función para buscar especies invasoras en los 
# conglomerados

### IMPORTANTE: hay un filtro ad-hoc en la tabla de conglomerado
especieInv <- function(noms){
  # noms: string que indica los nombres a agrupar separados por pipes, ej.
  #   "Bos|taurus|vaca"
  
  conglomerado <- tbl(base_input, "conglomerado_muestra") %>%
    collect() %>%
    filter(id != 638 & id != 632) %>% # quitamos cgls repetidos
    select(conglomerado_muestra_id = id, nombre, estado, municipio, 
      uso = uso_suelo_tipo)

  sitio <- tbl(base_input, "sitio_muestra") %>%
    collect() %>%
    select(conglomerado_muestra_id, sitio_muestra_id = id) %>%
    inner_join(conglomerado, by = "conglomerado_muestra_id") %>%
    select(sitio_muestra_id, conglomerado_muestra_id, nombre)
  
  tr_ei <- tbl(base_input, "transecto_especies_invasoras_muestra") %>%
    collect() %>%
    select(transecto_especies_invasoras_id = id, conglomerado_muestra_id) %>%
    left_join(conglomerado, by = "conglomerado_muestra_id")
  
  ei <- tbl(base_input, "especie_invasora") %>%
    collect() %>%
    mutate(
      ei_esp = grepl(noms, nombre_comun, ignore.case = TRUE) | 
        grepl(noms, nombre_cientifico, ignore.case = TRUE)
      ) %>%
    select(transecto_especies_invasoras_id, ei_esp) %>%
    left_join(tr_ei, by = "transecto_especies_invasoras_id") %>%
    group_by(nombre) %>%
    summarise(
      ei_esp = sum(ei_esp, na.rm = TRUE)
      ) %>%
    select(nombre, ei_esp)
  
  ei_ex <- tbl(base_input, "especie_invasora_extra") %>%
    collect() %>%
    mutate(
      ei_ex_esp = grepl(noms, nombre_comun, ignore.case = TRUE) | 
        grepl(noms, nombre_cientifico, ignore.case = TRUE)
      ) %>%
    select(conglomerado_muestra_id, ei_ex_esp) %>%
    right_join(conglomerado, by = "conglomerado_muestra_id") %>%
    group_by(nombre) %>%
    summarise(
      ei_ex_esp = sum(ei_ex_esp, na.rm = TRUE)
      )
  
  naZero <- function(x){
    ifelse(is.na(x), 0, (x > 0)*1)
  }
  desagregado <- conglomerado %>%
    left_join(ei) %>%
    left_join(ei_ex) %>%
    mutate_each(funs(naZero), contains("esp"))
  presencia <- desagregado %>%
    mutate(pres = (ei_esp + ei_ex_esp) > 0) %>%
    select(id = nombre, pres)
  list(desagregado = desagregado, presencia = presencia)
}

# extraemos los conglomerados con presencia de invasoras
especies_inv <- especieInv("Arundo donax|Bromus madritensis|Cyperus papyrus|Eichhornia crassipes|Hedera helix|Melinis minutiflora|Melinis repens|  Oeceoclades maculata|Pennisetum clandestinum|Rottboellia cochinchinensis|Salsola tragus")
especies_inv$desagregado
sum(especies_inv$presencia$pres)

# leemos la base con conglomerado e integridad producida por Julián
# vamos a favorecer esta base sobre la que creamos arriba, esto porque las
# coordenadas ingresadas manualmente tienen errores y Julián usó una malla
# teórica

ie <- read.csv("../data/congs_snmb_ie.csv")
head(ie)
head(sitio)
ggplot(ie, aes(x = ie)) + geom_histogram()

# creamos una variable categórica de integridad
ie_pastos <- ie %>%
  mutate(
    id = as.character(Cgl)
    ) %>%
  right_join(especies_inv$presencia) %>%
  left_join(select(sitio, id, ei, x_s, y_s)) %>%
  mutate(
    # si ie (base Julián) es NA tomamos el valor ei de sitio
    ie_join = ifelse(is.na(ie), ei, ie),
    ie_class = cut2(ie_join, cuts = c(0, 0.5, 0.75, 1)), 
    x = ifelse(is.na(x), x_s, x),
    y = ifelse(is.na(y), y_s, y)
    )

shape_invasoras <- select(ie_pastos, x, y, pres)
coordinates(shape_invasoras) <- ~ x + y
shape_invasoras@proj4string <- sitio_lcc@proj4string

writeOGR(shape_invasoras, "./shapes_invasoras", "invasoras", 
  driver = "ESRI Shapefile", verbose = FALSE, overwrite_layer = TRUE)


# vemos la proporción de inasoras en cada clase de integridad
ie_pastos %>%
  group_by(ie_class) %>%
  summarise(
    n = n(),
    num_invasoras = sum(pres, na.rm = TRUE), 
    prop_invasoras = round(mean(pres, na.rm = TRUE) * 100, 1), 
    se = round(sd(pres, na.rm = TRUE) / sqrt(n) * 100, 1)
  )
