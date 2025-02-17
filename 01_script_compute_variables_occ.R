#Script written by JP Vacher 5 March 2024 # 
#Updated 2 August 2024 #

#----------------------------------------#

#1. CHARGE LIBRARIES ####
library(sf) #to open and compute GIS data
library(tidyverse) #handle data tibbles

#2. CHARGE DATA ####
##2.1 monitoring grid ####
carres = st_read("couches_SIG_de_travail/Carres_200x200_SVF.shp")
carres

##2.2 occurrence points ####
###2.2.1 historical data ####
bombina1 = st_read("couches_SIG_de_travail/data_Bbombina_historiques.shp")
bombina1$ID_point = c(1:nrow(bombina1))
bombina1 = bombina1 %>%
  select(ID_point)
bombina1

###2.2.2 data collected in 2024 ####
bombina2 = st_read("couches_SIG_de_travail/data_Bbombina_JPV_2024.shp")
bombina2$ID_point = c(496:(495+nrow(bombina2)))
bombina2 = bombina2 %>%
  select(ID_point)
bombina2

###2.2.3 merge tables to get all occurrence points ####
bombina = bombina1 %>%
  bind_rows(bombina2)
bombina

#Check variables table in annex to get information and sources of the raw spatial data
##2.3 Forests: BD Foret V2 millesime 2014####
foret = st_read("couches_SIG_de_travail/boisements_moselle_BDFORETV2_2014.shp")

##2.4 water bodies: BD TOPO ####
planEau = st_read("couches_SIG_de_travail/plan_eau_moselle_BDTOPO2023.shp")

##2.5 rivers: BD TOPO ####
coursEau = st_read("couches_SIG_de_travail/reseau_hydro_albe_sarre.shp")

##2.6 meadows : RPG 2022 ####
#code groupe culture 16 ‘fourrage’
#code groupe culture 18 ‘prairies permanentes’
#code groupe culture 19 ‘prairies temporaires’
prairie = st_read("couches_SIG_de_travail/prairies_RPG2022.shp")


##2.7 crops: RPG 2022 ####
#code groupe culture 1, 2, 3, 4, 5, 6, 7, 8
cultures = st_read("couches_SIG_de_travail/cultures_RPG2022.shp")

#3. PREPARE CENTROIDS ####

##3.1 Centroids of cells of the grid ####
centroids = st_centroid(carres)

###let's check how it looks like:
ggplot() +
  geom_sf(data = carres) +
  geom_sf(data = centroids, size = .3)


#4. INTERSECT VARIABLES WITH CELLS AND COMPUTE VARIABLES ####

##4.1 DISTANCES ####
###4.1.1. Distance to nearest forest ####
#compute matrix of distances:

foret2 = st_geometry(obj = foret) %>%
  st_cast(to = "POLYGON") #ajust geometry from MULTIPOLYGON to POLYGON
distance.foret = st_geometry(obj = foret2) %>%
  st_cast(to = "POLYGON") %>%
  st_distance(y = carres) #compute distance

#select the minimum distance for each point:
distance.foret = tibble(ID_carre = centroids$ID_carre,
                        dist_foret = apply(distance.foret, 2, min))
distance.foret$dist_foret = round(distance.foret$dist_foret, 0)
distance.foret

###4.1.2. Distance to nearest river ####
#compute matrix of distances:

distance.coursEau = st_geometry(obj = coursEau) %>%
  st_distance(y = carres)

#select the minimum distance for each point:
distance.coursEau = tibble(ID_carre = centroids$ID_carre,
                        dist_coursEau = apply(distance.coursEau, 2, min))
distance.coursEau$dist_coursEau = round(distance.coursEau$dist_coursEau,0)
distance.coursEau

###4.1.3. Distance to the nearest water body ####
#compute matrix of distances:

distance.planEau = st_geometry(obj = planEau) %>%
  st_distance(y = carres)

#select the minimum distance for each point:
distance.planEau = tibble(ID_carre = centroids$ID_carre,
                           dist_planEau = apply(distance.planEau, 2, min))
distance.planEau$dist_planEau = round(distance.planEau$dist_planEau,0)
distance.planEau


###4.1.4 Distance to nearest occurrence point ####
#Select data before 2024 to take into account possible colonisation in 2024
#compute matrix of distances :

distance.bombina = st_geometry(obj = bombina1) %>%
  st_distance(y = carres)

#select the minimum distance for each point:
distance.bombina = tibble(ID_carre = centroids$ID_carre,
                          dist_bombina = apply(distance.bombina, 2, min))
distance.bombina$dist_bombina = round(distance.bombina$dist_bombina,0)
distance.bombina


##4.2. AREAS WITHIN CELLS ####
###4.2.1 forest area ####
intersect.foret = st_intersection(carres, foret2)
intersect.foret$area = st_area(intersect.foret)

foret.carres = tibble(intersect.foret) %>%
  group_by(ID_carre) %>%
  summarise(foret_surf = sum(area)) %>%
  mutate(foret_surf = as.double(foret_surf))


###beware because there are some cells without forests.
###therefore, we should add these cells in the table:
##first, extract cells with 0:
carres.zero = tibble(ID_carre = carres$ID_carre) %>%
  filter(!ID_carre %in% foret.carres$ID_carre) %>%
  mutate(foret_surf = 0)

###then append this table to the rest of the cells:
foret.carres = foret.carres %>%
  bind_rows(carres.zero)

foret.carres #this table should have 120 rows


###4.2.2 meadows area ####
intersect.prairie = st_intersection(carres, prairie)
intersect.prairie$area = st_area(intersect.prairie)

prairie.carres = tibble(intersect.prairie) %>%
  group_by(ID_carre) %>%
  summarise(prairie_surf = sum(area)) %>%
  mutate(prairie_surf = as.double(prairie_surf))


###beware because there are some cells without forests.
###therefore, we should add these cells in the table:
##first, extract cells with 0:
carres.zero = tibble(ID_carre = carres$ID_carre) %>%
  filter(!ID_carre %in% prairie.carres$ID_carre) %>%
  mutate(prairie_surf = 0)

###then append this table to the rest of the cells:
prairie.carres = prairie.carres %>%
  bind_rows(carres.zero)

prairie.carres #must contain 120 rows


###4.2.3 crops area ####
intersect.cultures = st_intersection(carres, cultures)
intersect.cultures$area = st_area(intersect.cultures)

cultures.carres = tibble(intersect.cultures) %>%
  group_by(ID_carre) %>%
  summarise(cultures_surf = sum(area)) %>%
  mutate(cultures_surf = as.double(cultures_surf))


###beware because there are some cells without forests.
###therefore, we should add these cells in the table:
##first, extract cells with 0:
carres.zero = tibble(ID_carre = carres$ID_carre) %>%
  filter(!ID_carre %in% cultures.carres$ID_carre) %>%
  mutate(cultures_surf = 0)

###then append this table to the rest of the cells:
cultures.carres = cultures.carres %>%
  bind_rows(carres.zero)

cultures.carres #must contain 120 rows


###4.2.4 water bodies area ####
intersect.planEau = st_intersection(carres, planEau)
intersect.planEau$area = st_area(intersect.planEau)

planEau.carres = tibble(intersect.planEau) %>%
  group_by(ID_carre) %>%
  summarise(planEau_surf = sum(area)) %>%
  mutate(planEau_surf = as.double(planEau_surf))

###beware because there are some cells without forests.
###therefore, we should add these cells in the table:
##first, extract cells with 0:
carres.zero = tibble(ID_carre = carres$ID_carre) %>%
  filter(!ID_carre %in% planEau.carres$ID_carre) %>%
  mutate(planEau_surf = 0)

###then append this table to the rest of the cells:
planEau.carres = planEau.carres %>%
  bind_rows(carres.zero)

planEau.carres #must contain 120 rows

###4.2.5 Number of small water bodies within cells and occurrence####
orniere = tibble(carres) %>%
  select(ID_carre, pres_2022, pres_2024, orniere) %>%
  rename(nb_planEau = orniere)

#5. OPEN DETECTION DATA TABLE AND ADD NUMBER OF WATER BODIES ####
det = read_tsv("table_var_detection.txt") %>% #open table with detection variables
  left_join(orniere) %>%
  select(ID_carre, temp1, temp2, temp3, temp4, nb_planEau) %>%
  rename(mare1 = nb_planEau) %>%
  mutate(mare2 = mare1, mare3 = mare1, mare4 = mare1) #replicate column as the values are the same at each survey occasion
  
det

#6. PREPARE WHOLE DATA TABLE ####

##6.1 Merge all data table by column ####
#This table will contain occupancy data (0/1), detection variables, and occupancy variables.

whole.data = det %>%
  left_join(distance.bombina) %>%
  left_join(distance.foret) %>%
  left_join(distance.coursEau) %>%
  left_join(distance.planEau) %>%
  left_join(foret.carres) %>%
  left_join(prairie.carres) %>%
  left_join(cultures.carres) %>%
  left_join(planEau.carres) %>%
  left_join(orniere) %>%
  arrange(ID_carre) %>%
  relocate(c(18,19), .after = 1)

whole.data #check how it look like. It should have 120 rows, and the first column should be Id_carre

##6.2 Save data table ####
write.table(data.frame(whole.data), "table_variables_SVF_2024.txt", sep = "\t", row.names = F)

#--------------END--------------#
