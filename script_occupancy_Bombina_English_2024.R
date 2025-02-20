#MULTI-SEASON OCCUPANCY MODEL
#FIRE-BELLIED TOAD IN MOSELLE PROJECT. YEARS 2022 AND 2024.
# Script written by Jean-Pierre VACHER, with the contribution of Vincent CLEMENT ###
#26 September 2024
#updated 16 January 2025

#------------------------------------------------------------------------------------#

#1. CHARGE LIBRARIES ####
library(tidyverse) #to manipulate data and tibbles
library(unmarked) #to fit occupancy models
library(AICcmodavg) #to estimate models (model averaging)
library(MuMIn) #to estimate models
library(sjPlot) #for tab_df function
library(gridExtra) #for grid.arrange function, to graph several figures within the same panel
library(sf) #to open, handle and compute spatial data, this will be useful for spatial autocorrelation
library(spdep) #for dnearneigh function (spatial autocorrelation)

#2. OPEN DATA TABLE ####
data <- read_tsv("dataset_for_model.txt") #open table
data #check how it looks like
summary(data) #check values
#Note: data were previously scaled.

#3. PREPARE DATA TABLES, MATRIX, AND LISTS ####
##3.1 MATRIX FOR OCCUPANCY DATA ####
datay <- data %>%
  select(occ1, occ2, occ3, occ4) %>% #select occupancy columns
  as.matrix() #transform to matrix
head(datay) #check how it looks like

##3.2 LIST FOR DETECTION VARIABLES ####
obscov <- data %>%
  select(temp1, temp2, temp3, temp4,
         mare1, mare2, mare3, mare4) %>% #select columns for detection variables
  as.matrix() #transform to matrix
head(obscov) #check how it looks like

###3.2.1 list of detection variable in 2022 ####
obscov2022 <- list(temp = obscov[,1:2], mares = obscov[, 5:6])

###3.2.2 list of detection variable in 2024 ####
obscov2024 <- list(temp = obscov[,3:4], mares = obscov[, 7:8])

###3.2.3 list of detection variable in 2022 and 2024 ####
obscov <- list(temp = obscov[,1:4], mares = obscov[, 5:8])

##3.3 TABLE OF OCCUPANCY VARIABLES ####
#Note:variables were previously scaled
#we transform values into log to avoid large gaps between values
covarsite <- data %>%
  select(dist_bombina, dist_coursEau, dist_planEau, prairie_surf, planEau_surf, nb_planEau) %>% #select occupancy variables
  mutate_at(c(1:6), list(~(c(.) + 2))) %>% #add 2 to avoid negative values
  mutate_at(c(1:6), list(~log(c(.)))) %>% #transforma into log to avoid distribution gaps
    data.frame() #transform to data frame

##3.4 MATRIX OF YEARS ####
year <- matrix(data = c("2022", "2024"),nrow = nrow(datay), ncol = 2, byrow = T) #create a matrix of years


#4. MULTI-SEASON OCCUPANCY MODEL 2022 - 2024 ####

##4.1 generate a table with all data gathered: ####
umf.multi <- unmarkedMultFrame(datay, siteCovs = covarsite, obsCovs = obscov,
                               yearlySiteCovs = list(Year = year),  numPrimary = 2) 
#numPrimary = number of years. In this case, there are 2 years of monitoring in 120 sample cells
summary(umf.multi)

#psi =  initial occupancy
#gamma = colonization
#epsilon = extinction
#p = detection


##4.2 FIRST STEP: SELECTION OF THE BEST TERM FOR DETECTION ####
###4.2.1 list of models ####
####Note: psi, gamma, and epsilon are kept constant.
p.null <- colext(psiformula = ~ 1, gammaformula = ~ 1, epsilonformula = ~ 1, pformula = ~ 1, data = umf.multi)
p.temp <- colext(psiformula = ~ 1, gammaformula = ~ 1, epsilonformula = ~ 1, pformula = ~ temp, data = umf.multi)
p.mares <- colext(psiformula = ~ 1, gammaformula = ~ 1, epsilonformula = ~ 1, pformula = ~ mares, data = umf.multi)
p.year <- colext(psiformula = ~ 1, gammaformula = ~ 1, epsilonformula = ~ 1, pformula = ~ Year, data = umf.multi)
p.tempmares <- colext(psiformula = ~ 1, gammaformula = ~ 1, epsilonformula = ~ 1, pformula = ~ temp + mares, data = umf.multi)
p.tempyear <- colext(psiformula = ~ 1, gammaformula = ~ 1, epsilonformula = ~ 1, pformula = ~ temp + Year, data = umf.multi)
p.mareyear <- colext(psiformula = ~ 1, gammaformula = ~ 1, epsilonformula = ~ 1, pformula = ~ mares+ Year, data = umf.multi)
p.tempmareyear <- colext(psiformula = ~ 1, gammaformula = ~ 1, epsilonformula = ~ 1, pformula = ~ temp + mares + Year, data = umf.multi)

modlist1 <- list(p.null = p.null, p.temp = p.temp, 
               p.mares = p.mares, p.year = p.year,
               p.tempmares = p.tempmares, p.tempyear = p.tempyear,
               p.mareyear = p.mareyear, p.tempmareyear = p.tempmareyear)

###4.2.2 model selection with AICc ####
aictab(modlist1)
#.              K   AICc Delta_AICc AICcWt Cum.Wt      LL
#p.mares        5 357.13       0.00   0.36   0.36 -173.30
#p.tempmares    6 357.57       0.44   0.29   0.66 -172.41
#p.mareyear     6 358.23       1.10   0.21   0.87 -172.74
#p.tempmareyear 7 359.17       2.04   0.13   1.00 -172.09
#p.temp         5 376.31      19.18   0.00   1.00 -182.89
#p.null         4 376.98      19.86   0.00   1.00 -184.32
#p.year         5 377.34      20.21   0.00   1.00 -183.41
#p.tempyear     6 377.85      20.72   0.00   1.00 -182.55

#The best model is the one with the pond (‘mares’) variable (= number of water bodies per cell).

confint(p.mares, type="det")
#no overlapping with 0
#we keep ‘mares’ as a variable associated with detection

####goodness of fit (gof) for the model p.mares : ####
system.time(mb.boot <- AICcmodavg::mb.gof.test(p.mares, nsim = 1000)) #1 minute on a MacBook Air 2020 M1 8Go
print(mb.boot, digit.vals = 4, digits.chisq = 4) # c-hat = 1.89  //   chi-s = 7.005, p-value = 0.221
#c-hat shows minor overdispersal
#gof seems ok!


###4.2.3 g1: Graph detection according to ponds ####

#attributes of scaling of ‘mares’
#attr(,"scaled:center")
#[1] 1.708333
#attr(,"scaled:scale")
#[1] 1.349733

mares2 <- round(obscov$mares*1.349733 + 1.708333, 0) #backtransform
nd.mares <- data.frame("mares" = (seq(min(mares2), max(mares2), length.out = 100))) #create a data frame with 100 values that range from min to max of ‘mares’
pred.p.mares <- predict(p.mares, type = 'det', newdata = nd.mares) #create a data frame with 100 values of predicts from the model
predictions.mares <- cbind(pred.p.mares, nd.mares) #bind tables
summary(predictions.mares) #check the values

#now graph number of ponds over detection probability:
g1 = ggplot(predictions.mares, aes(x = mares, y = Predicted)) +
  geom_line() +
  geom_ribbon(alpha = .3, aes(ymin = lower, ymax = upper)) +
  ylim(0, 1) +
  theme_bw() +
  labs(x="Number of water bodies per cell", y = "Detection probability") +
  ggtitle("(A)") +
  theme(axis.text = element_text(size = 14), 
        axis.title = element_text(size = 18),
        title = element_text(size = 18))
  
g1
#this shows that the higher the number of water bodies within a given sample cell, the higher the detection probability

#save the graph as a jpg file:
jpeg("figure1.jpg", width = 17, height = 17, res = 300, units = "cm")
g1
dev.off()


##4.3 SECOND STEP: SELECT THE BEST MULTI-SEASON OCCUPANCY MODEL ####

###4.3.1 select the best model ####

####first we build a model with all variables, then estimate AICc with dredge :####
psi.total  <- colext(psiformula = ~ dist_bombina + dist_coursEau + dist_planEau +
                       prairie_surf + planEau_surf + nb_planEau, gammaformula = ~ 1, 
                               epsilonformula = ~ 1, pformula = ~ mares, data = umf.multi)

modelList <- dredge(psi.total, rank = "AICc")
tab_df(modelList[modelList$delta<2,]) #select models with delta AICc <2
tab_df(modelList)
#the most parcimoniousmodel is psi(dist_bombina, dist_planEau et prairie_surf)gamma(.)epsilon(.)p(mares)

####let's fit this model with psi(dist_bombina, dist_planEau, prairie_surf) and p(mares) ####
psi.bombinadplanEauprairie  <- colext(psiformula= ~ dist_bombina + dist_planEau +
                       prairie_surf, gammaformula = ~ 1, 
                     epsilonformula = ~ 1, pformula = ~ mares, data = umf.multi)

summary(psi.bombinadplanEauprairie) #check the effect of landscape variables on occupancy
#we notice that dist_planEau (distance to water bodies) and prairie_surf (area of meadows) do not have a strong effect on estimation of occupancy probability

####goodness of fit (gof) for model psi.bombinadplanEauprairie: ####
system.time(mb.boot <- AICcmodavg::mb.gof.test(psi.bombinadplanEauprairie, nsim = 1000)) #1 min an a MacBook Air M1 2020 8Go 
print(mb.boot, digit.vals = 4, digits.chisq = 4) # c-hat = 1.46  // chi-sq = 5.68  p-value = 0.221
#c-hat with slight overdispersal
#gof not so bad!

####Fit a model with only on term for psi: dist_bombina ####
####this term has been selected because it inlfuences the most occupancy in the model with 3 terms
psi.bombina  <- colext(psiformula= ~ dist_bombina, 
                       gammaformula = ~ 1,
                       epsilonformula = ~ 1, 
                       pformula = ~ mares, 
                       data = umf.multi)
summary(psi.bombina)
AICc(psi.bombina)
#312.05

###goodness of fit (gof) for model  psi.bombina : ####
system.time(mb.boot <- AICcmodavg::mb.gof.test(psi.bombina, nsim = 1000)) #1 min on a MacBook Air M1 first generation
print(mb.boot, digit.vals = 4, digits.chisq = 4) # c-hat = 1.6  // chi-square = 6.0046  p-value = 0.164
#c-hat with slight overdispersal (<2)
#gof not so bad!

#c-hat is slightly higher than for the model with 3 terms, but still acceptable (<2)
#gof is ok
#Finaly, we select this parsimonious model with only one term

##4.4 Values of estimated parameters ####

###4.4.1 colonization (gamma) ####
col <- backTransform(linearComb(psi.bombina, 1, 'col'))# colonization
col
#Estimate     SE LinComb (Intercept)
#   0.181 0.0622   -1.48           1
round(confint(col), 2)
#     0.025     0.975
#      0.09 0.34

###4.4.2 extinction (epsilon) ####
ext <- backTransform(linearComb(psi.bombina, 1,'ext')) # extinction
ext
# Estimate     SE LinComb (Intercept)
#  0.000834 0.0279   -7.09          1
round(confint(ext),2)
#.      0.025 0.975
#.          0     1

###4.4.3 detection (p) ####
nd <- data.frame(mares = 0, temp = 0)
det <- round(predict(psi.bombina, type='det', newdata = nd, appendData = T),2)
det
#   Predicted   SE lower upper
#1        0.5 0.06  0.39  0.62



#### 4.4.4 Occupancy values for each year ###
## @smoothed = estimations of the bounded sample, these are the ones we will use for the article
## @projected = projected estimations at the scale of the population
mb <- nonparboot(psi.bombina, B = 500)
sample_occupancy <- data.frame(year = c("2022","2024"),
                                 smoothed_occ = round(smoothed(mb)[2,],2),
                                 SE = round(mb@smoothed.mean.bsse[2,],2))
sample_occupancy
#  year smoothed_occ   SE
#1 2022         0.3 0.06
#2 2024         0.43 0.07

round(mean(sample_occupancy$smoothed_occ), 2)
#[1] 0.36



##4.6 g2: graphing occupancy according to distance to the nearest occupied water body ####

#attributes of the scaling of the dist_bombina variable:
#attr(,"scaled:center")
#[1] 500.725
#attr(,"scaled:scale")
#[1] 526.7162

nd.bombina <- data.frame(dist_bombina = seq(min(covarsite$dist_bombina), #create a data frame with 100 values that range from min to max of dist_bombina
                                            max(covarsite$dist_bombina), length.out = 120))
pred.p.bombina <- predict(psi.bombina, type='psi', newdata = nd.bombina)#create a data frame with 100 values of predicts from the model
predictions.bombina <- cbind(pred.p.bombina, nd.bombina) #bind tables
predictions.bombina$dist_bombina2 = round((exp(predictions.bombina$dist_bombina)-2)*526.7162+500.725, 2) #backtransform
summary(predictions.bombina)

#now graph distance to the nearest occupied water body over occupancy probability:
g2 = ggplot(predictions.bombina, aes(x = dist_bombina2, y = Predicted)) +
  geom_line() +
  geom_ribbon(alpha=.3, aes(ymin = lower, ymax = upper)) +
  ylim(0, 1) +
  xlim(0,500) +
  labs(x="Distance to the nearest occupied pond (m)", y="Occupancy probability") +
  ggtitle("(B)") +
  theme_bw() +
  theme(axis.text = element_text(size = 14), 
        axis.title = element_text(size = 18),
        title = element_text(size = 18))
  
g2
#This shows that the further the cell from the nearest occupied water body, the less chance has a cell to be occupied


#save as jpg file :
jpeg("figure2.jpg", width = 17, height = 17, res = 300, units = "cm")
g2
dev.off()

##4.7 save both graphs g1 and g2 in the same panel in a jpg file ####
jpeg("figure5.jpg", height = 17, width = 17*2, units = "cm", res = 300) #we can downgrade res to 200, but 300 is better for printing
grid.arrange(g1, g2, ncol = 2, nrow = 1)
dev.off()

#5. SPATIAL AUTOCORRELATION####

##5.1 build a table with residuals of centroids coordinates of each 120 cells of the sample grid ####

###first, open the layer of the sample grid (120 cells of 200 m x 200 m):
carres <- st_read("GIS_Bombina_Moselle/Carres_200x200_SVF.shp")
carres
#carres2 <- carres %>% st_transform(4326) #transform to WGS84, not necessary here

###then add residuals of the model, cell IDs, and centroids coordinates:
residuals = tibble(res = rowMeans(residuals(psi.bombina))) %>%
  mutate(ID_carre = carres$ID_carre,
         lon = st_coordinates(st_centroid(carres))[,1],
         lat = st_coordinates(st_centroid(carres))[,2]) %>%
  arrange(ID_carre)
residuals

##5.2 parameters to compute Moran's I ####

dmin <- 0 #minimum distance in meters
dmax <- 10000 #maximum distance in meters
by <- 500 #distance gap for distance classes in meters
d0  <-  dmin + by #first distance class is defined from the minimum distance
dists  <-  seq(from = d0, to = dmax, by = by) #vector of distance classes
n  <-  length(dists) #number of distance classes
coords <- as.matrix(residuals[c("lon", "lat")]) #geographic coordinates of residuals
u  <-  dmin #minimum distance
MI  <-  numeric(0) #defines a numeric vector columns where we will store MI (Moran'sI)
MI.p  <-  numeric(0) #defines a numeric vector columns where we will store MI.p (p-value of Moran's I)
Crit.p  <-  c("*", rep("ns", n-1)) #column with significance

##5.3 Loop to computeMoran's I, p-value, and significance criteria for each distance class ####
for(i in 1:n){
  nb  <-  dnearneigh(coords, u, u + by, longlat = TRUE) # identifies neighbors
  nb.w  <-  nb2listw(nb, style="B", zero.policy = TRUE) # attributes weights
  mi <- moran.test(residuals$res, nb.w, randomisation = TRUE, 
                   na.action = na.exclude, zero.policy = TRUE) #compute Moran's I
  MI[i]  <-  mi$estimate #Moran's I value
  MI.p[i]  <-  mi$p.value#Moran's I p.value
  ifelse(MI.p[i] < 0.05, Crit.p[i] <- "*", Crit.p[i] <- "ns") #significance criteria
  u=u+by
}

##5.4 stock results in a table####
res  <-  as.data.frame(cbind(dists=dists, MI = round(MI,2), 
                             MI.p = round(MI.p,3), Crit.p.B = Crit.p))
res #check how it looks like
#In this case:
#no distance class displays any significant Moran's I.
#we can conclude that there is a weak effect of spatial autocorrelation on fitting of model psi.bombina

##5.5 save table of Moran's I results #####
write.table(res, "I_Moran_Bombina.txt", sep="\t", row.names=F)


#----------------END----------------#
