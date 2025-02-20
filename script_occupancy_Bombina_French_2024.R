#MODELE D'OCCUPATION MULTI-ANNÉES
#PROJET SONNEUR À VENTRE DE FEU EN MOSELLE. ANNÉES 2022 ET 2024.
# Script rédigé par Jean-Pierre VACHER, avec la contribution de Vincent CLEMENT ###
#26.9.2024
#mis à jour 16.1.2025

#------------------------------------------------------------------------------------#

#1. CHARGER LES BIBLIOTHÈQUES ####
library(tidyverse) #pour manipuler les données et les tables
library(unmarked) #pour les modèles d'occupation
library(AICcmodavg) #pour l'estimation des modèles (model averaging)
library(MuMIn) #pour l'estimation des modèles
library(sjPlot) #pour la fonction tab_df
library(gridExtra) #pour la fonction grid.arrange qui permet de placer plusieurs figure dans le même panneau
library(sf) #pour ouvrir et manipuler des données spatiales, ce sera utile pour estimer l'autocorrelation spatiale
library(spdep) #pour la fonction dnearneigh (autocorrelation spatiale)

#2. LIRE LA TABLE DE DONNÉES ####
data <- read_tsv("dataset_for_model.txt") #lire la table
data #on regarde à quoi ça ressemble
summary(data) #on vérifie les valeurs des données
#Note : les données ont été mises à l'échelle (scaled).

#3. PREPARATION DES TABLES DE DONNÉES, MATRICES, ET LISTES ####
##3.1 MATRICE POUR LES DONNÉES D'OCCURRENCE ####
datay <- data %>%
  select(occ1, occ2, occ3, occ4) %>% #selection des colonnes d'occurrence
  as.matrix() #transformer en matrice
head(datay) #on regarde à quoi ça ressemble

##3.2 LISTE POUR LES VARIABLES DE DÉTECTION ####
obscov <- data %>%
  select(temp1, temp2, temp3, temp4, #selection des colonnes pour les variables de détection
         mare1, mare2, mare3, mare4) %>%
  as.matrix() #transformer en matrice
head(obscov) #on regarde à quoi ça ressemble

###3.2.1 liste des variables de détection pour 2022 ####
obscov2022 <- list(temp = obscov[,1:2], mares = obscov[, 5:6])

###3.2.2 liste des variables de détection pour 2024 ####
obscov2024 <- list(temp = obscov[,3:4], mares = obscov[, 7:8])

###3.2.3 liste des variables de détection pour 2022 et 2024 ####
obscov <- list(temp = obscov[,1:4], mares = obscov[, 5:8])

##3.3 TABLE POUR LES VARIABLES D'OCCURRENCE ####
#Note : les variables sont mises à l'échelle (scaled)
#On va procéder à une transformation en log pour éviter les écarts importants entre les valeurs
covarsite <- data %>%
  select(dist_bombina, dist_coursEau, dist_planEau, #selection des variables d'occurrence
         prairie_surf, planEau_surf, nb_planEau) %>%
  mutate_at(c(1:6), list(~(c(.) + 2))) %>% #on ajoute 2 pour ne pas avoir de valeurs négatives
  mutate_at(c(1:6), list(~log(c(.)))) %>% #transformation en log pour éviter des écarts de distribution entre les données
    data.frame() #transformer en data frame

##3.4 MATRICE DES ANNÉES ####
year <- matrix(data = c("2022", "2024"),nrow = nrow(datay), ncol = 2, byrow = T) #on créé une matrice des années


#4. MODÈLE D'OCCUPATION MULTI-SAISON 2022 - 2024 ####

##4.1 On créé une table ave toutes les données rassemblées : ####
umf.multi <- unmarkedMultFrame(datay, siteCovs = covarsite, obsCovs = obscov,
                               yearlySiteCovs = list(Year = year),  numPrimary = 2) 
#numPrimary = nombre d'années. Dans le cas présent, il y a 2 années de suivi dans 120 unités d'échantillonnage
summary(umf.multi)

#psi = occupation initiale
#gamma = colonisation
#epsilon = extinction
#p = détection


##4.2 PREMIÈRE ÉTAPE : SELECTION DU MEILLEUR TERME POUR LA DÉTECTION ####
###4.2.1 liste des modèles ####
####Note : psi, gamma et epsilon sont gardés constants lors de cette phase.
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

###4.2.2 Selection des modèles avec AICc ####
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

#Le meilleur modèle est celui qui retient la variable mares (= nombre de mares par unité d'échantillonnage).

confint(p.mares, type="det")
#pas de chevauchement avec 0
#on garde mares comme variable associée à la détection

####qualité de l'ajustement (goodness of fit, gof) pour le modèle p.mares : ####
system.time(mb.boot <- AICcmodavg::mb.gof.test(p.mares, nsim = 1000)) #1 min sur un MacBook Air 2020 M1 8Go
print(mb.boot, digit.vals = 4, digits.chisq = 4) # c-hat = 1.89  //   chi-s = 7.005, p-value = 0.221
#c-hat montre un peu de surdispertion
#gof semble ok ! On continue.


###4.2.3 g1 : Graphe de la détection en fonction de mares ####

#attributs de la mise à l'échelle de mares
#attr(,"scaled:center")
#[1] 1.708333
#attr(,"scaled:scale")
#[1] 1.349733

mares2 <- round(obscov$mares*1.349733 + 1.708333, 0) #on revient aux valeurs d'origine (backtransform)
nd.mares <- data.frame("mares" = (seq(min(mares2), max(mares2), length.out = 100))) #on créé une table avec 100 valeurs comprises entre le min et le max du nombre de mares
pred.p.mares <- predict(p.mares, type = 'det', newdata = nd.mares) #on créé une table de 100 valeurs de predict issues du modèle
predictions.mares <- cbind(pred.p.mares, nd.mares) #on fusionne les deux tables (valeurs et predict)
summary(predictions.mares) #on regarde si les valeurs semblent correctes

#on génère le graphe avec ggplot :
g1 = ggplot(predictions.mares, aes(x = mares, y = Predicted)) +
  geom_line() +
  geom_ribbon(alpha = .3, aes(ymin = lower, ymax = upper)) +
  ylim(0, 1) +
  theme_bw() +
  labs(x="Nombre de mares par maille", y = "Probabilité de détection") +
  theme(axis.text = element_text(size = 14), 
        axis.title = element_text(size = 18))
  
g1 #on regarde le résultat
#Cela suggère que plus le nombre de pièces d'eau augmente dans l'unité d'échantillon, plus la détectabilité augmente

#on sauvegarde la figure au format jpg :
jpeg("figure1.jpg", width = 17, height = 17, res = 300, units = "cm")
g1
dev.off()



##4.3 SECONDE ÉTAPE : SÉLECTION DU MEILLEUR MODÈLE D'OCCUPATION MULTI-ANNÉES ####

###4.3.1 Selection du meilleur modèle ####

####on construit un modèle complet, puis on estime l'AICc avec dredge :####
psi.total  <- colext(psiformula = ~ dist_bombina + dist_coursEau + dist_planEau +
                       prairie_surf + planEau_surf + nb_planEau, gammaformula = ~ 1, 
                               epsilonformula = ~ 1, pformula = ~ mares, data = umf.multi)

modelList <- dredge(psi.total, rank = "AICc") #la fonction dredge permet d'évaluer l'AICc pour chaque combinaison de variable au sein des modèles
tab_df(modelList[modelList$delta<2,]) #on sélectionne les modèles retenus avec un delta d'AICc <2
tab_df(modelList)
#le modèle le plus parcimonieux est celui avec psi(dist_bombina, dist_planEau et prairie_surf) p(mares)

####on ajuste le modèle avec dist_bombina, dist_planEau, prairie_surf ####
psi.bombinadplanEauprairie  <- colext(psiformula= ~ dist_bombina + dist_planEau +
                       prairie_surf, gammaformula = ~ 1, 
                     epsilonformula = ~ 1, pformula = ~ mares, data = umf.multi)

summary(psi.bombinadplanEauprairie) #on regarde l'effet des covariables d'occupation
#On remarque que les covariables dist_planEau et prairie_surf n'ont pas un effet fort sur l'estimation de l'occupation

####qualité d'ajustement (gof) pour le modèle psi.bombinadplanEauprairie : ####
system.time(mb.boot <- AICcmodavg::mb.gof.test(psi.bombinadplanEauprairie, nsim = 1000)) #1 min sur un MacBook Air M1 2020 8Go 
print(mb.boot, digit.vals = 4, digits.chisq = 4) # c-hat = 1.46  // chi-sq = 5.68  p-value = 0.221
#c-hat avec une légère surdispertion
#gof pas mal !

####Test d'un modèle avec un seul terme pour psi : dist_bombina ####
####ce terme a été retenu car c'est celui qui inlfuence le plus l'occupation dans le modèle à 3 termes
psi.bombina  <- colext(psiformula= ~ dist_bombina, gammaformula = ~ 1, 
                                      epsilonformula = ~ 1, pformula = ~ mares, data = umf.multi)
summary(psi.bombina)
AICc(psi.bombina)
#312.05

###qualité d'ajustement (gof) pour le modèle  psi.bombina : ####
system.time(mb.boot <- AICcmodavg::mb.gof.test(psi.bombina, nsim = 1000)) #1 min on a MacBook Air M1 first generation
print(mb.boot, digit.vals = 4, digits.chisq = 4) # c-hat = 1.6  // chi-square = 6.0046  p-value = 0.164
#c-hat avec une légère surdispertion
#gof pas mal !

#le c-hat est légèrement plus élevé que le modèle avec 3 termes, mais il reste correct (<2)
#le gof semble correct
#Ce modèle parcimonieux avec un seul terme est retenu

##4.4 Values des paramètres estimatés ####

###4.4.1 colonisation (gamma) ####
col <- backTransform(linearComb(psi.bombina, 1, 'col'))# colonisation
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

###4.4.3 détection (p) ####
nd <- data.frame(mares = 0, temp = 0)
det <- round(predict(psi.bombina, type='det', newdata = nd, appendData = T),2)
det
#   Predicted   SE lower upper
#1        0.5 0.06  0.39  0.62



#### 4.4.4 Retrouver les valeurs de l'occupation pour chaque année ###
## @smoothed = estimations de notre échantillon "fini", ce qu'on utilisera dans le rapport/article
## @projected = estimations projetées à l'échelle de la population
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



##4.6 g2 : graphe de l'occupation selon la distance à la mare occupée la plus proche ####

#attributs de la mise à l'échelle de mares :
#attr(,"scaled:center")
#[1] 500.725
#attr(,"scaled:scale")
#[1] 526.7162

nd.bombina <- data.frame(dist_bombina = seq(min(covarsite$dist_bombina), max(covarsite$dist_bombina), length.out = 120)) #on créé une table de 120 données avec des valeurs de distance comprises entre le min et le max de distance
pred.p.bombina <- predict(psi.bombina, type='psi', newdata = nd.bombina) #on créé une table de 120 valeurs de predict issues du modèle retenu
predictions.bombina <- cbind(pred.p.bombina, nd.bombina) #on fusionne les deux tables
predictions.bombina$dist_bombina2 = round((exp(predictions.bombina$dist_bombina)-2)*526.7162+500.725, 2) #on revient aux valeurs d'origine (backtransform)
summary(predictions.bombina) #on vérifie si les valeurs semblent cohérentes

#on génère le graphe avec ggplot :
g2 = ggplot(predictions.bombina, aes(x = dist_bombina2, y = Predicted)) +
  geom_line() +
  geom_ribbon(alpha=.3, aes(ymin = lower, ymax = upper)) +
  ylim(0, 1) +
  xlim(0,500) +
  theme_bw() +
  theme(axis.text = element_text(size = 14), 
        axis.title = element_text(size = 18)) +
  labs(x="Distance à la mare occupée la plus proche (m)", y="Probabilité d'occupation")
g2
#Cela suggère que le plus loin on se trouve de la pièce d'eau occupée la plus proche, le moins de chances une unité d'échantillon a d'être occupée


#on sauvegarde la figure au format jpg :
jpeg("figure2.jpg", width = 17, height = 17, res = 300, units = "cm") #on peut mettre une résolution de 200 (ok pour Powerpoint ou sites internet), mais 300 est mieux pour l'impression
g2
dev.off()

##4.7 on sauvegarde les deux graphes dans le même panneau ####
jpeg("figure3.jpg", height = 17, width = 17*2, units = "cm", res = 300)
grid.arrange(g1, g2, ncol = 2, nrow = 1)
dev.off()

#5. AUTOCORRELATION SPATIALE ####

##5.1 on construit une table des residus avec les coordonnées des centroïdes de chacune des 120 unités d'échantillonnage ####

###D'abord, on ouvre la couche des unités d'échantillonnage (120 mailles de 200 x 200 m):
carres <- st_read("GIS_Bombina_Moselle/Carres_200x200_SVF.shp")
carres
#carres2 <- carres %>% st_transform(4326) #pour transformer en WGS84, mais pas nécessaire ici

###Ensuite on ajoute les résidus du modèle, les ID des mailles et les coordonnées des centroïdes :
residuals = tibble(res = rowMeans(residuals(psi.bombina))) %>% #on calcule la moyenne des résidus par ligne
  mutate(ID_carre = carres$ID_carre, #ID des mailles
         lon = st_coordinates(st_centroid(carres))[,1], #longitude des centroïdes des mailles
         lat = st_coordinates(st_centroid(carres))[,2]) %>% #latitude des centroïdes des mailles
  arrange(ID_carre) #on classe par ID de carré
residuals #on regarde à quoi ça ressemble

##5.2 paramètres pour calculer l'indice de Moran ####

dmin <- 0 #distance minimum en mètres
dmax <- 10000 #distance maximum en mètres
by <- 500 #fourchette de distance pour les classes de distance en mètres
d0  <-  dmin + by #on définit la première classe de distance à partir de la distance minimum
dists  <-  seq(from = d0, to = dmax, by = by) #vecteur de classes de distances
n  <-  length(dists) #nombre de classes de distance
coords <- as.matrix(residuals[c("lon", "lat")]) #coordonnées géographiques des résidus
u  <-  dmin #distance minimum
MI  <-  numeric(0) #on définit une colonne vecteur numérique où sera stocké MI (indice de Moran) 
MI.p  <-  numeric(0) #on définit une colonne vecteur numérique où sera stockée MI.p (p-value de l'indice de Moran)
Crit.p  <-  c("*", rep("ns", n-1)) #colonne qui indique la significativité

##5.3 Boucle pour calculer Moran's I, p-value et le critères de significativité pour chaque classe de distance ####
for(i in 1:n){
  nb  <-  dnearneigh(coords, u, u + by, longlat = TRUE) # identifie les voisins
  nb.w  <-  nb2listw(nb, style="B", zero.policy = TRUE) # attribue les poids
  mi <- moran.test(residuals$res, nb.w, randomisation = TRUE,  #calcul de l'indice de Moran
                   na.action = na.exclude, zero.policy = TRUE)
  MI[i]  <-  mi$estimate #valeur de l'indice de Moran
  MI.p[i]  <-  mi$p.value #p.value de l'indice de Moran
  ifelse(MI.p[i] < 0.05, Crit.p[i] <- "*", Crit.p[i] <- "ns") #critère de significativité
  u=u+by
}

##5.4 stocker les résultats dans une table ####
res  <-  as.data.frame(cbind(dists=dists, MI = round(MI,2), MI.p = round(MI.p,3), 
                             Crit.p.B = Crit.p))
res #on regarde à quoi ça ressemble
#aucune classe de distance ne montre un indice de Moran significatif.
#on peut en conclure qu'il y a un effet faible de l'autocorrélation spatiale sur l'ajustement du modèle psi.bombina

##5.5 sauvegarder la table des résultats #####
write.table(res, "I_Moran_Bombina.txt", sep="\t", row.names=F)


#-----FIN-----#
