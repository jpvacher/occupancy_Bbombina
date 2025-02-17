#-------------------------------------------------------#
#Script written by Jean-Pierre Vacher, 10 March 2024 #
#updated 12 September 2024 #
#-------------------------------------------------------#


#-------------------------------------------------------#
#1.  CHARGE PACKAGES ####
#-------------------------------------------------------#

library("tidyverse") #for tables and data manipulation
library("usdm") #for vif function
library("corrplot") #for correlation plots
library("FactoMineR") #for PCA
library("factoextra") #for PCA
library("car") #for Anova function
library("gridExtra") #to plot several graphics on the same panel
library("Hmisc") #for rcorr function



#-------------------------------------------------------#
#2. OPEN DATA TABLE ####
#-------------------------------------------------------#

occ.var = read_tsv("table_variables_SVF_2024.txt")
occ.var #check how it looks like
summary(occ.var) #check if values seem ok

#explore variable number of water bodies with histogram:
ggplot(occ.var, aes(x = nb_planEau)) +
  geom_histogram()

#there is one outlier with 20 water bodies.
#let's reduce to 5 water bodies, so 5 will be the maximum value for this variable:
occ.var.cor = occ.var %>%
  mutate(nb_planEau = ifelse(nb_planEau >5, 5, nb_planEau),
         mare1 = ifelse(mare1 >5, 5, mare1),
         mare2 = ifelse(mare2 >5, 5, mare2),
         mare3 = ifelse(mare3 >5, 5, mare3),
         mare4 = ifelse(mare4 >5, 5, mare4))

summary(occ.var.cor) 
#lets's see what it looks like now graphically:
g1 = ggplot(occ.var.cor, aes(x = nb_planEau)) +
  geom_histogram() +
  ggtitle("nbWaterbodies")
g1

#OK, no more outlier

#check the other variables:
g2 = ggplot(occ.var.cor, aes(x = dist_bombina)) +
  geom_histogram() +
  ggtitle("distBombina")

g3 = ggplot(occ.var.cor, aes(x = dist_foret)) +
  geom_histogram() +
  ggtitle("distForest")

g4 = ggplot(occ.var.cor, aes(x = dist_coursEau)) +
  geom_histogram() +
  ggtitle("distRivers")

g5 = ggplot(occ.var.cor, aes(x = dist_planEau)) +
  geom_histogram() +
  ggtitle("distWaterbodies")

g6 = ggplot(occ.var.cor, aes(x = foret_surf)) +
  geom_histogram() +
  ggtitle("surfForest")

g7 = ggplot(occ.var.cor, aes(x = prairie_surf)) +
  geom_histogram() +
  ggtitle("surfMeadow")

g8 = ggplot(occ.var.cor, aes(x = cultures_surf)) +
  geom_histogram() +
  ggtitle("surfCrops")

g9 = ggplot(occ.var.cor, aes(x = planEau_surf)) +
  geom_histogram() +
  ggtitle("surfWaterbodies")

g10 = ggplot(occ.var.cor, aes(x = nb_planEau)) +
  geom_histogram() +
  ggtitle("nbPonds")

grid.arrange(g1, g2, g3, g4, g5, g6, g7, g8, g9, g10)

#scale variables to get them all within the same range:
occ.var3 = occ.var.cor %>%
  mutate_at(c(4:20), list(~c(scale(.))))
occ.var3
summary(occ.var3)

#get the attributes of each variable, this will be usefull to backtransform
att.mare = scale(occ.var.cor$nb_planEau)
att.mare
#attr(,"scaled:center")
#[1] 1.708333
#attr(,"scaled:scale")
#[1] 1.349733

att.temp = scale(occ.var.cor$temp1)
att.temp
#attr(,"scaled:center")
#[1] 23.375
#attr(,"scaled:scale")
#[1] 4.15268

att.bombina = scale(occ.var.cor$dist_bombina)
att.bombina
#attr(,"scaled:center")
#[1] 500.725
#attr(,"scaled:scale")
#[1] 526.7162

att.prairie = scale(occ.var.cor$prairie_surf)
att.prairie
#attr(,"scaled:center")
#[1] 11038.48
#attr(,"scaled:scale")
#[1] 12707.66


#-------------------------------------------------------#
#3. EXPLORATION OF DATA ####
#-------------------------------------------------------#

colnames(occ.var3)

##3.1 Correlation with Variance Inflation Factor VIF ####
VIF = usdm::vif(data.frame(occ.var3[,c(12:20)]))
VIF$VIF = round(VIF$VIF, 2)
VIF
#.     Variables  VIF
#1  dist_bombina  1.40
#2    dist_foret  1.68
#3 dist_coursEau  1.54
#4  dist_planEau  1.95
#5    foret_surf 10.53 #VIF>5
#6  prairie_surf  5.66 #VIF>5
#7 cultures_surf  3.17
#8  planEau_surf  3.12
#9    nb_planEau  1.18

#forest (distance and area), distance and area to water bodies are correlated.
#let's discard them from the data set and perform VIF again:
VIF2 = usdm::vif(data.frame(occ.var3[,c(12:15, 17:20)]))
VIF2$VIF = round(VIF2$VIF, 2)
VIF2
#.     Variables  VIF
#1  ddist_bombina 1.39
#2    dist_foret 1.56
#3 dist_coursEau 1.47
#4  dist_planEau 1.78
#5  prairie_surf 1.72
#6 cultures_surf 1.38
#7  planEau_surf 1.85
#8    nb_planEau 1.17

#OK, all VIF values <5, let's keep this data set


##3.2 Correlation with Pearson's correlation  ####
###3.2.1 with all variables ####
cor = rcorr(as.matrix(occ.var3[,c(12:20)]), type = "pearson") #let's check the first matrix, correlated variables corrélées are the ones >0.7
corrplot(cor$r, method = "number",
         type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45, tl.cex = .5,
         number.cex = .5)

#forest and meadows area are slighty negatively correlated (-68%)
#forest area and distance to forest are slightly correlated (-59%)
# -> variables foret_surf can be discarded, it is in line with the VIF

##3.3 PCA visualization ####
res.pca <- PCA(occ.var3[,c(12:15, 17:20)], graph=F) # select variables, we discard forest surface
eig.val <- get_eigenvalue(res.pca) #get eigen values
eig.val #check eigen values
fviz_eig(res.pca, addlabels = T, ylim = c(0,40)) #screeplot
var <- get_pca_var(res.pca)
round(var$contrib,2) #table with the contribution of variable for each dimension
round(var$cor,2) #round to 2 digits

#dimensions 1 and 2 represent 49.65% of total variance

###Visualization of variables ####
fviz_pca_var(res.pca, label="var", col.var = "black", repel = T)


#area and distance to water bodies are inversely correlated on dimension 1.
#we can nonetheless keep these two variables because they bear different types of information

#------------------------------------------------------#
#4. SIMPLE LINEAR MODELS ON EACH LANDSCAPE VARIABLE ####
#------------------------------------------------------#

#first prepare the data table by adding occurrence data for both years of monitoring:
occ.var4 = occ.var3 %>%
  mutate(occurrence = pres_2022 + pres_2024,
         occurrence = ifelse(occurrence == 0, 0,
                             ifelse(occurrence == 1, 1, 1))) %>%
  relocate(occurrence, .after = pres_2024)
occ.var4


#4.1 Occurrence according to distance to nearest occupied pond ####
anova1 = lm(occurrence ~ dist_bombina, data = occ.var4)
summary(anova1)
#               Estimate Std. Error t value Pr(>|t|)    
#(Intercept)   0.34167    0.03761   9.084 2.94e-15 ***
#dist_bombina -0.24181    0.03777  -6.402 3.24e-09 ***

anova(anova1)
#              Df  Sum Sq Mean Sq F value Pr(>F)
#dist_foret     1  0.5945 0.59455  2.6577 0.1057
#Residuals  118 26.3971 0.22370 

#we keep the variable distance to the nearest occupied pond


#4.2 Occurrence according to distance to nearest forest ####
anova2 = lm(occurrence ~ dist_foret, data = occ.var4)
summary(anova2)
#             Estimate Std. Error t value Pr(>|t|)    
#(Intercept)   0.34167    0.04318   7.913 1.51e-12 ***
#dist_foret  -0.07068    0.04336  -1.630    0.106  

anova(anova2)
#            Df  Sum Sq Mean Sq F value Pr(>F)
#dist_foret   1  0.5945 0.59455  2.6577 0.1057
#Residuals  118 26.3971 0.22370  

#we discard distance to the nearest forest

#4.3 Occurrence according to distance to nearest river ####
anova3 = lm(occurrence ~ dist_coursEau, data = occ.var4)
plot(anova3)
summary(anova3)
#              Estimate Std. Error t value Pr(>|t|)    
#(Intercept)    0.34167    0.04292   7.961 1.18e-12 ***
#dist_coursEau -0.08747    0.04310  -2.030   0.0446 * 
  
anova(anova3)
#               Df  Sum Sq Mean Sq F value Pr(>F)
#dist_coursEau  1  0.9105 0.91047  4.1193 0.04465 *
#Residuals     118 26.0812 0.22103 

#we keep distance to the nearest river


#4.4 Occurrence according to distance to nearest water body ####
anova4 = lm(occurrence ~ dist_planEau, data = occ.var4)
summary(anova4)
#             Estimate Std. Error t value Pr(>|t|)    
#(Intercept)   0.34167    0.04174   8.185 3.61e-13 ***
#dist_planEau -0.13959    0.04192  -3.330  0.00116 ** 

anova(anova4)
#              Df  Sum Sq Mean Sq F value Pr(>F)
#dist_planEau   1  2.3186 2.31860  11.089 0.001159 **
#Residuals    118 24.6731 0.20909 

#we keep distance to the nearest water body


#4.5 Occurrence according to surface of meadows ####
anova5 = lm(occurrence ~ prairie_surf, data = occ.var4)
summary(anova5)
#             Estimate Std. Error t value Pr(>|t|)    
#(Intercept)   0.34167    0.04254   8.032  8.1e-13 ***
#prairie_surf -0.10737    0.04271  -2.514   0.0133 *  

anova(anova5)
#              Df  Sum Sq Mean Sq F value Pr(>F)
#prairie_surf   1  1.3719 1.37193  6.3189 0.01329 *
#Residuals    118 25.6197 0.21712  

#we keep meadows area


#4.6 Occurrence according to surface of crops ####
anova6 = lm(occurrence ~ cultures_surf, data = occ.var4)
summary(anova6)
#               Estimate Std. Error t value Pr(>|t|)    
#(Intercept)   0.34167    0.04334   7.884 1.76e-12 ***
#cultures_surf -0.05798    0.04352  -1.332    0.185   

anova(anova6)
#               Df  Sum Sq Mean Sq F value Pr(>F)
#cultures_surf   1  0.4001 0.40009  1.7754 0.1853
#Residuals     118 26.5916 0.22535 

#we discard crops area


#4.7 Occurrence according to surface of water bodies ####
anova7 = lm(occurrence ~ planEau_surf, data = occ.var4)
summary(anova7)
#             Estimate Std. Error t value Pr(>|t|)    
#(Intercept)   0.34167    0.04263   8.014 8.92e-13 ***
#planEau_surf  0.10264    0.04281   2.398   0.0181 * 

anova(anova7)
#              Df  Sum Sq Mean Sq F value Pr(>F)
#planEau_surf   1  1.2538 1.25376  5.7481 0.01808 *
#Residuals    118 25.7379 0.21812   

#we keep water bodies area


#4.8 Occurrence according to number of water bodies ####
anova8 = lm(occurrence ~ nb_planEau, data = occ.var4)
summary(anova8)
#                 Estimate Std. Error t value Pr(>|t|)    
#(Intercept)       0.34167    0.04127   8.279  2.2e-13 ***
#nb_planEau   0.15539    0.04144   3.749 0.000276 ***

anova(anova8)
#            Df  Sum Sq Mean Sq F value Pr(>F)
#nb_planEau   1  2.8734 2.87335  14.058 0.0002761 ***
#Residuals  118 24.1183 0.20439  

#we keep number of water bodies



#-------------------------------------------------------#
#5. FILTER THE DATASET ####
#-------------------------------------------------------#

#select the column that will be needed for the models:
occ.var.red = occ.var4 %>%
  select(ID_carre, temp1, temp2, temp3, temp4,
         mare1, mare2, mare3, mare4, 
         dist_bombina, dist_coursEau, dist_planEau, prairie_surf,
         planEau_surf, nb_planEau)


occ.var.red #check how it look like. Should have 120 rows.
summary(occ.var.red) #check the values

#log values (not necessary here):

#summary(occ.var.red)
#occ.var.red.log = occ.var.red %>%
#  mutate_at(c(4:7), funs((.)^2)) %>%
#  mutate_at(c(4:7), funs(log(.))) %>%
#  rename_at(c(4:7), function(x) paste0(x, "_log"))

#summary(occ.var.red.log)

#-------------------------------------#
#7. ADD OCCURRENCE DATA BY SESSION ####
#-------------------------------------#

##7.1 open occupancy table: ####
session = read_tsv("occ_sessions_SVF.txt")
session

##7.2 merge occupancy table with variables: ####
data.table = session %>%
  left_join(occ.var.red)
data.table #should have 120 rows
summary(data.table) #check the values

##7.3 save final data table that will be used for modelling: ####
write.table(data.frame(data.table), "dataset_for_model.txt", sep = "\t", row.names = F)


#-----------------END-----------------#