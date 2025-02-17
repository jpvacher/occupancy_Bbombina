# Occupancy modelling of introduced Bombina bombina Moselle - 2024.

#Jean-Pierre Vacher, Conservatoire d'espaces naturels de Lorraine, February 2025.


## ENGLISH

In this repo, you will find a script to fit a multi-year occupancy model on a two-year monitoring data set of Bombina bombina in Moselle (France), and to assess spatial autocorrelation.

###script_occupancy_Bombina_English_2024.R

Data collection is based on a standardized protocol on a spatial sample comprising of 120 200 x 200 cells randomly placed within and outside the known range of the species in the area of Albestroff (Moselle, France) at the beginning of the monitoring in 2022, to assess distribution shifts over time. 
Each cell was visited twice per year, and monitoring takes place every two years. The first year of monitoring was 2022, and it was carried out in 2024, making it a two-year data set.
We used the package unmarked to fit a multi-year occupancy model to this data set to evaluate detection, multi-year occupancy, colonisation, and extinction.
We provide a table of already-transformed data, that are ready to use with the script. We also provide a spatial layer of samples units that is useful to assess spatial autocorrelation (section #5 in the script).

###01_script_compute_variables_occ.R

We provide another script that has been used to compute landscape metrics used as landscape variables in the occupancy model. This script will not run as we do not provide raw spatial data. We nonetheless provide the table ‘table_variables_Bombina.pdf’ with the description of the raw spatial data that were used to compute spatial variables.

###02_script_selection_variables_occ.R

Finally, we provide a script to select landscape variables in several steps:
VIF, Pearson's correlation, PCA, and linear models for each variable. This script uses the table "table_variables_SVF_2024.txt” provided. This table has been generated with the previous script.

## FRANÇAIS

Dans cette archive, vous trouverez un script pour ajuster un modèle d'occupation pluriannuel sur un ensemble de données de suivi de deux ans de Bombina bombina en Moselle (France), et pour évaluer l'autocorrélation spatiale.

###script_occupancy_Bombina_French_2024.R

La collecte des données est basée sur un protocole standardisé sur un échantillon spatial comprenant 120 cellules de 200 x 200 placées au hasard à l'intérieur et à l'extérieur de l'aire de répartition connue de l'espèce en Moselle dans le secteur d'Albestroff au début du suivi, c'est-à-dire 2022, afin d'évaluer les changements de distribution au cours du temps. 
Chaque cellule a été visitée deux fois par an et le suivi se déroule tous les deux ans. La première année de surveillance était 2022, et il a été réalisé en 2024, ce qui en fait un ensemble de données sur deux ans.
Nous avons utilisé le package unmarked pour ajuster un modèle d'occupation pluriannuel à cet ensemble de données afin d'évaluer la détection, l'occupation pluriannuelle, la colonisation et l'extinction.
Nous donnons une table de données déjà transformées, qui est utilisée pour le script. Nous fournissons également une couche SIG des unités d'échantillon pour estimer l'autocorrélation spatiale (section #5 du script).

###01_script_compute_variables_occ.R

Nous fournissons également un script  qui a été utilisé pour calculer les métriques paysagères utilisée comme variables paysagères dans le modèle d'occupation. Ce script ne s'exécutera pas car les sources de données ne sont pas fournies. Nous fournissons le tableau ‘table_variables_Bombina.pdf’ qui décrit les données brutes spatiales que nous avons utilisé pour générer les variables paysagères.

###02_script_selection_variables_occ.R

Enfin, nous fournissons un script de la procédure de sélection des variables paysagères en plusieurs étapes : VIF, corrélation de Pearson, ACP et modèle linéaire sur chacune des variables. Ce script utilise la table "table_variables_SVF_2024.txt” issue de l'étape précédente de calcul des métriques paysagères. Cette table est fournie.

