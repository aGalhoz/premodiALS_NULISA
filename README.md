In this project it is implemented a workflow to analyse proteomic profiles of a NULISA assay. The data corresponds to samples from Germany and Turkey of the premodiALS project from the first visit. The data comprises two panels: a CNS (~120 proteins) and immune panel (~270 proteins), and of three fluids (plasma, serum and CSF).

The code available at this repository is divided into the following main steps:

* data mining: processing of proteomic and diseased profiles, analysis of covariates, exploration of normalisation and target detectability;
* group analysis: implementation of the Limma linear model to identifiy differentially expresseed proteins between groups, adjusted for age, sex and center. includes boxplots and PCA visualisations;
* differential expression analyses: volcano plots and signed p-value visualisations between all comparisons and for each fluid;
* analyses of proteins of interest with clinical variables (ALSFRS-R, disease duration, age, sex, site of onset). for significant correlations, a detailed exploration of the relationship between ALSFRS-R subscores and specific proteins was made;
* ML modelling between groups: Lasso modelling between ALS vs CTR, ALS vs PGMC and PGMC vs CTR, across 500 bootstrap iterations with an internal 5-fold cross-validation for hyperparameter tuning. extraction of relevant proteins, and definition of a protein signature based on ALS vs CTR model. Usage of the reeduced protein signature to predict an ALS risk score among PGMC individuals, with clustering analyses, PCA and heatmap visualisation.
