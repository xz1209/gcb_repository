The datasets and associated resources supporting the manuscript submitted to Global Change Biology


Title of manuscript:
"From deficiency to enrichment: diverging soil phosphorus trajectories across China’s croplands (1980–2018)"


The materials are organized into three major categories:
1. OVERVIEW
This repository contains the datasets, R scripts, and map outputs used in the study investigating long-term changes in soil available phosphorus (AP) across China’s croplands from 1980 to 2018.
The materials are organized into three major categories:
1) Raw Data
2) R Scripts
3) Created Maps
All files are provided to support transparency, reproducibility, and reuse of the published results.


2. ORIGINAL DATA
All Original datasets are compiled in the Excel file: Supplementary_data.xlsx. 
This file contains multiple worksheets corresponding to different datasets used in this study.

2.1 National Soil Survey Data (Three Periods)
The following worksheets contain observed AP data from three national survey periods:
1) Sheet: National soil surveys_1980
Soil AP observations representing the historical baseline period (1980)
2) Sheet: National soil surveys_2012
Soil AP observations representing the intermediate survey period (2012)
3) Sheet: National soil surveys_2018
Soil AP observations representing the intermediate survey period (2018)

2.2 Long-Term Fertilization Experiment Data
This worksheet contains observations from long-term fertilization experiments used to support the interpretation of national phosphorus trajectories.
Sheet: Long-term experiments

2.3 Random Forest Modeling Dataset
This worksheet contains the integrated dataset used for machine-learning model development.
Sheet: APrate_modeling_data


3. R SCRIPTS
This folder contains the R code used for machine learning modeling, prediction, and interpretation.

3.1 Cross-Validation + Bayesian Optimization + Sampling Density Correction
File: R1.R
Main functions: Five-fold cross-validation; Bayesian hyperparameter optimization; Sampling density correction; Weighted model evaluation

3.2 Random Forest Model Building
File: R2.R
Main functions: Data preprocessing; Train/test split; Variable standardization; Random forest model fitting; Prediction; Performance evaluation

3.3 Random Forest + SHAP Analysis
File: R3.R
Main functions: Random forest model fitting; SHAP value calculation; Ariable importance; Local interpretation of model predictions


4. CREATED MAPS
This folder contains the GeoTIFF raster maps generated in this study.
The GeoTIFF files include:
1) APrate_1980_2018.tif
National map of the annual rate of change in soil available P (APrate) from 1980 to 2018.
2) APrate_prediction.tif
National map of APrate predicted by the random forest model.
3) APrate_dominant_drivers.tif
Spatial distribution of the dominant drivers of APrate.
Note: Raster values indicate dominant driver categories:
  1 = Environment
  2 = Climate
  3 = Vegetation
  4 = Management
  5 = Climate + Environment
  6 = Environment + Vegetation
  7 = Environment + Management
  8 = Climate + Vegetation
  9 = Climate + Management
  10 = Management + Vegetation
4) APrate_sensitivity_risk.tif
Spatial distribution of cropland sensitivity-risk areas classified according to the relative effects of management practices and climate factors on APrate.
Note: Raster values indicate sensitivity-risk categories:
  1 = Strongly management-buffered
  2 = Weakly management-buffered
  3 = Weakly climate-constrained
  4 = Strongly climate-constrained


6. SOFTWARE REQUIREMENTS
Recommended software: R (version 4.4.1), ArcGIS (version 10.2)
Main R packages: data.table; ranger; mlr3; mlr3tuning; fastshap


7. CONTACT
For questions regarding the datasets, scripts, or related materials, please contact the corresponding author:
Email: xuzhen1209@126.com



