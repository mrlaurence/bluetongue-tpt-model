# README

This repository contains all the code used in the R-Project "A Compartmental Disease Model to Describe Long-Term Bluetongue Virus Behaviour, with Analysis of Overwintering in France".

To run a simulation of the Bluetongue_TPT_Model, configure the simulation parameters in Bluetongue_TPT_Model/Run.m (and Bluetongue_TPT_Model.Define_Within_Farm_Parameters() if you wish to change within-farm parameters) and then run Bluetongue_TPT_Model.Run() in MATLAB 2023a.

All files needed for the compartmental model are contained in Bluetongue_TPT_Model. Files in Other were used for smaller experiments in the course of the project and are provided for completeness but without documentation.

# Licenses and Attribution

All .m files were written by Laurence Dhonau, with the exception of haversine.m, which was written by Josiah Renfree, obtained from MathWorks MATLAB Central File Exchange and included under the license provided in the file.

The compartmental model is based on the model presented in the following paper:
Sumner T, Orton RJ, Green DM, Kao RR, Gubbins S. Quantifying the roles of host movement and vector dispersal in the transmission of vector-borne diseases of livestock. PLoS Comput Biol. 2017 Apr 3;13(4):e1005470. doi: 10.1371/journal.pcbi.1005470. PMID: 28369082; PMCID: PMC5393902.

The code in the Bluetongue_TPT_Model folder largely implements the model described by Sumner et al. but with (very substantial) modifications, which are described in the project report. No code was used directly and the moedl was adapted under Attribution 4.0 International (CC BY 4.0) (see https://creativecommons.org/licenses/by/4.0/).

Climate data in Datasets/Temperature_Raw were extracted from the European Climate Assessment & Dataset (EC&D) and used under the following license:

EUROPEAN CLIMATE ASSESSMENT & DATASET (ECA&D), file created on 22-03-2023
THESE DATA CAN BE USED FREELY PROVIDED THAT THE FOLLOWING SOURCE IS ACKNOWLEDGED:

Klein Tank, A.M.G. and Coauthors, 2002. Daily dataset of 20th-century surface
air temperature and precipitation series for the European Climate Assessment.
Int. J. of Climatol., 22, 1441-1453.
Data and metadata available at http://www.ecad.eu

We constructed the following files in the Datasets folder based on the climate dataset: Farm_Square_Nearest_Station_Dictionary.mat, Station_Maximum_Temperature_Table.mat, Station_Minimum_Temperature_Table.mat, Station_Table.mat, Station_Temperature_Dictionary.mat.

The files GLW3_CATTLE_DA_FRANCE_POINTS.geojson and GLW3_SHEEP_DA_FRANCE_POINTS.geojson were produced using the QGIS software packaged from the Gridded Livestock of the World - 2010 (GLW 3) dataset. This is downloadable from the Harvard Dataverse here: https://dataverse.harvard.edu/dataverse/glw. The dataset is in the public domain (https://creativecommons.org/publicdomain/zero/1.0/). Details of the GLW dataset may be found in the following paper:
Gilbert M, Nicolas G, Cinardi G, Van Boeckel TP, Vanwambeke SO, Wint GRW, Robinson TP. Global distribution data for cattle, buffaloes, horses, sheep, goats, pigs, chickens and ducks in 2010. Sci Data. 2018 Oct 30;5:180227. doi: 10.1038/sdata.2018.227. PMID: 30375994; PMCID: PMC6207061.
