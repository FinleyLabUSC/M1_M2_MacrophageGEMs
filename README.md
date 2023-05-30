# M1_M2_MacrophageGEMs
Genome-scale Metabolic Modeling of CRC-associated Macrophage subtypes
Existing Toolboxes Needed:
Cobra Toolbox (https://github.com/opencobra) 
REDGEM and REDHUMAN (https://github.com/EPFL-LCSB) 
PageRank Analysis (https://www.nature.com/articles/s41540-018-0067-y)
REDGEM pipeline seems to work best with TOMLAB, but is not necessary (https://tomopt.com/tomlab/) 

Use "GenerateCohort_MacrophageModels" to make the set of models, using the different techniques.
All but CORDA are included with the CobraToolbox (https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004808)

The "M1_xxxx" and "M2_xxxx" files are the resulting models
The "redHUMAN_xxx" files are the consolidated M1 and M2 models, built using REDHUMAN

Model analyses can be performed with the "AnalysesMacrophageModels" file

