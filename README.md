# DRdriver
This repository contains the files and scripts described in "DRdriver: identifying drug resistance driver genes using individual-specific gene regulatory network". DRdriver is a method to identify personalized drug resistance driver gene. Through integrating mutation, expression and gene regulatory network, the genes with differential mutations which regulated the most differentially expressed genes and the least non-differentially expressed genes were predicted as driver genes for each resistant patient.

DRdriver.r: The main function for identifying drug resistance driver genes, which have already included all scripts for running the whole pipeline. We also provided the example.r as an example to run DRdriver.r. In example.r, we took the condition LGG_Temozolomide as an example to show how to get the driver genes based on DRdriver.r. Users can change the condition to obtain the driver genes of other conditions.

differential_exp.r: The function for identifying DEGs.

differential_mutation.r: The function for identifying differential mutations.

genetic_algorithm.r: The pipeline for genetic algorithm used in this study.

specific_network.r: The function for constructing individual-specific network.

GRN.txt: The file contains gene regulatory network. 

patient_list.txt: The file contains the resistant sample IDs that can used for identifying driver genes. Several patients which who have too few overlapped genes in candidate genes and regulators or whose individual-specific network do not contain DEGs were deleted.

Except the data of LGG_Temozolomide, the folder data contained all input data of 13 conditions. Because there were several random parameters in genetic algorithm, such as mutation and crossover rate, the results of different runs were not the same. The similarity was almost over 90%.
