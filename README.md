# DRdriver
This repository contains the files and scripts described in "DRdriver: identifying drug resistance driver genes using individual-specific gene regulatory network".

DRdriver.r: The main function for identifying drug resistance driver genes.

differential_mutation.r: The function for identifying differential mutations.

differential_exp.r: The function for identifying DEGs.

specific_network.r: The function for constructing individual-specific network.

genetic_algorithm.r: The pipeline for genetic algorithm used in this study.

example.r: Taking the condition LGG_Temezolomide as an example to identify driver genes based on DRdriver.r, users can changing the condition to obtain the driver genes of other conditions.

patient_list.txt: The file contains the resistant sample IDs that can used for identifying driver genes. Several patients which who have too few overlapped genes in candidate genes and regulators or whose individual-specific network do not contain DEGs were deleted.

GRN.txt: The file contains gene regulatory network. 

Except the data of LGG_Temozolomide, the folder data contained all input data of 13 conditions.

Because there were several random parameters in genetic algorithm, such as mutation and crossover rate, the results of different runs were not the same. The similarity was almost over 90%.
