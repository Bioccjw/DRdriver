# DRdriver
This repository contains the files and scripts described in "DRdriver: identifying drug resistance driver genes using individual-specific gene regulatory network".

differential_mutation.r: The script for identifying differential mutations.

differential_expression.r: The script for identifying differential expressed genes.

specific_network.r: The script for constructing individual-specific network.

genetic_algorithm.r: The pipline for genetic algorithm used in this study.

DRdriver: The main script for identifying drug resistance driver genes.


Taking the condition of LGG_Temezolomide as an example:

sample_list.txt: The file contains the resistant sample IDs that can used for identifying driver genes. 8 samples (TCGA-VM-A8CE, TCGA-DB-A4XG, TCGA-DB-A64P, TCGA-FG-A710, TCGA-HT-7884, TCGA-HW-8320, TCGA-DU-A76R, TCGA-TM-A84L) who have too few intersected genes in candidate genes and regulators or whose individual-specific network do not contain DEGs were deleted. The users can get the corresponding information by providing the sample ID in corresponding script.

sensitive_mutation.txt and resistant_mutation.txt provide the mutation information for the condition. Before the users run the script of differential_mutation.r, they must load these two files first.

sensitive_exp.txt and resistant_exp.txt provide the expression profiles for the condition. Before the users run the script of differential_expression.r, they must load these two files first.

GRN.txt: The file contains gene regulatory network. This file must load when the users run the script specific_network.r. Additionally, differential_mutation.r and differential_expression.r should be loaded before running specific_network.r.

After loading all files and scripts provided, users can provide the sample ID to get the driver genes by run DRdriver.r.

Since there were several random parameters in genetic algorithm, the results of different run of the script might not be the same. But the similarity of them almost over 80%.
