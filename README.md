# cnvPRED

CNV annotation and interpretation tool

## Introduction
Copy-number variants (CNVs) are defined as segments of the genome larger than 1kb that are either deleted or duplicated. CNVs can be benign polymorphic variations or can significantly affect phenotypic variability. CNVs have been associated with a variety of diseases, especially cancer, autism, schizophrenia, and developmental delay. With the rapid evolution of different CNV detection methodologies,  increased resolutions of CNV detection, and large sample collections,   tools for the computational prioritization of pathogenic CNVs will be necessary. Therefore, I developed a python3-based tool (cnvPRED) to improve the clinical utility of computational CNV prediction to identify the pathogenic CNVs from large sample collections. 

cnvPRED consists of two stages: CNV annotation and CNV interpretation. The annotation stage is designed to annotate predicted CNVs with functional information using different source datasets, including ExAC, OMIM, HGMD, and 1000G (cnv_annotation.py). The CNV interpretation stage is designed to use the Machine learning methods and use the annotated features to predict is the CNV is pathogenic or benign. The training data were downloaded from ClinGen CNV database. Different classifiers were tested using training/training_pick_classifiers.py and training/generate_roc.py. 

## Basic Manual 
1. CNV annotation
python3 cnv_annotation.py  -i input_file -o your_output
2. CNV interpretation
python3 cnv_prediction.py -i (the output of the annoation stage)

The input file is a bed file (tab delimited) with CNV position, sample ID and the estimated length of CNVs:
eg:

Chr     Start           End           CNV_type    SampleID Length

1       1635261         1670548       DEL         A01       35287

I'll continue to improve cnvPRED. If you have good suggestions please let me know! 
