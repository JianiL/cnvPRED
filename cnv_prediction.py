#!/usr/bin/env python3

import sys
import pickle
import os
import numpy as np
import pandas as pd


from sklearn.preprocessing import LabelEncoder

#root = '../data'

input = sys.argv[1]

def load_data( standardize=False ):
    data = pd.read_table(input)
    id = pd.concat([data.pop('Chr'), data.pop('Start'), data.pop('End')], axis = 1)
    genes = data.pop('genes')
    data.drop([col for col in ['gene_breaked_start', 'gene_breaked_end']], axis = 1, inplace = True)
    samples = data.pop('SampleID')
    cnv_type = data['VarType']
    data.loc[data['VarType'] == 'DUP', 'VarType'] = 0
    data.loc[data['VarType'] == 'DEL', 'VarType'] = 1
    data.loc[data['num_disease_genes'] > 0, 'num_disease_genes'] = 1
    data.loc[data['num_disease_genes'] == 0, 'num_disease_genes'] = 0
    if standardize:
        data = StandardScaler().fit(data).transform(data)

    return id, data, samples, cnv_type, genes

cnv_id, cnv_data, sample_id, cnv_type, gene = load_data( standardize=False )

feature_list = ['num_disease_genes', 'pli_max', 'mis_z_max', 'end_break_gene_mis_z', 'VarType', 'num_gene', 'num_largeZmis', 'Length', 'end_break_gene_pli', 'num_largePLI', 'start_break_gene_mis_z', 'start_break_gene_pli']

pickle_in = open('resources/RandomForestClassifier.pickle', 'rb')
clf = pickle.load(pickle_in)

predictions = clf.predict(cnv_data[feature_list])

col = [ 'predictions']

pred_result = pd.DataFrame(predictions, columns=col)
raw_output = pd.read_table(input)
output = pd.concat([raw_output, pred_result], axis=1)
output.reset_index()

output_file = sys.argv[2]

output.to_csv(output_file, index = False)
