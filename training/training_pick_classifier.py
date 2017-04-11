#!/usr/bin/python3

import sys
import os
import numpy as np
import pandas as pd

from sklearn.preprocessing import LabelEncoder
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split


from sklearn.metrics import accuracy_score, log_loss, recall_score, precision_score
from sklearn.neighbors import KNeighborsClassifier
from sklearn.svm import SVC, LinearSVC, NuSVC
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.naive_bayes import GaussianNB

from sklearn.feature_selection import f_classif
from sklearn.feature_selection import SelectKBest
import pickle

root = '../merge_data'


def load_data( standardize=False ):

    '''
    load the data to the pandas dataframe and returns
    a tuple of the cnv ids, the data, and the labels.
    '''
    data = pd.read_table(os.path.join( root, 'train_cnv_w_gene_score.txt' ))
    #data.loc[data['ClinicalSignificance'] == "Benign", 'ClinicalSignificance'] = 0
    #data.loc[data['ClinicalSignificance'] == 'Pathogenic', 'ClinicalSignificance'] = 1
    id = pd.concat([data.pop('Chr'), data.pop('Start'), data.pop('End')], axis = 1)
    y = data.pop('ClinicalSignificance')
    y = LabelEncoder().fit(y).transform(y)
    data.drop([col for col in ['genes', 'gene_breaked_start', 'gene_breaked_end', 'disease_genes']], axis = 1, inplace = True)
    data.loc[data['VarType'] == 'gain', 'VarType'] = 0
    data.loc[data['VarType'] == 'loss', 'VarType'] = 1
    data.loc[data['num_disease_genes'] > 0, 'num_disease_genes'] = 1
    data.loc[data['num_disease_genes'] == 0, 'num_disease_genes'] = 0
    if standardize:
        data = StandardScaler().fit(data).transform(data)

    return id, y, data


cnv_ID, data_labels, data_content = load_data( standardize=False )

'''
sellect the most important features for training
'''
#print(data_content.head(20))

features = list(data_content.columns.values)
num_feature = int(sys.argv[1])

selector = SelectKBest(f_classif, k=num_feature)
features_new = selector.fit_transform(data_content, data_labels)
print(selector.scores_)


'''
print out the features ranked by more much variants of the
results the could be explained by the features.
'''
zipped = zip(features, selector.scores_)
zipped_sort = sorted(zipped, key = lambda t:t[1], reverse = True)

new_features_list = []

i = 0
while i<= num_feature -1:
    new_features_list.append(zipped_sort[i][0])
    i += 1
print(new_features_list)

'''
split the data into training and testing dataset
'''

#features_train, features_test, labels_train, labels_test = train_test_split(
#                    data_content, data_labels, test_size = 0.3, random_state = 42)



def try_different_clf():
    '''
    Test the performance of different classifier
    '''

    classifiers = [KNeighborsClassifier(3),
              SVC(kernel='rbf', C=0.025, probability=True),
              NuSVC(probability=True),
              DecisionTreeClassifier(),
              RandomForestClassifier(),
              GradientBoostingClassifier(),
              GaussianNB()]
    log_cols = ['Classifier', 'Accuracy', 'Recall', 'Precision', 'log_loss']
    log = pd.DataFrame(columns=log_cols)

    for clf in classifiers:
        clf.fit(features_train, labels_train)
        name = clf.__class__.__name__
        print('='*30)
        print(name)
        print('******Result*******')
        predictions = clf.predict(features_test)
        acc = accuracy_score(labels_test, predictions)
        print("Accuracy:{:.4%}".format(acc))
        recall = recall_score(labels_test, predictions)
        print("Recall:{:.4%}".format(recall))
        prec = precision_score(labels_test, predictions)
        print("Precision:{:.4%}".format(prec))

        predictions = clf.predict_proba(features_test)
        ll = log_loss(labels_test, predictions)
        print("Log Loss: {}".format(ll))

    print('='*30)



'''
After the initial selection on the classifier, perform KFold for cross_validation
to narrow down to the classifier that perform best.
'''

from sklearn.model_selection import StratifiedKFold


classifiers = [DecisionTreeClassifier(),
                RandomForestClassifier(),
                GradientBoostingClassifier()]

accuracys = []
precisions = []
recalls = []
logs = []

skf = StratifiedKFold(n_splits=15, random_state = 1, shuffle = True)

temp_feature = ['pli_max', 'mis_z_max', 'num_gene', 'num_disease_genes', 'Length', 'VarType', 'end_break_gene_mis_z', 'end_break_gene_pli', 'start_break_gene_mis_z']
temp_feature2 = ['pli_max', 'mis_z_max', 'num_gene', 'num_largeZmis', 'num_disease_genes', 'Length', 'VarType', 'num_largePLI']
for clf in classifiers:
    for train_index, test_index in skf.split(data_content, data_labels):
        train_predictors, test_predictors = data_content[new_features_list].iloc[train_index], data_content[new_features_list].iloc[test_index]
        train_labels, test_labels = data_labels[train_index], data_labels[test_index]
        clf.fit(train_predictors, train_labels)
        predictions = clf.predict(test_predictors)
        acc_score = accuracy_score(test_labels, predictions)
        precision = precision_score(test_labels, predictions)
        recall = recall_score(test_labels, predictions)
        accuracys.append(acc_score)
        precisions.append(precision)
        recalls.append(recall)
        predictions = clf.predict_proba(test_predictors)
        ll = log_loss(test_labels, predictions)
        logs.append(ll)
    name = clf.__class__.__name__

    print('='*30)
    print(name)
    print('******Result*******')
    print("Accuracy:{:.4%}".format(sum(accuracys)/len(accuracys)))
    print("Precisions:{:.4%}".format(sum(precisions)/len(precisions)))
    print("Recall:{:.4%}".format(sum(recalls)/len(recalls)))
    print("Log loss: {}".format(sum(logs)/len(logs)))

print('='*30)

clf = RandomForestClassifier()
clf2 = GradientBoostingClassifier()


#for train_index, test_index in skf.split(data_content, data_labels):
#    train_predictors, test_predictors = data_content[new_features_list].iloc[train_index], data_content[new_features_list].iloc[test_index]
#    train_labels, test_labels = data_labels[train_index], data_labels[test_index]
#    clf.fit(train_predictors, train_labels)

clf.fit(data_content[new_features_list], data_labels)
clf2.fit(data_content[new_features_list], data_labels)
with open('RandomForestClassifier.pickle', 'wb') as f:
    pickle.dump(clf, f)
with open('GradientBoostingClassifier.pickle', 'wb') as f:
    pickle.dump(clf2, f)
