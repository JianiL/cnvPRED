import numpy as np
import os
import sys
import pandas as pd

from sklearn.preprocessing import LabelEncoder
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split


from sklearn.metrics import accuracy_score, log_loss, recall_score, precision_score, roc_curve, roc_auc_score
from sklearn.neighbors import KNeighborsClassifier
from sklearn.svm import SVC
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import OneHotEncoder
from sklearn.ensemble import VotingClassifier


import matplotlib.pyplot as plt
from matplotlib import style

style.use('ggplot')


root = '../data'


def load_data( standardize=False ):

    '''
    load the data to the pandas dataframe and returns 
    a tuple of the cnv ids, the data, and the labels.
    '''
    data = pd.read_table(os.path.join( root, 'CNV_w_gene_w_score_wBreak.txt' ))
    #data.loc[data['ClinicalSignificance'] == "Benign", 'ClinicalSignificance'] = 0
    #data.loc[data['ClinicalSignificance'] == 'Pathogenic', 'ClinicalSignificance'] = 1
    id = pd.concat([data.pop('Chr'), data.pop('Start'), data.pop('End')], axis = 1)
    y = data.pop('ClinicalSignificance')
    y = LabelEncoder().fit(y).transform(y)
    data.drop([col for col in ['genes', 'gene_breaked_start', 'gene_breaked_end']], axis = 1, inplace = True)
    data.loc[data['VarType'] == 'gain', 'VarType'] = 0
    data.loc[data['VarType'] == 'loss', 'VarType'] = 1
    if standardize:
        data = StandardScaler().fit(data).transform(data)
    
    return id, y, data


cnv_ID, data_labels, data_content = load_data( standardize=False )


features_train, features_test, labels_train, labels_test = train_test_split(
                    data_content, data_labels, test_size = 0.3, random_state = 42)


features_train1, features_train2, labels_train1, labels_train2 = train_test_split(features_train, 
                                                                                    labels_train,
                                                                                    test_size = 0.5)


knn = KNeighborsClassifier(3)
knn.fit(features_train, labels_train)

knn_prediction =  knn.predict_proba(features_test)[:, 1]
fpr_knn, tpr_knn, _ = roc_curve(labels_test, knn_prediction)
knn_auc = roc_auc_score(labels_test, knn_prediction)

lr = LogisticRegression()
lr.fit(features_train, labels_train)

lr_prediction = lr.predict_proba(features_test)[:, 1]
fpr_lr, tpr_lr, _ = roc_curve(labels_test, lr_prediction)
lr_auc = roc_auc_score(labels_test, lr_prediction)


rf = RandomForestClassifier()
rf.fit(features_train, labels_train)

rf_prediction = rf.predict_proba(features_test)[:, 1]
fpr_rf, tpr_rf, _ = roc_curve(labels_test, rf_prediction)
rf_auc = roc_auc_score(labels_test, rf_prediction)
rf_prediction = rf.predict(features_test)
rf_acc = accuracy_score(labels_test, rf_prediction)

#eclf = VotingClassifier(estimators=[('lr', lr), ('rf', rf), ('knn', knn)], voting = 'hard')
#eclf.fit(features_train, labels_train)

#eclf_predictions = eclf.predict(features_test)
#eclf_acc = accuracy_score(labels_test, eclf_predictions)
#fpr_eclf, tpr_eclf, _ = roc_curve(labels_test, eclf_predictions)

rf_enc = OneHotEncoder()
rf_lm = LogisticRegression()
rf.fit(features_train1, labels_train1)
rf_enc.fit(rf.apply(features_train1))
rf_lm.fit(rf_enc.transform(rf.apply(features_train2)), labels_train2)

rf_lm_predictions = rf_lm.predict_proba(rf_enc.transform(rf.apply(features_test)))[:, 1]
fpr_rf_lm, tpr_rf_lm, _ = roc_curve(labels_test, rf_lm_predictions)
rf_lm_auc = roc_auc_score(labels_test, rf_lm_predictions)
rf_lm_predictions = rf_lm.predict(rf_enc.transform(rf.apply(features_test)))
rf_lm_acc = accuracy_score(labels_test, rf_lm_predictions)



svc = SVC(kernel='rbf', C=0.025, probability=True)
svc.fit(features_train, labels_train)

svc_prediction = svc.predict_proba(features_test)[:, 1]
fpr_svc, tpr_svc, _ = roc_curve(labels_test, svc_prediction)
svc_auc = roc_auc_score(labels_test, svc_prediction)


gb = GradientBoostingClassifier()
gb.fit(features_train, labels_train)

gb_prediction = gb.predict_proba(features_test)[:, 1]
fpr_gb, tpr_gb, _ = roc_curve(labels_test, gb_prediction)
gb_auc = roc_auc_score(labels_test, gb_prediction)

#gb_enc = OneHotEncoder()
#gb_lm = LogisticRegression()
#gb.fit(features_train1, labels_train1)
#gb_enc.fit(gb.apply(features_train1)[:, :, 0])
#gb_lm.fit(gb_enc.transform(gb.apply(features_train2)[:, :, 0]), labels_train2)

#gb_lm_predictions = gb_lm.predict_proba(gb_enc.transform(gb.apply(features_test)[:, :, 0]))[:, 1]
#fpr_gb_lm, tpr_gb_lm, _ = roc_curve(labels_test, gb_lm_predictions)
#gb_lm_auc = roc_auc_score(labels_test, gb_lm_predictions)

#print('KNN AUC: {}'.format(knn_auc))
#print("Random Forest AUC: {}".format(rf_auc))
#print('Logestic regression AUC: {}'.format(lr_auc))
#print('RF + LR AUC: {}'.format(rf_lm_auc))
#print('GradientBoosting AUC: {}'.format(gb_auc))
#print('GB + LR AUC: {}'.format(gb_lm_auc))
#print('Accuracy ensamble: {}'.format(eclf_acc))
#print('Accuracy RF: {}'.format(rf_acc))
#print('Accuracy RF+LM: {}'.format(rf_lm_acc))


plt.figure()
plt.plot([0, 1], [0, 1], 'k--')
plt.plot(fpr_knn, tpr_knn, label='KNN (area = %0.3f)' % knn_auc)
plt.plot(fpr_lr, tpr_lr, label='LR (area = %0.3f)' % lr_auc )
plt.plot(fpr_rf, tpr_rf, label='RF (area = %0.3f)' % rf_auc )
plt.plot(fpr_svc, tpr_svc, label="SVC (area = %0.3f)" % svc_auc )
#plt.plot(fpr_rf_lm, tpr_rf_lm, label='RF + LR (area = %0.3f)' % rf_lm_auc)
plt.plot(fpr_gb, tpr_gb, label='GB (area = %0.3f)' % gb_auc)
#plt.plot(fpr_gb_lm, tpr_gb_lm, label='GB + LR')
plt.xlabel('False postivie rate')
plt.ylabel('True positive rate')
plt.title('Roc curve')
plt.legend(loc='best')
plt.show()


