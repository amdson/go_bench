import numpy as np
from scipy.sparse import csc_matrix, csr_matrix, dok_matrix, hstack, vstack
from sklearn.metrics import multilabel_confusion_matrix, precision_recall_fscore_support


import math
from collections import Counter

#Calculates relative information gain of each term found in set of protein annotations.
#Based on ratio between specificity of most specific parent, and specificity of term. 
def calculate_ic(annots, ontology):
    cnt = Counter()
    for x in annots.values():
        cnt.update(x)
    all_terms = set()
    for t in annots.values():
        all_terms.update(t)
    ic = {}
    for go_id, n in cnt.items():
        if(not go_id in ontology):
            ic[go_id] = 0
            continue
        parents = ontology[go_id]._parents
        parents = set(parents).intersection(all_terms)
        if len(parents) == 0:
            min_n = n
        else:
            min_n = max(min([cnt[x] for x in parents]), n)
        ic[go_id] = math.log(min_n / n, 2)
    return ic

#Convert ic dict to ic matrix for calculating s_metric. 
def ic_mat(terms, ic_dict):
    term_ic = np.zeros(len(terms))
    for i, term in enumerate(terms):
        term_ic[i] = ic_dict[term]
    return term_ic

def s_metric(testing_matrix, prediction_matrix, test_ia):
    conf_mat = multilabel_confusion_matrix(testing_matrix, prediction_matrix)
    fp = conf_mat[:, 0, 1]
    fn = conf_mat[:, 1, 0]
    mi = np.dot(fp, test_ia) / testing_matrix.shape[0]
    ru = np.dot(fn, test_ia) / testing_matrix.shape[0]
    return mi, ru, np.sqrt(mi*mi + ru*ru)

def threshold_stats(testing_matrix, prediction_matrix, test_ia):
    precs = []
    recs = []
    f_scores = [] 
    mis = []
    rus = []
    s_vals = []
    rms = []
    for threshold in np.linspace(0.001, 1, 100):
        preds = prediction_matrix.copy()
        preds.data = np.where(preds.data >= threshold, 1, 0)
        preds.eliminate_zeros()
        p, r, f, support = precision_recall_fscore_support(testing_matrix, preds, average='micro')
        mi, ru, s = s_metric(testing_matrix, preds, test_ia)
        precs.append(p)
        recs.append(r)
        f_scores.append(f)
        mis.append(mi)
        rus.append(ru)
        s_vals.append(s)
        rms.append(r*r * preds.shape[0] * preds.shape[1] / preds.sum())
    return precs, recs, f_scores, rms, mis, rus, s_vals
