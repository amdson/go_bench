import numpy as np
from scipy.sparse import csc_matrix, csr_matrix, dok_matrix, hstack, vstack
from sklearn.metrics import multilabel_confusion_matrix, precision_recall_fscore_support


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
