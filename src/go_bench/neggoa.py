from go_bench.utils import list_ancestors
import numpy as np

godag = None

def inclusive_ancestors(term):
    ancestors = set(godag[term].get_all_parents())
    ancestors.add(term)
    return ancestors
        
def lca(t1, t2):
    s1 = set(list_ancestors(t1))
    s2 = set(list_ancestors(t2))
    shared_ancestors = s1.intersection(s2)
    max_depth_term = None
    max_depth = 0
    for a in shared_ancestors:
        if(godag[a].depth >= max_depth):
            max_depth_term = a
            max_depth = godag[a].depth
    return max_depth_term

def calc_info_content(selGOs, nTerms):
    log_terms = np.log(nTerms)
    GO_ICs =np.ones(len(selGOs))
    for i, term in enumerate(selGOs): 
        descendants = godag[term].get_all_children()
        if(len(descendants) > 0):
            GO_ICs[i] = 1-np.log(len(descendants) + 1)/log_terms
    return GO_ICs

go_terms = [x for x in godag]
go_info_content = calc_info_content(go_terms, len(godag))

#Connect terms only to parents or children in sparse matrix. 
def get_connected_term_similarities(selGOs, nTerms):
    ancestor_terms = {x:inclusive_ancestors(x) for x in selGOs}
    term_col_mapping = {term:i for i, term in enumerate(selGOs)}
    go_info_content = calc_info_content(selGOs, nTerms)
    sim_mat = lil_matrix((len(selGOs), len(selGOs)))
    sel_term_set = set(selGOs)
    for i, term in enumerate(selGOs):
        term_ancestors = ancestor_terms[term]
        parents = godag[term]._parents
        for parent in parents:
            j = term_col_mapping[parent]
            parent_ancestors = ancestor_terms[parent]
            shared_ancestors = term_ancestors.intersection(parent_ancestors)
            if(len(shared_ancestors) == 0):
                print(term, parent, "terms have no common ancestor")
            pca = max(go_info_content[term_col_mapping[x]] for x in shared_ancestors)
            if(pca != go_info_content[j]):
                print("term {} has weird shared ancestors {}".format(term, parent))
            sim = 2*pca / (go_info_content[i] + go_info_content[j])
            sim_mat[i, j] = sim
            sim_mat[j, i] = sim
        sim_mat[i, i] = 1
    return sim_mat.tocsr()

#Connect terms only to parents or children in sparse matrix. 
def get_connected_term_similarities(selGOs, nTerms):
    ancestor_terms = {x:inclusive_ancestors(x) for x in selGOs}
    term_col_mapping = {term:i for i, term in enumerate(selGOs)}
    go_info_content = calc_info_content(selGOs, nTerms)
    sim_mat = dok_matrix((len(selGOs), len(selGOs)))
    sel_term_set = set(selGOs)
    for i, term in enumerate(selGOs):
        term_ancestors = ancestor_terms[term]
        parents = godag[term]._parents
        for parent in parents:
            j = term_col_mapping[parent]
            parent_ancestors = ancestor_terms[parent]
            shared_ancestors = term_ancestors.intersection(parent_ancestors)
            if(len(shared_ancestors) == 0):
                print(term, parent, "terms have no common ancestor")
            pca = max(go_info_content[term_col_mapping[x]] for x in shared_ancestors)
            if(pca != go_info_content[j]):
                print("term {} has weird shared ancestors {}".format(term, parent))
            sim = 2*pca / (go_info_content[i] + go_info_content[j])
            sim_mat[i, j] = sim
            sim_mat[j, i] = sim
        sim_mat[i, i] = 1
    return sim_mat.tocsr()

def RWR(transMat, alpha, maxIter):
    num_terms = transMat.shape[0]
    eyeMat = sparse.eye(num_terms)
    oldW = eyeMat
    newW = oldW
    threshold = eps
    i = 0
    delta = 10*threshold
    while delta>threshold and i < maxIter:
        newW = eyeMat*alpha + (oldW@transMat)*(1-alpha)
        delta = np.sum(newW-oldW)
        i += 1
    return newW

def getProbSim(protein_annotation_dict, prot_ids, term_col_mappings):
    prot_mat = dok_matrix((len(prot_ids), len(term_col_mappings)))
    count = 0
    for i, prot_id in enumerate(prot_ids):
        terms = protein_annotation_dict[prot_id]
        for term in terms:
            if(term in term_col_mappings):
                prot_mat[i, term_col_mappings[term]] += 1
        count += 1
        if(count % 20000 == 0):
            print(f"Percentage progress: {count / len(protein_annotation_dict) * 100}%")
            print(len(terms))    
    prot_mat = prot_mat.tocsc()
    prot_mat_T = prot_mat.transpose()
    coocurrence_mat = prot_mat_T @ prot_mat
    return coocurrence_mat

def normalize_trans_matrix(prob_sim):
    row_sums = prob_sim.sum(axis=1)
    row_sums[row_sums < epsilon] = 1
    row_norm = sparse.diags(1/row_sums.A.ravel())
    return row_norm@prob_sim