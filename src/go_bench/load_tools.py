import os
from collections import defaultdict

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix, dok_matrix, lil_matrix

#from Bio.UniProt.GOA import GAF20FIELDS
GAF20FIELDS = ['DB', 'DB_Object_ID', 'DB_Object_Symbol', 
            'Qualifier', 'GO_ID', 'DB:Reference', 'Evidence', 
            'With', 'Aspect', 'DB_Object_Name', 'Synonym', 
            'DB_Object_Type', 'Taxon_ID', 'Date', 'Assigned_By', 
            'Annotation_Extension', 'Gene_Product_Form_ID']

def to_data_num(year, month, day):
    return int(year)*10000+int(month)*100+int(day)

#The annotations dict maps protein ids to a list of (GO annotation, evidence_code) tuples. 

def read_table_annotations(tab_path, protein_annotations_dict, code):
    for df in pd.read_table(tab_path, sep='\t', chunksize=4000, 
                usecols=["Entry", "Gene ontology IDs"], keep_default_na=False, dtype=str):
        for row in df.itertuples():
            name = row[1]
            if(len(name) <= 2):
                print(name)
            annotations = row[2].split("; ")
            annotations = [(annotation, code) for annotation in annotations]
            if(row[1] in protein_annotations_dict):
                protein_annotations_dict[row[1]].extend(annotations)
            else:
                protein_annotations_dict[row[1]] = list(annotations)


#Load all protein annotations based on select annotation codes, from a specific time-range. 
def load_protein_annotations(goa_path, annotation_codes, min_date=None, max_date=None):
    annotation_codes = set(annotation_codes)
    annot_dict = defaultdict(set)
    df_iter = pd.read_csv(goa_path, dtype=str,
                          sep='\t',
                          comment='!',
                          names=GAF20FIELDS,
                          chunksize=int(1e6))    
    if(min_date is None and max_date is None):
        for zdf in df_iter:
            # For now, remove all with a qualifier
            zdf = zdf[zdf.Evidence.isin(annotation_codes) & zdf.Qualifier.isnull()]
            for tup in zdf.itertuples():
                annot_dict[tup[2]].add(tup[5])
    else:
        for zdf in df_iter:
            # For now, remove all with a qualifier
            dates = zdf.Date.astype(int)
            
            zdf = zdf[zdf.Evidence.isin(annotation_codes) & zdf.Qualifier.isnull() & dates.ge(min_date) & dates.le(max_date)]
            for tup in zdf.itertuples():
                annot_dict[tup[2]].add(tup[5])
    return annot_dict

#Data management methods and classes
def load_GO_tsv_file(path):
    prot_dict = {}
    df_iter = pd.read_csv(path, dtype=str,
                          sep='\t',
                          header=0,
                          chunksize=int(1e4))
    for zdf in df_iter:
        for tup in zdf.itertuples():
            prot_id = tup[1]
            annotations = tup[2]
            prot_dict[prot_id] = list(annotations.split(", ")) if isinstance(annotations, str) else []
    return prot_dict

#Convert a protein annotation dict to a sparse matrix. 
def convert_to_sparse_matrix(protein_annotation_dict, term_list, prot_id_list):
    term_col_mappings = {term:i for i, term in enumerate(term_list)}
    prot_row_mappings = {prot:i for i, prot in enumerate(prot_id_list)}

    labels = lil_matrix((len(prot_id_list), len(term_list)), dtype=np.int8)

    for place, prot_id in enumerate(prot_id_list):
        for go_id in protein_annotation_dict[prot_id]:
            if(go_id in term_col_mappings):
                labels[place, term_col_mappings[go_id]] = 1
    labels = labels.tocsr()
    return labels

def read_sparse(fn, prot_rows, GO_cols):
    prm = {prot:i for i, prot in enumerate(prot_rows)}
    tcm = {term:i for i, term in enumerate(GO_cols)}
    sparse_probs = dok_matrix((len(prot_rows), len(GO_cols)))
    df = pd.read_csv(fn)
    for (i, prot, go_id, prob) in df.itertuples():
        if(prot in prm and go_id in tcm):
            sparse_probs[prm[prot], tcm[go_id]] = prob
    return csr_matrix(sparse_probs)
