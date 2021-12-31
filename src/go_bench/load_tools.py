import os
from collections import defaultdict
import gzip

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix, dok_matrix, lil_matrix
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

#from Bio.UniProt.GOA import GAF20FIELDS
GAF20FIELDS = ['DB', 'DB_Object_ID', 'DB_Object_Symbol', 
            'Qualifier', 'GO_ID', 'DB:Reference', 'Evidence', 
            'With', 'Aspect', 'DB_Object_Name', 'Synonym', 
            'DB_Object_Type', 'Taxon_ID', 'Date', 'Assigned_By', 
            'Annotation_Extension', 'Gene_Product_Form_ID']

def to_data_num(year, month, day):
    return int(year)*10000+int(month)*100+int(day)

#The annotations dict maps protein ids to a list of (GO annotation, evidence_code) tuples. 

neg_qualifiers = {'NOT|colocalizes_with', 'NOT|acts_upstream_of_or_within', 
                  'NOT|acts_upstream_of', 'NOT|acts_upstream_of_or_within_negative_effect', 
                  'NOT|enables', 'NOT|part_of', 'NOT|is_active_in', 'NOT|located_in', 
                  'NOT|contributes_to', 'NOT|involved_in', 'acts_upstream_of_or_within_negative_effect', 
                  'acts_upstream_of_negative_effect'}
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
            zdf = zdf[zdf.Evidence.isin(annotation_codes) & ~zdf.Qualifier.isin(neg_qualifiers)]
            for tup in zdf.itertuples():
                annot_dict[tup[2]].add(tup[5])
    else:
        for zdf in df_iter:
            # For now, remove all with a qualifier
            dates = zdf.Date.astype(int)
            zdf = zdf[zdf.Evidence.isin(annotation_codes) & ~zdf.Qualifier.isin(neg_qualifiers) & dates.ge(min_date) & dates.le(max_date)]
            for tup in zdf.itertuples():
                annot_dict[tup[2]].add(tup[5])
    return annot_dict

#loads proteins sequences from a tab or fasta file and filters by a whitelist. Returns sequences and 
#"prot_row_mappings" giving the entry in sequences that corresponds with each protein id. 
def load_protein_sequences(path, prot_whitelist=None):
    prot_whitelist = set(prot_whitelist)

    sequences = []
    prot_ids = []
    if(path[-3:] == "tab" or path[-6:] == "tab.gz"):
        
        df_iter = pd.read_csv(path, dtype=str,
                                sep='\t',
                                header=0,
                                chunksize=int(1e5)) 
        print("Loading sequences")
        count = 0
        for df in df_iter:
            print(f"{count*1e5} rows read")
            if(prot_whitelist):
                df = df[df["Entry"].isin(prot_whitelist)]
            sequences.extend([s.lower() for s in df["Sequence"]])
            prot_ids.extend(df["Entry"])
            count += 1
    elif(path[-8:] == "fasta.gz"):
        with gzip.open(path, 'rt') as fp:
            seqs = SeqIO.parse(fp, 'fasta')
            count = 0
            for rec in seqs:
                seq_id = rec.id
                if('|' in seq_id):
                    seq_id = rec.id.split('|')[1]
                seq = rec.seq
                if(seq_id in prot_whitelist):
                    sequences.append(str(seq.lower()))
                    prot_ids.append(seq_id)
                if(count % 100000 == 0):
                    print(f"{count} proteins processed")
                count += 1
    elif(path[-5:] == "fasta"):
        with open(path, 'r') as fp:
            seqs = SeqIO.parse(fp, 'fasta')
            count = 0
            for rec in seqs:
                seq_id = rec.id
                if('|' in seq_id):
                    seq_id = rec.id.split('|')[1]
                seq = rec.seq
                if(prot_whitelist is None or seq_id in prot_whitelist):
                    sequences.append(str(seq.lower()))
                    prot_ids.append(seq_id)
                if(count % 10000 == 0):
                    print(f"{count} proteins processed")
                count += 1
    return sequences, prot_ids

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

def load_tsv_as_sparse(path):
    annotation_dict = load_GO_tsv_file(path)


def read_sparse(fn, prot_rows, GO_cols):
    prm = {prot:i for i, prot in enumerate(prot_rows)}
    tcm = {term:i for i, term in enumerate(GO_cols)}
    sparse_probs = dok_matrix((len(prot_rows), len(GO_cols)))
    df = pd.read_csv(fn)
    for (i, prot, go_id, prob) in df.itertuples():
        if(prot in prm and go_id in tcm):
            sparse_probs[prm[prot], tcm[go_id]] = prob
    return csr_matrix(sparse_probs)

def write_sparse(fn, preds, prot_rows, GO_cols, go, min_certainty):
    with open(fn, mode='w') as f:
        f.write("g,t,s\n")
        for row, col in zip(*preds.nonzero()):
            prot_id = prot_rows[row]
            go_id = GO_cols[col]
            val = preds[row, col]
            if(val > min_certainty and go_id in go):
                f.write(f"{prot_id},{go_id},{val}\n")
