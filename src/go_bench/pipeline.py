import json
import logging
import os
import pickle

import pandas as pd
from go_bench.processing import (enforce_count, enforce_threshold,
                                       filter_dict, get_counts_dict, get_namespace_terms,
                                       invert_protein_annotation_dict,
                                       propogate_annotations)
from go_bench.load_tools import load_protein_annotations
from go_bench.utils import namespaces, experimental_codes

def construct_tsv(path, prot_dict, prot_ids, term_set):
    print(path, len(prot_dict), len(prot_ids), len(term_set))
    columns = ["Uniprot ID", "GO Associations"]
    with open(path, mode='w') as f:
        f.write("\t".join(columns) + "\n")
        for prot_id in prot_ids:
            if(prot_id in prot_dict):
                terms = term_set.intersection(prot_dict[prot_id])
                if(len(terms) > 0):
                    f.write("{}\t{}\n".format(prot_id, ", ".join(terms)))

set_types = ("training", "validation", "testing")
def pipeline(goa_path, split_path, save_dir, godag, codes=experimental_codes, namespaces=namespaces, set_types=set_types,
                propogate_terms=True, min_date=None, max_date=None,
                filter_type=('min_samples', 100), filter_testing=False):
    """
    Inputs:
        goa_path: Path to a GOA formatted tsv file. 
        split_path: Path to a directory containing training_ids.txt, validation_ids.txt, and testing_ids.txt. 
        Each file should contain a newline separated list of UniProt protein identifiers representing those IDs allowed in
        each portion of the train test split. 
        save_dir: Path to the directory where results should be stored. 
        godag: A goatools.obo_parser.GODag object, loaded with prefered version of the gene ontology. Example link in ipynb. 
        experimental_codes: A set of the ontology evidence quality identifiers which should be allowed into the model. 
        namespaces: List of namespaces that should be generated in output dataset. 
        Any set of ("molecular_function", "biological_process", "cellular_component") is valid. 
        set_types: List of set types to be included in output. Any set of ("training", "validation", "testing") is valid. 
        propogate_terms: Control whether protein gene ontology term annotations are augmented by inluding ancestors in the 
        gene ontology tree. 
        min_date: Min date of annotations to include from GOA file. Date must be formatted by to_data_num from load_tools.py. 
        max_date: Max date of annotations to include from GOA file. 
    """
    print("loading annotations")
    prot_dict = load_protein_annotations(goa_path, codes, min_date=min_date, max_date=max_date) 
    print("initial len", sum(len(x) for x in prot_dict.values()))
    print("filtering annotations to godag")
    prot_dict = filter_dict(prot_dict, godag)
    print("filtered len", sum(len(x) for x in prot_dict.values()))
    #Read in terms
    term_list = [term for term in godag]
    #Count occurences of each GO term
    print("propogating annotations")
    if(propogate_terms):
        propogate_annotations(prot_dict, term_list, godag)
    print("propogated len", sum(len(x) for x in prot_dict.values()))
    print("getting counts")
    annotation_dict = invert_protein_annotation_dict(prot_dict) #Maps GO IDs to lists of protein IDs. 
    annotation_counts_dict = get_counts_dict(annotation_dict)
    
    #Filter terms for each namespace
    filter_type, f_k = filter_type
    if(filter_type == "min_samples"):
        filter_method = lambda filter_set: enforce_threshold(annotation_counts_dict, filter_set, f_k)
    elif(filter_type == "top_k"):
        filter_method = lambda filter_set: enforce_count(annotation_counts_dict, filter_set, f_k)

    analysis_content = {}
    print("generating datasets")
    for namespace in namespaces:
        namespace_terms = get_namespace_terms(godag, namespace)

        print(f"{len(namespace_terms)} terms in {namespace}")
        print("filtering namespace:", namespace)

        filtered_list = filter_method(namespace_terms)
        print("filtered_list length", len(filtered_list))

        #Term counts for return        
        namespace_term_counts = [annotation_counts_dict[term] for term in filtered_list]
        namespace_df = pd.DataFrame(list(zip(filtered_list, namespace_term_counts)), columns=["GO_term", "count"])
        analysis_content[namespace] = namespace_df
        
        json_path = f"{save_dir}/{namespace}_terms.json"
        with open(json_path, "w") as f:
            json.dump(filtered_list, f)
            
        for prot_set_type in set_types:
            if(prot_set_type == "testing" and not filter_testing):
                write_list = enforce_threshold(annotation_counts_dict, namespace_terms, 0) 
                print("testing list length", len(write_list))
                json_path = f"{save_dir}/testing_{namespace}_terms.json"
                with open(json_path, "w") as f:
                    json.dump(list(write_list), f)
            else:
                write_list = filtered_list
            with open(f"{split_path}/{prot_set_type}_ids.txt", "r") as f:
                prot_ids = set([x[:-1] for x in f.readlines()])
            print("prot_ids:", len(prot_ids))
            print("saving results")
            path = f"{save_dir}/{prot_set_type}_{namespace}_annotations.tsv"
            print(f"saving to {path}")
            construct_tsv(path, prot_dict, prot_ids, set(write_list))
    return analysis_content
