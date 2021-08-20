import pandas as pd
import random

def nested_len(cluster_list):
    return sum(len(x) for x in cluster_list)

def load_clusters(clusters_path, protein_ids):
    cluster_list = []
    protein_count = 0
    reader = pd.read_csv(clusters_path, sep='\t', lineterminator='\n', chunksize = 10000)
    chunk_count = 0
    for chunk in reader:
        chunk_count += 1
        column = chunk["Cluster members"]
        for id_list in column:
            unfiltered_cluster = id_list.split("; ")
            cluster = []
            for prot_id in unfiltered_cluster:
                if(prot_id in protein_ids):
                    cluster.append(prot_id)
            if(len(cluster) > 0):
                cluster_list.append(cluster)
    return cluster_list

def build_random_split(cluster_list, train_fraction=0.9):
    protein_count = nested_len(cluster_list)
    random.shuffle(cluster_list)
    i = 0
    training_size = 0
    training_set = []
    while (training_size < train_fraction*protein_count):
        cluster = cluster_list[i]
        training_set.extend(cluster)
        training_size += len(cluster)
        i += 1
    testing_set = []
    while (i < len(cluster_list)):
        cluster = cluster_list[i]
        testing_set.extend(cluster)
        i += 1
    return training_set, testing_set

def build_random_cluster_split(cluster_list, train_fraction=0.9):
    protein_count = nested_len(cluster_list)
    random.shuffle(cluster_list)
    i = 0
    training_size = 0
    training_set = []
    while (training_size < train_fraction*protein_count):
        cluster = cluster_list[i]
        training_set.append(cluster)
        training_size += len(cluster)
        i += 1
    testing_set = []
    while (i < len(cluster_list)):
        cluster = cluster_list[i]
        testing_set.append(cluster)
        i += 1
    return training_set, testing_set