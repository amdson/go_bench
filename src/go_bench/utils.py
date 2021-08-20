namespaces = ("molecular_function", "biological_process", "cellular_component")

evidence_codes = ('ISO', 'IGI', 'ISA', 'IGC', 'IEP', 
                  'RCA', 'HDA', 'HGI', 'IKR', 'TAS', 'HEP', 
                  'ND', 'IBA', 'IMP', 'EXP', 'IDA', 'IC', 'ISM', 
                  'ISS', 'NAS', 'IRD', 'IEA', 'IPI', 'HMP')


def extract_code_acronyms(code_description):
    return tuple([x[1] for x in code_description.values()])

def eca(code_description):
    return tuple([x[1] for x in code_description.values()])

experimental_descr = {"experiment":("ECO_0000269", 'EXP'), 
                        "direct_assay":("ECO_0000314", 'IDA'), 
                        "physical_interaction":("ECO_0000353", 'IPI'), 
                        "mutant_phenotype":("ECO_0000315", 'IMP'),
                        "genetic_interaction":("ECO_0000316", 'IGI'),
                        "expression_pattern":("ECO_0000270", 'IEP')}
experimental_codes = eca(experimental_descr)

high_throughput_descr = {"high_throughput_experiment": ("ECO_0006056", 'HTP'),
                        "high_throughput_direct_assay": ("ECO_0007005", 'HDA'),
                        "high_throughput_mutant_phenotype":("ECO_0007001", 'HMP'),
                        "high_throughput_genetic_interaction": ("ECO_0007003", 'HGI'),
                        "high_throughput_expression_pattern":("ECO_0007007", 'HEP')
                        }
high_throughput_codes = eca(high_throughput_descr)

phylogenetic_descr = {"biological_aspect_of_ancestor":("ECO_0000318", 'IBA'),
                        "biological_aspect_of_descendant":("ECO_0000319", 'IBD'),
                        "key_residues":("ECO_0000320", 'IKR'),
                        "rapid_divergence":("ECO_0000320", 'IRD')
}
phylogenetic_codes = eca(phylogenetic_descr)

reviewed_computational_descr = {"sequence_or_structural_similarity":("ECO_0000250", 'ISS'),
                        "sequence_orthology":("ECO_0000266", 'ISO'),
                        "sequence_alignment":("ECO_0000247", 'ISA'),
                        "sequence_model":("ECO_0000255", 'ISM'),
                        "genomic_context":("ECO_0000317", 'IGC'),
                        "reviewed_computational_analysis":("ECO_0000245", 'RCA')
}
reviewed_computational_codes = eca(reviewed_computational_descr)

unreviewed_computational_descr = {
    "inferred from electronic annotation": ('IEA'),
}
unreviewed_computational_codes = eca(unreviewed_computational_descr)

reviewed_misc_descr = {
    "tracable author statement": ("_", 'TAS'),
    "inferred by curator": ("_", 'IC')
}
reviewed_misc_codes = eca(reviewed_misc_descr)

unreviewed_misc_descr = {
    "non-tracable author statement": ("_", 'NAS'),
    "no biological data available": ("_", 'ND')
}
unreviewed_misc_codes = eca(unreviewed_misc_descr)

exp_group = experimental_codes + reviewed_misc_codes
cafa_group = exp_group
exp_phy_ht_group = exp_group + high_throughput_codes + phylogenetic_codes
non_iea_group = exp_phy_ht_group + reviewed_computational_codes
all_group = evidence_codes



def list_ancestors(term, godag):
    """Generator yielding all ancestors of a GO term."""
    ancestors = list(godag[term]._parents)
    seen = set()
    while ancestors:
        ancestor = ancestors.pop()
        if(ancestor in seen):
            continue
        seen.add(ancestor)
        ancestors.extend(godag[ancestor]._parents)
        yield ancestor