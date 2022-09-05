# GO Bench
go_bench is a python package for generating Gene Ontology datasets from the SwissProt and GOA databases. It contains functions for downloading and processing protein annotations to generate training, validation, and testing sets from predetermined random splits. 

General usage, with download links for required datasets, is shown in gen_datasets.ipynb. To replicate results of [gobench.org](gobench.org), it is important to use the same train/validation/test split used the internal server, instead of the script included in this package. These are included in the data_split directory. 

# Install
go_bench can be installed from git with `pip install git+https://github.com/amdson/GO_pipeline.git`. It requires python >= 3.6, and lists numpy, pandas, scipy, and goatools as dependencies. 

# Data
Most go_bench functions parse or manipulate gene ontology data. Raw data can be downloaded from SwissProt, the GOA database, and the Gene Ontology, with links for all requirements included in the go_bench_usage.ipynb notebook. Reproducing the datasets from [gobench.org](gobench.org) requires an additional set of train/validation/test splits included in the data_splits directory. 

# Contact
Please send questions or suggestions to amdickson (at) berkeley.edu. Also see our lab website at https://biomechanics.berkeley.edu/. 

