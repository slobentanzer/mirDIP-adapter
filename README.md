# BioCypher mirDIP adapter
This is an adapter for the mirDIP dataset, which is a collection of miRNA-target
interactions. The dataset is available 
[here](https://ophid.utoronto.ca/mirDIP/download.jsp). The adapter is based on
version 5.2 of the dataset, which it assumes to be located in the
`data/mirDIP_Bidirectional_search_v_5_2` folder.

## Usage
To create a standalone BioCypher database of the mirDIP dataset, run the 
`create_mirDIP.py` script. This will use the adapter in `mirDIP_adapter.py` to
create protein and miRNA nodes and relationships between them, ready to load
into a Neo4j DBMS in the `biocypher-out` folder. To integrate the mirDIP
interaction data with other data in an extended BioCypher database, follow the
instructions for adapter usage on
[https://biocypher.org](https://biocypher.org).

