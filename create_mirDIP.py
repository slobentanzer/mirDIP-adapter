import biocypher
from tqdm import tqdm

from mirDIP_adapter import (
    mirDIPAdapter,
    mirDIPAdapterNodeType,
    mirDIPAdapterEdgeType,
)

# Instantiate the BioCypher driver
# You can use `config/biocypher_config.yaml` to configure the driver or supply
# settings via parameters below
driver = biocypher.Driver(
    user_schema_config_path="config/schema_config.yaml",
    skip_bad_relationships=True,  # Neo4j admin import option
    skip_duplicate_nodes=True,  # Neo4j admin import option
    wipe=True,  # Neo4j admin import option
    strict_mode=True,
)

# Take a look at the ontology structure of the KG according to the schema
driver.show_ontology_structure()

# Choose node types to include in the knowledge graph.
# These are defined in the adapter (`adapter.py`).
node_types = [
    mirDIPAdapterNodeType.PROTEIN,
    mirDIPAdapterNodeType.MICRORNA,
]

edge_types = [
    mirDIPAdapterEdgeType.MIRNA_PROTEIN_INTERACTION,
]

# Create a protein adapter instance
adapter = mirDIPAdapter(
    node_types=node_types,
    edge_types=edge_types,
    # we can leave edge fields empty, defaulting to all fields in the adapter
    test_mode=True,
    clear_cache=True,
)

adapter.read_data()


# Create a knowledge graph from the adapter
driver.write_nodes(adapter.get_nodes())

batches = adapter.get_edge_batches()
for batch in tqdm(batches):
    driver.write_edges(adapter.get_edges(batch=batch))

# Write admin import statement
driver.write_import_call()

# Check output
driver.log_duplicates()
driver.log_missing_bl_types()
