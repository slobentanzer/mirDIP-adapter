from enum import Enum, auto
from itertools import chain
import polars as pl
import os
from pypath.utils import mapping
from tqdm import tqdm
from biocypher._logger import logger

logger.debug(f"Loading module {__name__}.")


class mirDIPAdapterNodeType(Enum):
    """
    Define types of nodes the adapter can provide.
    """

    PROTEIN = auto()
    # pypath mapping to ensg fails, use uniprot in the meantime
    MICRORNA = auto()


class mirDIPAdapterEdgeType(Enum):
    """
    Enum for the types of the protein adapter.
    """

    MIRNA_PROTEIN_INTERACTION = auto()


class mirDIPAdapterMirnaGeneEdgeField(Enum):
    """
    Define possible fields the adapter can provide for protein-protein edges.
    """

    GENE_SYMBOL = "GENE_SYMBOL"
    GENE_UNIPROT_ID = "GENE_UNIPROT_ID"  # mapping though pypath
    MICRORNA = "MICRORNA"
    RANK = "RANK"
    SCORE = "SCORE"
    SOURCE = "SOURCE_NAME"
    ORIGINAL_SOURCE_GENE_SYMBOL = "GENE_SYMBOL_ORI"
    ORIGINAL_SOURCE_MICRORNA = "MICRORNA_ORI"
    SCORE_CLASS = "SCORE_CLASS"


class mirDIPAdapter:
    """
    BioCypher adapter for mirDIP. Generates miRNA and protein nodes and
    miRNA-protein edges for creating a knowledge graph.

    Args:
        fields (list): List of fields to include in the node.
    """

    def __init__(
        self,
        node_types: str = None,
        edge_types: str = None,
        edge_fields: str = None,
        test_mode: bool = False,
    ):
        self._set_types_and_fields(node_types, edge_types, edge_fields)

        self.data_source = "mirDIP"
        self.data_version = "5.2"
        self.data_licence = "free to use, copy, and modify for academic and non-commercial purposes"

        self.test_mode = test_mode

    def get_nodes(self):
        """
        Returns a generator of node tuples for node types specified in the
        adapter constructor.
        """

        if mirDIPAdapterNodeType.PROTEIN in self.node_types:

            logger.debug("Generating gene nodes.")

            self.unmapped_gene_symbols = set()

            gene_symbols = set(self.data.get_column("GENE_SYMBOL").to_list())

            for gene_symbol in tqdm(gene_symbols):

                # get ensg id from pypath
                uniprot_id = mapping.map_name(
                    name=gene_symbol,
                    id_type="genesymbol",
                    target_id_type="uniprot",
                    ncbi_tax_id=9606,
                )

                if not uniprot_id:

                    # get original gene symbol from self.data
                    original_gene_symbol = (
                        self.data.filter(pl.col("GENE_SYMBOL") == gene_symbol)
                        .get_column("GENE_SYMBOL_ORI")
                        .to_list()[0]
                    )

                    if original_gene_symbol == gene_symbol:
                        self.unmapped_gene_symbols.add(gene_symbol)
                        continue

                    uniprot_id = mapping.map_name(
                        name=original_gene_symbol,
                        id_type="genesymbol",
                        target_id_type="uniprot",
                        ncbi_tax_id=9606,
                    )

                if not uniprot_id:
                    self.unmapped_gene_symbols.add(gene_symbol)
                    continue

                props = {
                    "gene_symbol": gene_symbol,
                    "source": self.data_source,
                    "version": self.data_version,
                    "licence": self.data_licence,
                }

                for id in uniprot_id:
                    yield (id, "protein", props)

            print(self.unmapped_gene_symbols)
            print(len(self.unmapped_gene_symbols))

        if mirDIPAdapterNodeType.MICRORNA in self.node_types:

            logger.debug("Generating miRNA nodes.")

            mirs = set(self.data.get_column("MICRORNA").to_list())

            for mir in mirs:

                props = {
                    "source": self.data_source,
                    "version": self.data_version,
                    "licence": self.data_licence,
                }

                yield (mir, "mirna", props)

    def get_edges(self):
        """
        Returns a generator of edge tuples for edge types specified in the
        adapter constructor.
        """

        logger.debug("Generating miRNA-gene edges.")

        pass

    def _set_types_and_fields(self, node_types, edge_types, edge_fields):
        if node_types:
            self.node_types = node_types
        else:
            self.node_types = [type for type in mirDIPAdapterNodeType]

        if edge_types:
            self.edge_types = edge_types
        else:
            self.edge_types = [type for type in mirDIPAdapterEdgeType]

        if edge_fields:
            self.edge_fields = edge_fields
        else:
            self.edge_fields = [field for field in chain()]

    def _read_data(self):
        """
        Read mirDIP data from the data directory.
        """

        logger.info("Reading mirDIP data from disk.")

        path = os.path.join("data", "mirDIP_Bidirectional_search_v_5_2")

        # column names: load README.txt using polars
        # each column name is on a separate line, skip the first line
        columns = pl.read_csv(os.path.join(path, "README.txt"))

        # first column to list
        columns = columns[
            "Columns for file: mirDIP_Bidirectional_search_v.5.txt"
        ].to_list()

        # read data from mirDIP_Bidirectional_search_v.5.txt, using columns as
        # column names
        if not self.test_mode:

            self.data = pl.read_csv(
                os.path.join(path, "mirDIP_Bidirectional_search_v.5.txt"),
                has_header=False,
                new_columns=columns,
            )

        else:

            self.data = pl.read_csv(
                os.path.join(path, "mirDIP_Bidirectional_search_v.5.txt"),
                has_header=False,
                new_columns=columns,
                n_rows=200,
            )

        # # show data where GENE_SYMBOL and GENE_SYMBOL_ORI are not the same
        # print(
        #     self.data.filter(pl.col("GENE_SYMBOL") != pl.col("GENE_SYMBOL_ORI"))
        # )

        # # show data where MICRORNA and MICRORNA_ORI are not the same
        # print(self.data.filter(pl.col("MICRORNA") != pl.col("MICRORNA_ORI")))
