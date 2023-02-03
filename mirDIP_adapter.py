from enum import Enum, auto
import hashlib
from itertools import chain
import pickle
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

    GENE_SYMBOL = 0
    MICRORNA = 1
    RANK = 2
    SCORE = 3
    SCORE_CLASS = 7
    SOURCE = 4
    ORIGINAL_SOURCE_GENE_SYMBOL = 5
    ORIGINAL_SOURCE_MICRORNA = 6

    GENE_UNIPROT_ID = auto()  # mapping though pypath


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
        clear_cache: bool = False,
    ):
        self._set_types_and_fields(node_types, edge_types, edge_fields)

        self.data_source = "mirDIP"
        self.data_version = "5.2"
        self.data_licence = "free to use, copy, and modify for academic and non-commercial purposes"

        self.test_mode = test_mode
        self.clear_cache = clear_cache

    def get_nodes(self):
        """
        Returns a generator of node tuples for node types specified in the
        adapter constructor.
        """

        if mirDIPAdapterNodeType.PROTEIN in self.node_types:

            logger.debug("Generating gene nodes.")

            for gene_symbol, uniprots in self.symbol_to_uniprot.items():

                props = {
                    "gene_symbol": gene_symbol,
                    "source": self.data_source,
                    "version": self.data_version,
                    "licence": self.data_licence,
                }

                for id in uniprots:
                    yield (f"uniprot:{id}", "protein", props)

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

    def get_edges(self, batch):
        """
        Returns a generator of edge tuples for edge types specified in the
        adapter constructor.
        """

        logger.debug("Generating miRNA-gene edges.")

        for row in batch.iter_rows():
            gene_symbol = row[mirDIPAdapterMirnaGeneEdgeField.GENE_SYMBOL.value]
            mir_name = row[mirDIPAdapterMirnaGeneEdgeField.MICRORNA.value]

            label = (
                "mirna_protein_interaction_"
                + row[mirDIPAdapterMirnaGeneEdgeField.SCORE_CLASS.value]
            )

            props = {
                "source": row[mirDIPAdapterMirnaGeneEdgeField.SOURCE.value],
                "version": f"{self.data_source} version {self.data_version}",
                "licence": self.data_licence,
            }

            if mirDIPAdapterMirnaGeneEdgeField.RANK in self.edge_fields:
                props["rank"] = row[mirDIPAdapterMirnaGeneEdgeField.RANK.value]

            if mirDIPAdapterMirnaGeneEdgeField.SCORE in self.edge_fields:
                props["score"] = row[
                    mirDIPAdapterMirnaGeneEdgeField.SCORE.value
                ]

            if mirDIPAdapterMirnaGeneEdgeField.SCORE_CLASS in self.edge_fields:
                props["score_class"] = row[
                    mirDIPAdapterMirnaGeneEdgeField.SCORE_CLASS.value
                ]

            uniprots = self.symbol_to_uniprot.get(gene_symbol, [])

            # create md5 hash of row to use as edge id
            _id = hashlib.md5(str(row).encode()).hexdigest()

            for uniprot in uniprots:

                yield (
                    _id,
                    mir_name,
                    f"uniprot:{uniprot}",
                    label,
                    props,
                )

    def get_edge_batches(self, batch_size: int = int(1e6)):
        return self.data.iter_slices(n_rows=batch_size)

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
            self.edge_fields = [
                field for field in mirDIPAdapterMirnaGeneEdgeField
            ]

    def read_data(self):
        """
        Read mirDIP data from the data directory.
        """

        logger.info("Reading mirDIP data from disk.")

        path = os.path.join("data", "mirDIP_Bidirectional_search_v_5_2")

        # column names: load README.txt using polars
        # each column name is on a separate line, skip the first line
        self.columns = pl.read_csv(os.path.join(path, "README.txt"))

        # first column to list
        self.columns = self.columns[
            "Columns for file: mirDIP_Bidirectional_search_v.5.txt"
        ].to_list()

        # read data from mirDIP_Bidirectional_search_v.5.txt, using columns as
        # column names
        if not self.test_mode:

            self.data = pl.read_csv(
                os.path.join(path, "mirDIP_Bidirectional_search_v.5.txt"),
                has_header=False,
                new_columns=self.columns,
            )

        else:

            self.data = pl.read_csv(
                os.path.join(path, "mirDIP_Bidirectional_search_v.5.txt"),
                has_header=False,
                new_columns=self.columns,
                n_rows=20_000_000,
            )

        # show dataframe description
        logger.info("Printing mirDIP dataframe description.")
        print(self.data.describe())

        # # show data where GENE_SYMBOL and GENE_SYMBOL_ORI are not the same
        # print(
        #     self.data.filter(pl.col("GENE_SYMBOL") != pl.col("GENE_SYMBOL_ORI"))
        # )

        # # show data where MICRORNA and MICRORNA_ORI are not the same
        # print(self.data.filter(pl.col("MICRORNA") != pl.col("MICRORNA_ORI")))

        self._translate_gene_symbols()

    def _translate_gene_symbols(self):
        """
        Translate gene symbols to uniprot ids using pypath.
        """

        # load from pickle if available
        if (
            os.path.exists("data/symbol_to_uniprot.pickle")
            and not self.clear_cache
        ):

            logger.info("Loading symbol_to_uniprot from pickle.")

            with open("data/symbol_to_uniprot.pickle", "rb") as f:
                self.symbol_to_uniprot = pickle.load(f)

            with open("data/unmapped_gene_symbols.pickle", "rb") as f:
                self.unmapped_gene_symbols = pickle.load(f)

            return

        # else use pypath
        logger.info("Translating gene symbols to uniprot ids using pypath.")

        self.unmapped_gene_symbols = set()
        self.symbol_to_uniprot = {}

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

            self.symbol_to_uniprot[gene_symbol] = uniprot_id

        # pickle to file in data/
        with open("data/symbol_to_uniprot.pickle", "wb") as f:
            pickle.dump(self.symbol_to_uniprot, f)

        # pickle unmapped gene symbols to file in data/
        with open("data/unmapped_gene_symbols.pickle", "wb") as f:
            pickle.dump(self.unmapped_gene_symbols, f)
