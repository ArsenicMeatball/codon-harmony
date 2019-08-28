import json
import python_codon_tables as pct
from collections import namedtuple
from Bio import Entrez, Restriction
from ..util import logging

# Species-specifc data can be found on the Codon Usage Database
# (http://www.kazusa.or.jp/codon) using the NCBI Taxonomy database
# (http://www.ncbi.nlm.nih.gov/taxonomy) id (e.g. 413997)
# or the organism's Latin name (e.g. Escherichia coli B).
# Mapping species names to Taxonomy IDs can be done here:
_tax_id_url = "https://www.ncbi.nlm.nih.gov/Taxonomy/TaxIdentifier/tax_identifier.cgi"

logger = logging.getLogger(__name__)
Entrez.email = "login@anonymous.com"


def _tax_id_from_species(species, exc=None):
    """Map the name of a species from a string to the NCBI taxonomy ID and
    return it.

    Args:
        species (str): Name of the species to map.

    Raises:
        ValueError: If the NCBI taxonomy ID cannot be determined, raise a ``ValueError``
        informing the user and directing them to the NCBI Taxonomy Browser.

    Returns:
        int: The NCBI taxonomy ID for the supplied species.
    """
    search_species = species.replace(" ", "+").replace("_", "+").strip()
    handle = Entrez.esearch(term=search_species, db="taxonomy")
    record = Entrez.read(handle)
    try:
        taxid = int(record["IdList"].pop())
    except IndexError:
        raise ValueError(
            '"{}" is not a valid host id. '.format(species)
            + "Supported hosts (Latin and NCBI taxonomy IDs) can be found at "
            + _tax_id_url
        ) from exc

    logger.info(
        'Mapped host name "{}" to NCBI taxonomy ID "{}".'.format(species, taxid)
    )
    return taxid


def codon_tables(taxid, table_path=None):
    """Download the codon use table for the given species and return it as
    a dictionary.

    Returns:
        int: The NCBI taxonomy ID for the supplied species.

    Args:
        taxid (int): NCBI taxonomy ID for the desrired species.
        table_path (str): Defaults to None. Path to a JSON-formatted file representing the
            codon usage to consider. If None, the table is fetched from the internet.

    Raises:
        ValueError: If the NCBI taxonomy ID is not associated with a codon
        usage table, raise a ``ValueError`` informing the user and directing
        them to the NCBI Taxonomy Browser.

    Returns:
        dict{str, float}: A dictionary with codons as keys and the frequency
        that the codon is used to encode its amino acid as values.
    """
    if table_path is None:
        try:
            taxid = int(taxid)
        except ValueError as exc:
            taxid = _tax_id_from_species(taxid, exc)
        logger.info("Downloading host table for NCBI taxonomy ID {}...".format(taxid))
        codon_table_by_aa = pct.download_codons_table(taxid)
    else:
        # load table from disk -- JSON format
        with open(table_path, "r") as table:
            codon_table_by_aa = json.load(table)

    return_dict = {}
    for _, codon_dict in codon_table_by_aa.items():
        for codon, frequency in codon_dict.items():
            return_dict[codon] = frequency
    if not return_dict:
        raise ValueError(
            '"{}" is not a valid host id. '.format(taxid)
            + "Supported hosts (Latin and NCBI taxonomy IDs) can be found at "
            + _tax_id_url
        )

    return return_dict


def RestrictionEnzymes(restriction_enzymes):
    """Create a RestrictionBatch instance to search for sites for a supplied
    list of restriction enzymes.

    Args:
        restriction_enzymes (list[str], optional): List of restriction
            enzymes to consider. Defaults to ["NdeI", "XhoI", "HpaI", "PstI",
            "EcoRV", "NcoI", "BamHI"].

    Returns:
        Bio.Restriction.Restriction.RestrictionBatch: RestrictionBatch instance
        configured with the input restriction enzymes.
    """
    return Restriction.RestrictionBatch(
        [Restriction.AllEnzymes.get(enz) for enz in restriction_enzymes]
    )


class GCParams(namedtuple("GCParams", "name window_size low high")):
    """High and low values for GC-content within a specified window size.

    Attributes:
        name (str): Name of the parameter set.
        window_size (int): Number of nucleotides over which the GC content
            will be calculated.
        low (float): The minimum fraction of GC in the window.
        high (float): The maximum fraction of GC in the window.
    """


GC_content = [
    GCParams("Custom", 20, 0.15, 0.60),
    GCParams("IDT", 20, 0.15, 0.90),
    GCParams("twist", 50, 0.15, 0.80),
    GCParams("IDT_long", 100, 0.28, 0.68),
    GCParams("twist_long", "x3", 0.3, 0.65),
]

RibosomeBindingSites = {
    "rbs_0": "GGGGG",
    "rbs_1": "GGGGA",
    "rbs_2": "GGGAG",
    "rbs_3": "GGGAA",
    "rbs_4": "GGAGG",
    "rbs_5": "GGAGA",
    "rbs_6": "GGAAG",
    "rbs_7": "GGAAA",
    "rbs_8": "GAGGG",
    "rbs_9": "GAGGA",
    "rbs_10": "GAGAG",
    "rbs_11": "GAGAA",
    "rbs_12": "GAAGG",
    "rbs_13": "GAAGA",
    "rbs_14": "GAAAG",
    "rbs_15": "GAAAA",
    "rbs_16": "AGGGG",
    "rbs_17": "AGGGA",
    "rbs_18": "AGGAG",
    "rbs_19": "AGGAA",
    "rbs_20": "AGAGG",
    "rbs_21": "AGAGA",
    "rbs_22": "AGAAG",
    "rbs_23": "AGAAA",
    "rbs_24": "AAGGG",
    "rbs_25": "AAGGA",
    "rbs_26": "AAGAG",
    "rbs_27": "AAGAA",
    "rbs_28": "AAAGG",
    "rbs_29": "AAAGA",
    "rbs_30": "AAAAG",
    "rbs_31": "AAAAA",
}
