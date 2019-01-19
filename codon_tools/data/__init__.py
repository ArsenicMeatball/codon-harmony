from collections import namedtuple
from Bio import Restriction


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


GCParams = namedtuple("GCParams", "name window_size low high")
GC_content = [
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
