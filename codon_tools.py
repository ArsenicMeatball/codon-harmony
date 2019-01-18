#! /usr/bin/env python
import argparse
import random
import re
import sys

import Bio
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.Data import CodonTable
from Bio.Restriction import Analysis
from Bio.SeqUtils import CodonUsage, GC, seq3

from codon_tools.util import Seq, MutableSeq, codon_use, logging, log_levels
from codon_tools.data import GC_content, RibosomeBindingSites, get_restriction_enzymes


# options
parser = argparse.ArgumentParser(
    description="Optimize your AA or DNA sequence to harmonize with a host's condon usage.",
    epilog="2017-12-04, v0.47 (contact yhsia@uw.edu if stuff does not work)",
)
parser.add_argument("--input", type=str, required=True, help="input file with sequence")
parser.add_argument(
    "--host",
    type=str,
    default="413997",
    help='host table code: http://www.kazusa.or.jp/codon/, default is "Escherichia coli B"',
)
parser.add_argument(
    "--host_threshold",
    type=float,
    default="0.10",
    help="lowest codon fraction per AA in the host that is allowed",
)
parser.add_argument(
    "--verbose",
    type=int,
    default=0,
    choices=[0, 1, 2, 3],
    help="verbose output level (0=only result, 1=standard output, 2=extra output 3=debugging)",
)
parser.add_argument(
    "--local_homopolymer_threshold",
    type=int,
    default="4",
    help="number of consecutive NT repeats allowed",
)
parser.add_argument(
    "--cycles",
    type=int,
    default=1000,
    help="max number of cycles to run optimization, 0=unlimited",
)
parser.add_argument(
    "--max_relax",
    type=float,
    default="0.1",
    help="maximum % deviation from host profile",
)

args = parser.parse_args()


logging.basicConfig(level=log_levels[args.verbose])
logger = logging.getLogger(__name__)
random.seed()


def compare_profiles(codons_count, host, relax):
    """Compute the deviation from the expected codon usage based on a host
    codon usage profile.

    Note:
        The `relax` parameter uniformly increases the host codon usage that
        is used to estimate the number of times each codon should appear in
        the sequence. These values are rounded and then iteratively adjusted
        to be consistent with the length of the sequence of interest.
        Increasing this parameter further distorts the codon use distribution
        from the host.

    Args:
        codons_count (dict({str : int})): A dictionary with each codon as
            keys and the number of times it appears in a gene as values.
        host (dict({str : foat})): A dictionary with each codon as keys and the
            frequency of its use in the host organism as values.
        relax (float): The maximum deviation from the host profile to tolerate.

    Returns:
        dict({str : dict({str : int})}): A dictionary with
            each codon as keys, and dictionaries of the difference between
            the observed and expected codon usage.
        float: The number of mutations per residue that are needed to make
        the sequence match the host codon usage.
    """
    logger.info("===== COMPARING PROFILES =====")
    table = {}
    # loop AAs
    for AA, synonymous_codons in CodonUsage.SynonymousCodons.items():
        logger.detail(AA)
        temp_table = {}

        # calculate total usage of codon in input
        tot_usage = sum([codons_count[codon] for codon in synonymous_codons])

        # calculate ideal usage of codon in host
        tot_ideal = 0
        for codon in synonymous_codons:
            ideal_usage_abs = int(round(host[codon] * tot_usage, 0))
            ideal_usage = int(round(host[codon] * relax * tot_usage, 0))
            logger.detail("{0}: {1}".format(codon, ideal_usage))
            tot_ideal += ideal_usage
            temp_table[codon] = {
                "input_count": codons_count[codon],
                "input_perc": codons_count[codon],
                "ideal_usage_abs": ideal_usage_abs,
                "ideal_usage": ideal_usage,
                "host_perc": host[codon],
            }

        # account for rounding issues and relaxation of the host profile
        # by adjusting the most- and least-used codons as necessary

        # if the calculated usage exceeds the actual usage, subtract one from
        # the "ideal_use" of the least-used host codon that appears in the sequence
        while tot_ideal > tot_usage:
            codon = min(
                [c for c, d in temp_table.items() if d["ideal_usage"] > 0],
                key=lambda codon: temp_table[codon]["host_perc"],
            )
            temp_table[codon]["ideal_usage"] -= 1
            tot_ideal -= 1

        # if the calculated usage is less than the actual usage, add one to
        # the "ideal_use" of the most-used host codon
        while tot_ideal < tot_usage:
            codon = max(temp_table, key=lambda codon: temp_table[codon]["host_perc"])
            temp_table[codon]["ideal_usage"] += 1
            tot_ideal += 1

        # populate return table
        for codon, usage in temp_table.items():
            table[codon] = {
                "input_count": usage["input_count"],
                "ideal_usage_abs": usage["ideal_usage_abs"],
                "difference": usage["input_count"] - usage["ideal_usage"],
                "difference_abs": usage["input_count"] - usage["ideal_usage_abs"],
            }

    # calculate difference value
    number_of_residues, diff_total = 0, 0
    for _, synonymous_codons in CodonUsage.SynonymousCodons.items():
        for codon in synonymous_codons:
            number_of_residues += table[codon]["ideal_usage_abs"]
            diff_total += abs(table[codon]["difference_abs"])

    diff_total /= 2  # convert to the number of mutations needed
    diff = diff_total / number_of_residues
    logger.info(
        "CURRENT DIFFERENCE TOTAL: {0} of {1}".format(
            int(diff_total), number_of_residues
        )
    )
    logger.info("CURRENT DIFFERENCE %: {0}".format(diff))

    return table, diff


# returns mutated sequence after given a difference profile
def harmonize_codon_use_with_host(dna_sequence, mutation_profile):
    """Adjust the codon usage in the DNA sequence to be consistent with
    the host profile.

    Args:
        dna_sequence (Bio.Seq.Seq): A read-only representation of
            the DNA sequence.
        mutation_profile (dict({str : dict({str : int})})): A dictionary
            with each codon as keys, and dictionaries of the difference
            between the observed and expected codon usage.

    Returns:
        Bio.Seq.Seq: A read-only representation of the new DNA sequence.
    """
    logger.info("===== OPTIMIZING SEQENCE =====")

    mutable_seq = dna_sequence.tomutable()
    for _, synonymous_codons in CodonUsage.SynonymousCodons.items():
        # get the index of relevant codons in the sequence
        mutation_table = {}
        for codon in synonymous_codons:
            mutation_table[codon] = {
                "difference": mutation_profile[codon]["difference"]
            }

            pos_list = []
            for i in range(0, len(mutable_seq), 3):
                codon_idx = slice(i, i + 3)
                if mutable_seq[codon_idx] == codon:
                    pos_list.append(codon_idx)

            mutation_table[codon]["pos"] = pos_list

        # check if this AA even needs to be adjsuted
        tot_diff = sum(
            abs(mutation_table[codon]["difference"]) for codon in synonymous_codons
        )

        while tot_diff > 0:
            # randomly select a pair of codons to adjust
            codon_to_remove = random.choice(
                [c for c, d in mutation_table.items() if d["difference"] > 0]
            )
            codon_to_add = random.choice(
                [c for c, d in mutation_table.items() if d["difference"] < 0]
            )

            # randomly select the position to update
            codon_idx = random.choice(mutation_table[codon_to_remove]["pos"])

            # remove from sequence
            mutation_table[codon_to_remove]["pos"].remove(codon_idx)
            mutation_table[codon_to_remove]["difference"] -= 1

            # add to sequence
            mutation_table[codon_to_add]["pos"].append(codon_idx)
            mutation_table[codon_to_add]["difference"] += 1

            mutable_seq[codon_idx] = codon_to_add

            # update difference
            tot_diff = sum(
                abs(mutation_table[codon]["difference"]) for codon in synonymous_codons
            )
        # mutation_table now has difference = 0 for all codons

    return mutable_seq.toseq()


def mutate_codon(codon_in, codon_use_table):
    """Select a synonymous codon in accordance with the frequency of use
    in the host organism.

    Args:
        codon_in (Bio.Seq.Seq): A single codon.
        codon_use_table (dict({str : list[list, list]})): A dictionary with
            each amino acid three-letter code as keys, and a list of two
            lists as values. The first list is the synonymous codons that
            encode the amino acid, the second is the frequency with which
            each synonymous codon is used.

    Returns:
        Bio.Seq.Seq: A new codon.
    """
    AA = seq3(CodonTable.standard_dna_table.forward_table[str(codon_in)]).upper()

    synonymous_codons, codon_use_freq = codon_use_table[AA]
    if len(synonymous_codons) == 1:
        return codon_in

    # pick new codon
    codon_out = codon_in
    while codon_in == codon_out:
        codon_out = random.choices(synonymous_codons, codon_use_freq).pop()

    logger.detail(
        "mutating [{0}] codon from {1} to {2}".format(AA, codon_in[1], codon_out)
    )

    return codon_out


def resample_codons(dna_sequence, codon_use_table):
    """Generate a new DNA sequence by swapping synonymous codons.
    Codons are selected in accordance with their frequency of occurrence in
    the host organism.

    Args:
        dna_sequence (Bio.Seq.Seq): A read-only representation of
            the DNA sequence.
        codon_use_table (dict({str : list[list, list]})): A dictionary with
            each amino acid three-letter code as keys, and a list of two
            lists as values. The first list is the synonymous codons that
            encode the amino acid, the second is the frequency with which
            each synonymous codon is used.

    Returns:
        Bio.Seq.Seq: A read-only representation of the new DNA sequence.
    """
    resampled_dna = "".join(
        [
            random.choices(*codon_use_table[seq3(AA).upper()]).pop()
            for AA in dna_sequence.translate()
        ]
    )

    return Seq(resampled_dna, IUPAC.unambiguous_dna)


def gc_scan(dna_sequence, window_size, low, high):
    """Scan across a sequence and replace codons to acheive a desired GC
    content within the window.

    Args:
        dna_sequence (Bio.Seq.Seq): A read-only representation of
            the DNA sequence.
        window_size (int): Size of sliding window (in nucelotides) to
            examine for GC content. Window sizes can also be expressed
            as factors of the length of `dna_sequence` by passing a string
            that begins with "x" (e.g. "x0.5").
        low (float): Minimum GC content in window.
        high (float): Maximum GC content in window.

    Returns:
        Bio.Seq.Seq: A read-only representation of the new DNA sequence.
    """
    logger.info(
        "===== GC CONTENT SCAN IN WINDOW: {0} bps, threshold: {1} < x < {2}=====".format(
            window_size, low, high
        )
    )

    # some windows may be expressed as function of the sequence length
    if isinstance(window_size, str) and window_size.startswith("x"):
        window_size = int(float(window_size[1:]) * len(dna_sequence))

    # iterate across overlapping chunks of complete codons
    codon_window = window_size // 3
    overlap = codon_window // 2
    mutable_seq = dna_sequence.tomutable()

    # iterate by codon, but map back to sequence-based indices
    for i in range(0, len(mutable_seq) // 3, (codon_window - overlap) * 3):
        window = slice(
            i * 3,
            (i + codon_window) * 3
            if (i + codon_window) * 3 < len(mutable_seq)
            else len(mutable_seq),
        )
        logger.debug("Current segment: {0}".format(mutable_seq[window]))

        gc_percent = GC(mutable_seq[window]) / 100
        count = 0  # counter to prevent infinite loop
        # check gc_percent of current segment
        while (gc_percent < low or gc_percent > high) and count < codon_window * 2:
            position = random.randrange(0, len(mutable_seq[window]), 3)
            codon_idx = slice((i * 3) + position, ((i + 1) * 3) + position)

            init_codon = mutable_seq[codon_idx]
            new_codon = mutate_codon(init_codon, codon_use_table)

            if (GC(new_codon) < GC(init_codon) and gc_percent > high) or (
                GC(new_codon) > GC(init_codon) and gc_percent < low
            ):
                mutable_seq[codon_idx] = new_codon
                logger.debug("Mutating position: {0}".format(position))
                gc_percent = GC(mutable_seq[window]) / 100

            count += 1

    return mutable_seq.toseq()


def remove_restriction_sites(dna_sequence, restrict_sites):
    """Identify and remove seuences recognized by a set of restriction
    enzymes.

    Args:
        dna_sequence (Bio.Seq.Seq): A read-only representation of
            the DNA sequence.
        restrict_sites (Bio.Restriction.RestrictionBatch): RestrictionBatch
            instance configured with the input restriction enzymes.

    Returns:
        Bio.Seq.Seq: A read-only representation of the new DNA sequence.
    """

    logger.info("===== REMOVE RESTRICTION SITES =====")

    # check each unwanted restriction site
    analysis = Analysis(restrictionbatch=restrict_sites, sequence=dna_sequence)
    result = analysis.full()

    mutable_seq = dna_sequence.tomutable()
    for enz, cuts in result.items():
        for cut in cuts:
            logger.info(
                "The restriction enzyme {0} can cut the sequence before position {1}!".format(
                    str(enz), cuts
                )
            )
            # map sequence position to codon position
            # subtract 1 from `cut` to convert from sequence to string indexing
            codon_pos, offset = divmod((cut - 1) - (enz.size // 2), 3)

            # ensure the whole codon we mutate is being recognized by the restriction enzyme
            if offset:
                codon_pos += 1
            codon_idx = slice(codon_pos * 3, (codon_pos + 1) * 3)

            new_codon = mutate_codon(mutable_seq[codon_idx], codon_use_table)
            mutable_seq[codon_idx] = new_codon

    return mutable_seq.toseq()


def remove_start_sites(dna_sequence, ribosome_binding_sites, table_name="Standard"):
    """Identify and remove alternate start sites using a supplied set of
    ribosome binding sites and a codon table name.

    Args:
        dna_sequence (Bio.Seq.Seq): A read-only representation of
            the DNA sequence.
        ribosome_binding_sites (dict({str : str})): A dictionary with named
            ribosome binding sites as keys and the corresponding sequences
            as values.
        table_name (str, optional): Name of a registered NCBI table. See
            `Bio.Data.CodonTable.unambiguous_dna_by_name.keys()` for
            options. Defaults to "Standard".

    Returns:
        Bio.Seq.Seq: A read-only representation of the new DNA sequence.
    """
    codon_table = CodonTable.unambiguous_dna_by_name[table_name]
    logger.info(
        "===== REMOVE START SITES: {0} =====".format(
            ", ".join(codon_table.start_codons)
        )
    )

    # find all start codon sites (xTG)
    start_codon_positions = [
        m.start()
        for start_codon in codon_table.start_codons
        for m in re.finditer(start_codon, str(dna_sequence))
    ]

    if not len(start_codon_positions):
        logger.info("No start codons found in sequence")
        return dna_sequence

    logger.info(
        "Found {0} start codons. Checking for upstream RBSs...".format(
            len(start_codon_positions)
        )
    )

    # check each start site for RBS
    # 18 base pairs upstream of each xTG, ignore 3 bp closest to xTG
    _rbs_offset = 18
    rbs_positions = [
        pos - _rbs_offset for pos in start_codon_positions if pos >= _rbs_offset
    ]
    mutable_seq = dna_sequence.tomutable()

    for rbs_start in rbs_positions:
        # ignore 3 bp closest to xTG
        rbs_stop = rbs_start + _rbs_offset - 3
        rbs_query_seq = str(mutable_seq[rbs_start:rbs_stop])

        logger.detail(
            "checking sequence: {0}.{1}".format(
                rbs_query_seq, mutable_seq[rbs_stop : rbs_stop + 6]
            )
        )

        # check each unwanted RBS in each potential fragment
        for rbs, site in ribosome_binding_sites.items():
            logger.detail("checking start site: {0}, {1}".format(rbs, site))
            search = rbs_query_seq.find(site)

            count = 0  # counter to prevent infinite loop
            while search != -1 and count < 10:
                # mutate residues if site is found
                codon_pos = (search + rbs_start) // 3
                for i in range(2):
                    codon_idx = slice((codon_pos + i) * 3, (codon_pos + i + 1) * 3)
                    new_codon = mutate_codon(mutable_seq[codon_idx], codon_use_table)
                    mutable_seq[codon_idx] = new_codon

                # reset sequence and search again
                rbs_query_seq = str(mutable_seq[rbs_start : rbs_stop + 3])
                search = rbs_query_seq.find(site)
                count += 1

    return mutable_seq.toseq()


def remove_repeating_sequences(dna_sequence, window_size):
    """Idenify and remove repeating sequences of codons or groups of
    codons within a DNA sequence.

    Args:
        dna_sequence (Bio.Seq.Seq): A read-only representation of
            the DNA sequence.
        window_size (int): Size the window (in nucleotides) to examine.
            Window sizes are adjusted down to the nearest multiple of 3 so
            windows only contain complete codons.

    Returns:
        Bio.Seq.Seq: A read-only representation of the new DNA sequence.
    """
    logger.info(
        "===== REPEAT FRAGMENT SCAN FOR SIZE: {0} bps =====".format(window_size)
    )

    def _mutate_and_keep_looping(mutable_seq, window, offset):
        num_mutations = random.randint(1, 2)
        logger.debug("Mutating {0} codons".format(num_mutations))
        for _ in range(num_mutations):
            position = random.randrange(0, len(mutable_seq[window]), 3)
            codon_idx = slice(offset + position, (offset + 3) + position)
            new_codon = mutate_codon(mutable_seq[codon_idx], codon_use_table)
            mutable_seq[codon_idx] = new_codon

        return True

    # iterate across overlapping chunks of complete codons
    codon_window = window_size // 3
    overlap = codon_window - 1
    mutable_seq = dna_sequence.tomutable()

    current_cycle = 0  # prevent infinite loops (caused by poly-TRP or poly-MET)
    keep_looping = True
    # `keep_looping` if any mutation is made,
    # i.e. continue until both checks pass without mutations
    while keep_looping and (current_cycle < (codon_window * 10)):

        keep_looping = False

        # iterate by codon, but map back to sequence-based indices
        for i in range(0, len(mutable_seq) // 3, (codon_window - overlap) * 3):
            window = slice(
                i * 3,
                (i + codon_window) * 3
                if (i + codon_window) * 3 < len(mutable_seq)
                else len(mutable_seq),
            )

            # make each mutable codon immutable so it can be hashed later
            codons = [
                str(mutable_seq[window][i : i + 3])
                for i in range(0, len(mutable_seq[window]), 3)
            ]

            # check if all codons in the window are identical
            if len(set(codons)) == 1:
                logger.detail("All codons in window are identical: {0}".format(codons))
                keep_looping = _mutate_and_keep_looping(mutable_seq, window, (i * 3))

            # check if the segment is found in the full sequence
            non_overlapping_matches = re.findall(
                str(mutable_seq[window]), str(mutable_seq)
            )
            if len(non_overlapping_matches) > 1 and len(mutable_seq[window]) > 3:
                logger.debug("Current window is repeated in the sequence")
                keep_looping = _mutate_and_keep_looping(mutable_seq, window, (i * 3))

        current_cycle += 1

    return mutable_seq.toseq()


def remove_local_homopolymers(dna_sequence, n_codons=2):
    """Identify and remove consecutive streches of the same nucleotides
    using a sliding window of a fixed number of codons.

    Args:
        dna_sequence (Bio.Seq.Seq): A read-only representation of
            the DNA sequence.
        n_codons (int, optional): Size of window (in codons) to examine.
            Defaults to 2.

    Returns:
        Bio.Seq.Seq: A read-only representation of the new DNA sequence.
    """
    logger.info("===== REMOVE LOCAL HOMOPOLYMERS =====")
    mutable_seq = dna_sequence.tomutable()

    # look at each 6-mer
    keep_looping = True
    while keep_looping:
        for i in range(0, len(mutable_seq), 3):
            window = slice(
                i,
                i + (n_codons * 3)
                if i + (n_codons * 3) < len(mutable_seq)
                else len(mutable_seq),
            )

            seq = str(mutable_seq[window])
            nt_counts = {letter: seq.count(letter) for letter in set(seq)}
            letter = max(nt_counts, key=lambda letter: nt_counts[letter])

            if nt_counts[letter] <= args.local_homopolymer_threshold:
                keep_looping = False
                continue

            logger.detail("position: {0}: {1}".format(i, seq))
            logger.detail("{0}, count={1}".format(letter, nt_counts[letter]))

            for j in range(n_codons):
                codon_idx = slice(i + (j * 3), i + ((j + 1) * 3))
                mutable_seq[codon_idx] = mutate_codon(
                    mutable_seq[codon_idx], codon_use_table
                )
            keep_looping = True

    return mutable_seq.toseq()


##########################################################
#
# 	ACTUAL SCRIPT STARTS HERE
#
##########################################################

logger.info("===== SCRIPT START =====")

# read the input sequence file and parse lines
records = list(SeqIO.parse(args.input, "fasta", IUPAC.protein))

logger.info("Total number of sequences: {0}".format(len(records)))
[logger.detail(record) for record in records]

# generate host profile
codon_use_table, host_profile, codoon_relative_adativeness = codon_use.host_codon_usage(
    args.host, args.host_threshold
)

# process through all supplied sequences
for count, record in enumerate(records):
    logger.info("===== PROCESSING SEQUENCE {0} ===== {1}".format(count + 1, record.id))

    # check input seq style
    logger.info("INPUT IS AA SEQUENCE")
    dna = record.seq.back_translate()

    logger.detail("===== DUMPING SEQUENCE =====")
    logger.detail(str(dna))

    # intialize bookkeeping variables
    difference, current_cycle, relax, best_cai, best_dna_seq = 1.0, 1, 1.0, 0.0, Seq("")

    # keep running while there are cycles AND difference between current and host is less than the % relax allowed
    while ((current_cycle <= args.cycles) or (args.cycles == 0)) and (
        difference >= (relax - 1)
    ):
        logger.info(
            "~~~~~~~~~~ Current cycle: {0}/{1} ~~~~~~~~~~".format(
                current_cycle, args.cycles
            )
        )

        dna = resample_codons(dna, codon_use_table)

        # identify and remove undesirable features
        for _, gc_content in GC_content.items():
            # check various GC content requirements
            dna = gc_scan(dna, **gc_content)
        dna = remove_restriction_sites(dna, get_restriction_enzymes())
        dna = remove_start_sites(dna, RibosomeBindingSites, "Standard")
        dna = remove_repeating_sequences(dna, 9)
        dna = remove_local_homopolymers(dna)

        # count codon usage in the sequence
        count_table = codon_use.count_codons(dna)

        # adjust harmonization relax
        relax = 1 + (args.max_relax * ((current_cycle - 1) / args.cycles))
        logger.info("Relax coeff: {0}".format(relax))

        # measure the deviation from the host profile and adjust accordingly
        mutation_table, difference = compare_profiles(count_table, host_profile, relax)
        dna = harmonize_codon_use_with_host(dna, mutation_table)

        # if the codon adaptation index is better than what we've
        # seen so far, store this sequence
        cai = codoon_relative_adativeness.cai_for_gene(str(dna))
        if cai > best_cai:
            best_dna_seq = dna

        # tick cycle
        current_cycle += 1

    # hit the max number of cycles?
    if current_cycle > args.cycles:
        logger.info("You hit the max number of cycles: {0}".format(args.cycles))

    # check GC content
    logger.info("===== GC CONTENT =====")
    gc_percent = round(GC(best_dna_seq) / 100, 3)
    if gc_percent < 0.3 or gc_percent > 0.65:
        logger.warning("Overall GC content is {0}!".format(gc_percent))

    # display result name and sequence
    logger.output("===== SEQUENCE NAME =====")
    logger.output("{0}".format(record.id))

    logger.output(
        "Final codon-use difference between host and current sequence "
        + "(0.00 is ideal): {0}".format(round(difference, 2))
    )
    logger.output("===== DUMPING SEQUENCE =====")
    logger.output(str(best_dna_seq))
