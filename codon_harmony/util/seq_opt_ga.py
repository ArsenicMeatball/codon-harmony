import logging
import random
import re

from itertools import product

from Bio.Alphabet import IUPAC

from Bio.Data import CodonTable
from Bio.Restriction import Analysis
from Bio.SeqUtils import CodonUsage, GC, seq3

from . import Seq, MutableSeq, codon_use

logger = logging.getLogger(__name__)


def the_sequence_optimizer(dna_sequence, codon_use_table, window_size):
    """Genetic algorithm that creates a population of sequences,
    does checks on them, mutates and repeats until powerrrrr.

    Args:
        dna_sequence (Bio.Seq.Seq): A read-only representation of
            the DNA sequence. initial population
        codon_use_table (dict{str, list[list, list]}): A dictionary with
            each amino acid three-letter code as keys, and a list of two
            lists as values. The first list is the synonymous codons that
            encode the amino acid, the second is the frequency with which
            each synonymouscodon is used.
        window_size (int): Size the window (in nucleotides) to examine.
            Window sizes are adjusted down to the nearest multiple of 3 so
            windows only contain complete codons.

    Returns:
        Bio.Seq.Seq: A read-only representation of the new DNA sequence.
    """
    logger.info("Scanning for repeated stretches of {} nucleotides".format(window_size))

    def _eval_host():
        return 0.5

    def _eval_restriction_sites(restrict_sites):
        """Identify sequences with restriction enzymes
        and give a score based on the number of sequences afflicted

        TODO: crossbreeding function

        Args:
            restrict_sites (Bio.Restriction.RestrictionBatch): RestrictionBatch
                instance configured with the input restriction enzymes.

        Returns:
            score: the score for the fitness function
        """

        logger.info("Testing restriction sites")

        # check unwanted restriction sites
        analysis = Analysis(restrictionbatch=restrict_sites, sequence=dna_sequence)
        result = analysis.full()

        # score the sequence based on the number of restriction sites
        score = 0
        for enz, cuts in result.items():
            for cut in cuts:
                logger.info(
                    "Restriction enzyme ({}) cut site detected at position {}.".format(
                        str(enz), cuts
                    )
                )
                score += 1
        return score

    def _eval_start_sites(ribosome_binding_sites, table_name="Standard"):
        """Identify and remove alternate start sites using a supplied set of
        ribosome binding sites and a codon table name.

        Args:
            ribosome_binding_sites (dict{str, str}): A dictionary with named
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
            "Testing alternate start sites: {}".format(", ".join(codon_table.start_codons))
        )

        # find all start codon sites (xTG)
        start_codon_positions = [
            m.start()
            for start_codon in codon_table.start_codons
            for m in re.finditer(start_codon, str(dna_sequence))
        ]

        # None found
        if not len(start_codon_positions):
            logger.info("No start codons found in sequence")
            return 0


        logger.info(
            "Found {} start codons. Checking for upstream RBSs...".format(
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
                'checking for sequence "{}" in "{}"'.format(
                    rbs_query_seq, mutable_seq[rbs_stop: rbs_stop + 6]
                )
            )

            # check each unwanted RBS in each potential fragment
            for rbs, site in ribosome_binding_sites.items():
                logger.detail('checking for start site "{}" in "{}"'.format(rbs, site))
                search = rbs_query_seq.find(site)

                count = 0  # counter to prevent infinite loop
                while search != -1 and count < 10:
                    # mutate residues if site is found
                    codon_pos = (search + rbs_start) // 3
                    for ii in range(2):
                        codon_idx = slice((codon_pos + ii) * 3, (codon_pos + ii + 1) * 3)
                        new_codon = mutate_codon(mutable_seq[codon_idx], codon_use_table)
                        mutable_seq[codon_idx] = new_codon

                    # reset sequence and search again
                    rbs_query_seq = str(mutable_seq[rbs_start: rbs_stop + 3])
                    search = rbs_query_seq.find(site)
                    count += 1

        return mutable_seq.toseq()

    def _eval_repeats(mutable_seq, window, offset):
        return 0.5

    def _eval_homopolymers():
        return 0.5

    def _eval_splice_sites():
        return 0.5

    def _eval_gc_content():
        return 0.5

    def _eval_sequence():
        return 0.5

    def _mutate_and_keep_looping(mutable_seq, window, offset):
        num_mutations = random.randint(1, 2)
        logger.debug("Mutating {} codons".format(num_mutations))
        for _ in range(num_mutations):
            position = random.randrange(0, len(mutable_seq[window]), 3)
            codon_idx = slice(offset + position, (offset + 3) + position)
            new_codon = mutate_codon(mutable_seq[codon_idx], codon_use_table)
            mutable_seq[codon_idx] = new_codon

        return True

    # iterate across overlapping chunks of complete codons
    codon_window = window_size // 3
    mutable_seq = dna_sequence.tomutable()

    current_cycle = 0  # prevent infinite loops (caused by poly-TRP or poly-MET)
    keep_looping = True
    # `keep_looping` if any mutation is made,
    # i.e. continue until both checks pass without mutations
    while keep_looping and (current_cycle < (codon_window * 10)):

        keep_looping = False

        # iterate by codon, but map back to sequence-based indices
        for i in range(len(mutable_seq) // 3):
            window = slice(
                i * 3,
                (i + codon_window) * 3
                if (i + codon_window) * 3 < len(mutable_seq)
                else len(mutable_seq),
            )

            # make each mutable codon immutable so it can be hashed later
            codons = [
                str(mutable_seq[window][i: i + 3])
                for i in range(0, len(mutable_seq[window]), 3)
            ]

            # check if all codons in the window are identical
            if len(set(codons)) == 1:
                logger.detail("All codons in window are identical: {}".format(codons))
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