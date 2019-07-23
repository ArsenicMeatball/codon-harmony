#! /usr/bin/env python
import argparse
import itertools
import multiprocessing
import random
import sys

from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import GC

from .util import codon_use, logging, log_levels, seq_opt
from .data import GC_content, RibosomeBindingSites, RestrictionEnzymes

logger = logging.getLogger(__name__)


def get_parser():
    parser = argparse.ArgumentParser(
        description="Reverse translate your amino acid sequence harmoniously with "
        + "a host's condon usage.",
        epilog="v1.0.0 (contact bweitzner@lyellbio.com if you encounter errors)",
    )
    parser.add_argument(
        "--input", type=str, required=True, help="input file with sequence"
    )
    parser.add_argument(
        "--output", type=str, required=True, help="output file to write DNA sequence(s)"
    )
    parser.add_argument(
        "--host",
        type=str,
        default="413997",
        help="host table code: http://www.kazusa.or.jp/codon/, default is "
        + '"Escherichia coli B"',
    )
    parser.add_argument(
        "--host-threshold",
        type=float,
        default="0.10",
        help="lowest codon fraction per AA in the host that is allowed",
    )
    parser.add_argument(
        "--local-host-profile",
        type=str,
        default=None,
        help="path to host codon usage table as JSON file",
    )
    parser.add_argument(
        "--verbose",
        type=int,
        default=0,
        choices=[0, 1, 2, 3],
        help="verbose output level (0=only result, 1=standard output, "
        + "2=extra output 3=debugging)",
    )
    parser.add_argument(
        "--local-homopolymer-threshold",
        type=int,
        default="4",
        help="number of consecutive NT repeats allowed",
    )
    parser.add_argument(
        "--cycles",
        type=int,
        default=10,
        help="number of independent codon samples to run. 0 means 1 pass",
    )
    parser.add_argument(
        "--inner-cycles",
        type=int,
        default=10,
        help="number of times to iteratively optimize each independent codon"
        + " sample. 0 means 1 pass",
    )
    parser.add_argument(
        "--max-relax",
        type=float,
        default="0.1",
        help="maximum percent deviation from host profile",
    )
    parser.add_argument(
        "--restriction-enzymes",
        nargs="*",
        type=str,
        default="NdeI XhoI HpaI PstI EcoRV NcoI BamHI".split(),
        help="list of restriction enzyme sites to remove "
        + "(e.g. --restriction_enzymes NdeI XhoI HpaI). ",
    )

    splice_sites = parser.add_mutually_exclusive_group(required=False)
    splice_sites.add_argument(
        "--remove-splice-sites",
        dest="splice_sites",
        action="store_true",
        help="Remove splice sites. Use for mammalian hosts.",
    )
    splice_sites.add_argument(
        "--no-remove-splice-sites",
        dest="splice_sites",
        action="store_false",
        help="Do not remove splice sites.",
    )
    parser.set_defaults(splice_sites=True)

    start_sites = parser.add_mutually_exclusive_group(required=False)
    start_sites.add_argument(
        "--remove-start-sites",
        dest="start_sites",
        action="store_true",
        help="Remove alternate start sites. Use for bacterial hosts.",
    )
    start_sites.add_argument(
        "--no-remove-start-sites",
        dest="start_sites",
        action="store_false",
        help="Do not remove alternate start sites.",
    )
    parser.set_defaults(start_sites=True)
    parser.add_argument(
        "--one-line-fasta",
        dest="one_line_fasta",
        action="store_true",
        help="Remove splice sites. Use for mammalian hosts.",
    )

    return parser


def _harmonize_sequence(
    seq_record,
    args,
    codon_use_table,
    host_profile,
    codon_relative_adativeness,
    rest_enz,
    seq_no=0,
    one_line_fasta=False,
):
    """Convert an amino acid sequence to DNA and optimize it to be synthesizable
    and to match a host codon usage profile.

    Args:
        seq_record (Bio.SeqRecord.SeqRecord): The amino acid sequence to be
            reverse translated and optimized.
        args ([type]): A namespace populated with attributes from the
            command line argument parser.
        codon_use_table (dict{str, list[list, list]}): A dictionary with
            each amino acid three-letter code as keys, and a list of two
            lists as values. The first list is the synonymous codons that
            encode the amino acid, the second is the frequency with which
            each synonymous codon is used.
        host_profile (dict{str, int}): A dictionary with codons as keys
            and the corresponding frequency of occurences as values.
        codon_relative_adativeness (Bio.SeqUtils.CodonUsage.CodonAdaptationIndex):
            A `CodonAdaptationIndex` instance configured to calculate CAI
            for a target gene.
        rest_enz (Bio.Restriction.Restriction.RestrictionBatch):
            RestrictionBatch instance configured with the input restriction
            enzymes.
        seq_no (int, optional): The number of the sequence in the input file.
            Used for logging purposes. Defaults to 0.
        one_line_fasta (bool, optional): Toggles the format of the output DNA
            sequence. If True, the sequence will be on one line without any
            line wrapping. Defaults to False.

    Returns:
        dict{str, str}: A dictionary with "protein" or "dna" as a key, and the
        input protein- or optimized DNA-sequence in FASTA format as a value. The
        protein sequence is returned if the function cannot satisify the
        harmonization requirements.
    """

    out_format = "fasta-2line" if one_line_fasta else "fasta"
    logger.info(
        "Processing sequence number {}:\n{}".format(
            seq_no + 1, seq_record.format(out_format)
        )
    )

    dna = seq_record.seq.back_translate()
    logger.detail(
        "Initial DNA sequence:\n{}".format(
            SeqRecord(dna, id=seq_record.id).format(out_format)
        )
    )

    # intialize bookkeeping variables
    best_cai, best_dna = 0.0, ""

    args.cycles = 1 if args.cycles == 0 else args.cycles
    args.inner_cycles = 1 if args.inner_cycles == 0 else args.inner_cycles

    # run `args.cycles` independent trials
    for sample_no in range(args.cycles):
        logger.info("Current sample no: {}/{}".format(sample_no + 1, args.cycles))

        # relax harmonization requirement
        relax = 1 + (args.max_relax * ((sample_no + 1) / (args.cycles)))
        dna = seq_opt.resample_codons_and_enforce_host_profile(
            dna, codon_use_table, host_profile, relax
        )

        # go through a few cycles with the same starting sequence to
        # allow iterative improvements to the same sample of codons
        for _ in range(args.inner_cycles):
            # identify and remove undesirable features
            for gc_content in GC_content:
                # check various GC content requirements
                dna = seq_opt.gc_scan(dna, codon_use_table, gc_content)

            if args.start_sites:
                dna = seq_opt.remove_start_sites(
                    dna, codon_use_table, RibosomeBindingSites
                )

            dna = seq_opt.remove_repeating_sequences(dna, codon_use_table, 9)
            dna = seq_opt.remove_local_homopolymers(
                dna,
                codon_use_table,
                n_codons=2,
                homopolymer_threshold=args.local_homopolymer_threshold,
            )
            dna = seq_opt.remove_hairpins(dna, codon_use_table, stem_length=10)
            if len(rest_enz):
                dna = seq_opt.remove_restriction_sites(dna, codon_use_table, rest_enz)

        # measure the deviation from the host profile post-cleanup
        # only move forward if we haven't deviated too much from host
        _, difference = seq_opt.compare_profiles(
            codon_use.count_codons(dna), host_profile, relax
        )
        if difference >= args.max_relax:
            continue

        # if the codon adaptation index is better than what we've
        # seen so far, store this sequence
        cai = codon_relative_adativeness.cai_for_gene(str(dna))
        if cai > best_cai:
            best_cai = cai
            best_dna = SeqRecord(
                dna,
                id=seq_record.id,
                name=seq_record.name,
                description=seq_record.description,
            )

    logger.info(
        "Completed {} independent codon samples with optimization!".format(args.cycles)
    )

    if isinstance(best_dna, str):
        logger.warning(
            "Unable to create suitable DNA sequence for input sequence {}.\n{}".format(
                seq_no + 1, seq_record.format(out_format)
            )
        )
        return {"protein": seq_record.format(out_format)}

    if args.splice_sites:
        logger.output("Detecting and removing splice sites before outputting.")
        best_dna.seq = seq_opt.remove_splice_sites(best_dna.seq, codon_use_table)

    logger.output("Optimized gene metrics and sequence")
    # check GC content
    gc_frac = GC(best_dna.seq) / 100
    logger.output("Final overall GC content is {:.0%}".format(gc_frac))
    if gc_frac < 0.3 or gc_frac > 0.65:
        logger.warning(
            "The sequence's GC content ({:.2f}) is beyond normal ranges (0.3 > GC < 0.65)!".format(
                gc_frac
            )
        )

    # measure the final deviation from the host profile
    _, difference = seq_opt.compare_profiles(
        codon_use.count_codons(best_dna.seq), host_profile, relax
    )

    logger.output(
        "Final codon-use difference between host and current sequence: {:.2f}".format(
            difference
        )
    )

    best_dna_fasta = best_dna.format(out_format)
    logger.output(
        "The designed gene's CAI is: {:.2f}\n{}".format(best_cai, best_dna_fasta)
    )
    return {"dna": best_dna_fasta}



class AAForm:
    '''
    Class to take in arguments
    (dictionary with at least inputfile specified) 
    not from command-line,
    used by runner.py(the_dict) in scripts
    '''
    def __init__(self, argv):
        self.host = "413997"
        self.host_threshold = 0.1
        self.local_homopolymer_threshold = 4
        self.cycles = 10
        self.inner_cycles = 10
        self.max_relax = 0.1
        resties = "NdeI XhoI HpaI PstI EcoRV NcoI BamHI".split()
        self.restriction_enzymes = resties
        self.splice_sites = False
        self.start_sites = False
        self.local_host_profile = None
        self.verbose = 1
        self.one_line_fasta = False
        self.output = "output.fasta"
        self.run = True
        if argv.get("input"):
            self.input = argv.get('input')
        if argv.get("host") is not None:
            self.host = argv.get('host')
        if argv.get("host_threshold") is not None:
            self.host_threshold = float(argv.get('host_threshold'))
        if argv.get("local_homopolymer_threshold") is not None:
            resties = int(argv.get('local_homopolymer_threshold'))
            self.local_homopolymer_threshold = resties
        if argv.get("cycles") is not None:
            self.cycles = int(argv.get('cycles'))
        if argv.get("inner_cycles") is not None:
            self.inner_cycles = int(argv.get('inner_cycles'))
        if argv.get("max_relax") is not None:
            self.max_relax = float(argv.get('max_relax'))
        if argv.get("restriction_sites") is not None:
            self.restriction_enzymes = argv.get('restriction_sites').split()
        if argv.get("remove_splice_sites") is not None:
            self.splice_sites = argv.get('remove_splice_sites')
        if argv.get("remove_alternate_start_site") is not None:
            self.start_sites = argv.get('remove_alternate_start_site')
        if argv.get("host_profile") is not None:
            self.local_host_profile = argv.get('host_profile')
        # gc_richness_threshold = argv['gc_richness_threshold']
        # gc_richness_chunk_size = argv['gc_richness_chunk_size']
        # title = argv['title']

    def get_host(self):
        return self.host


def main(argv=None):
    """Read in a fasta-formatted file containing amino acid sequences and
    reverse translate each of them in accordance with a specified host's
    codon usage frequency. The DNA sequence is then processed to remove
    unwanted features.
    """
    if argv and isinstance(argv, AAForm) and argv.get("input"):
        args = AAForm(argv)
    else:
        args = get_parser().parse_args(argv)

    logging.basicConfig(level=log_levels[args.verbose])

    random.seed()
    logger.info("Beginning codon use optimization")

    # stores the final sequences
    out_seqs = []

    # generate host profile
    codon_use_table, host_profile, codon_relative_adativeness = codon_use.host_codon_usage(
        args.host, args.host_threshold, args.local_host_profile
    )

    # initialize the restriction sites of interest
    rest_enz = RestrictionEnzymes(args.restriction_enzymes)

    def _input_prep():
        """Read in all supplied sequences and prepare args for the
        ``_harmonize_sequence`` function
        """
        for seq_no, record in enumerate(
            SeqIO.parse(args.input, "fasta", IUPAC.protein)
        ):
            yield (
                record,
                args,
                codon_use_table,
                host_profile,
                codon_relative_adativeness,
                rest_enz,
                seq_no,
                args.one_line_fasta,
            )

    with multiprocessing.Pool() as pool:
        input_iter = _input_prep()
        # split input into chunks no larger than 3 x Ncpu
        N = multiprocessing.cpu_count() * 3

        # in principle, the number of sequences are in the fasta file is unknown, so
        # go through it with an islice and break out when there are no more sequences
        while True:
            seqs = pool.starmap(_harmonize_sequence, itertools.islice(input_iter, N))
            if seqs:
                out_seqs.extend(seqs)
            else:
                break

    dna_seqs = "".join([seq["dna"] for seq in out_seqs if "dna" in seq])
    failed_seqs = "".join([seq["protein"] for seq in out_seqs if "protein" in seq])

    # write sequences to file
    with open(args.output, "w") as f:
        f.write(dna_seqs)

    if failed_seqs:
        logger.warning(
            "Unable to create suitable DNA sequence for one or more input sequences.\n{}".format(
                failed_seqs
            )
        )
        logger.warning(
            'Writing failed sequences to "failed_seqs.fasta". '
            "Try re-running with a higher value for --max-relax"
        )
        with open("failed_seqs.fasta", "w") as f:
            f.write(failed_seqs)
