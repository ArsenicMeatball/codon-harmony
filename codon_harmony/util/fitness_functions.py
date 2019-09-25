import re
from Bio.Data import CodonTable
from Bio.Restriction import Analysis
from Bio.SeqUtils import GC


from test_functions import codon_use_table
from . import codon_use, seq, MutableSeq
# consensus donor seq is "GGGTRAGT"
# below is all possible versions with the first GT fixed, and at
# least 2 other NTs from the consensus seq
_splice_donors = [
    re.compile(r"GGGT\wAGT", re.UNICODE),
    re.compile(r"\wGGT\w[AT]GT", re.UNICODE),
    re.compile(r"G\wGT\wAGT", re.UNICODE),
    re.compile(r"\w\wGT[AG]AGT", re.UNICODE),
    re.compile(r"GGGT[AG]\w[AG]T", re.UNICODE),
    re.compile(r"GGGT[AG]\wG\w", re.UNICODE),
    re.compile(r"GGGT[AG]A[AG]\w", re.UNICODE),
    re.compile(r"GGGT[AG]AG\w", re.UNICODE),
    re.compile(r"GGGT[AG]A[AG]T", re.UNICODE),
    re.compile(r"GGGT[AG]\wGT", re.UNICODE),
    # re.compile(r"\wGGT[AG]A[AG]\w", re.UNICODE), # redundant with below
    re.compile(r"\wGGT[AG][AT][ATG]\w", re.UNICODE),
    re.compile(r"\wGGT[AG]\wGT", re.UNICODE),
    re.compile(r"G\wGT[AG]A[AG]\w", re.UNICODE),
    re.compile(r"G\wGT[AG]A[AG]\w", re.UNICODE),
    re.compile(r"G\wGT[AG]AG\w", re.UNICODE),
    re.compile(r"G\wGT[AG]A[AG]T", re.UNICODE),
    re.compile(r"G\wGT[AG]\wGT", re.UNICODE),
    re.compile(r"\w\wGT[AG]AG\w", re.UNICODE),
    re.compile(r"\w\wGT[AG]\wGT", re.UNICODE),
    re.compile(r"\w\wGT[AG]\wG\w", re.UNICODE),
]

# consensus branch seq is "YTRAC"
# ignore branch points (for now) because they are small
# and occur 20-50 NTs upstream of acceptor -- not specific enough

# consensus acceptor seq is "YYYYYNCAGG"
# below are all sequences ending in NCAGG, NNAGG and NCAGN
# where at least 3 of the 5 upstream NTs are pyrimidines (Y, [TC])
_splice_acceptors = [
    re.compile(r"[TC][TC][TC]\w\w[ATCG]CAG\w", re.UNICODE),
    re.compile(r"[TC][TC]\w[TC]\w[ATCG]CAG\w", re.UNICODE),
    re.compile(r"[TC]\w[TC][TC]\w[ATCG]CAG\w", re.UNICODE),
    re.compile(r"\w[TC][TC][TC]\w[ATCG]CAG\w", re.UNICODE),
    re.compile(r"\w[TC][TC]\w[TC][ATCG]CAG\w", re.UNICODE),
    re.compile(r"\w[TC]\w[TC][TC][ATCG]CAG\w", re.UNICODE),
    re.compile(r"\w\w[TC][TC][TC][ATCG]CAG\w", re.UNICODE),
    re.compile(r"[TC]\w[TC]\w[TC][ATCG]CAG\w", re.UNICODE),
    re.compile(r"[TC]\w\w[TC][TC][ATCG]CAG\w", re.UNICODE),
    re.compile(r"[TC][TC]\w\w[TC][ATCG]CAG\w", re.UNICODE),
]


def eval_host(sequence, ancestor_sequence):
    return sum(['host'] for a, b in zip(sequence, ancestor_sequence) if a != b)


def eval_restriction_sites(sequence, restrict_sites):
    """Identify sequences with restriction enzymes
    and give a score based on the number of restriction sites found
    """
    # check unwanted restriction sites
    analysis = Analysis(restrictionbatch=restrict_sites, sequence=sequence)
    result = analysis.full()
    # score the sequence based on the number of restriction sites
    score = 0
    for enz, cuts in result.items():
        for cut in cuts:
            score += 1
    return score


def eval_start_sites(sequence, ribosome_binding_sites=rbs, table_name="Standard"):
    """Identify alternate start sites using a supplied set of
    ribosome binding sites and a codon table name.
    Then gives a score based on the number of alternate start sites found

    TODO: crossbreeding function

    Args:
        ribosome_binding_sites (dict{str, str}): A dictionary with named
            ribosome binding sites as keys and the corresponding sequences
            as values.
        table_name (str, optional): Name of a registered NCBI table. See
            `Bio.Data.CodonTable.unambiguous_dna_by_name.keys()` for
            options. Defaults to "Standard".

    Returns:
        The score for the fitness function
    """
    codon_table = CodonTable.unambiguous_dna_by_name[table_name]

    # find all start codon sites (xTG)
    start_codon_positions = [
        m.start()
        for start_codon in codon_table.start_codons
        for m in re.finditer(start_codon, str(sequence))
    ]

    # None found
    if not len(start_codon_positions):
        return 0

    # check each start site for RBS
    # 18 base pairs upstream of each xTG, ignore 3 bp closest to xTG
    _rbs_offset = 18
    rbs_positions = [
        pos - _rbs_offset for pos in start_codon_positions if pos >= _rbs_offset
    ]
    mutable_seq = sequence.tomutable()

    score = 0
    for rbs_start in rbs_positions:
        # ignore 3 bp closest to xTG
        rbs_stop = rbs_start + _rbs_offset - 3
        rbs_query_seq = str(mutable_seq[rbs_start:rbs_stop])

        # check each unwanted RBS in each potential fragment
        for rbs, site in ribosome_binding_sites.items():
            search = rbs_query_seq.find(site)

            count = 0  # counter to prevent infinite loop
            while search != -1 and count < 10:
                # mutate residues if site is found
                codon_pos = (search + rbs_start) // 3
                for ii in range(2):
                    codon_idx = slice((codon_pos + ii) * 3, (codon_pos + ii + 1) * 3)
                    score += 1

                # reset sequence and search again
                rbs_query_seq = str(mutable_seq[rbs_start: rbs_stop + 3])
                search = rbs_query_seq.find(site)
                count += 1
    return score


def eval_repeats(sequence, window_size=9):
    """Identify repeating sequences of nucleotides within a DNA sequence.
    assign a score based on the number of repeats

    Args:
        window_size (int): Size the window (in nucleotides) to examine.
            Window sizes are adjusted down to the nearest multiple of 3 so
            windows only contain complete codons.

    Returns:
        score: (int) representing the repeats
        :param sequence:
    """

    # iterate across overlapping chunks of complete codons
    codon_window = window_size // 3
    mutable_seq = sequence.tomutable()

    score = 0

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
                score += 1
            # check if the segment is found in the full sequence
            non_overlapping_matches = re.findall(
                str(mutable_seq[window]), str(mutable_seq)
            )
            if len(non_overlapping_matches) > 1 and len(mutable_seq[window]) > 3:
                score += len(non_overlapping_matches)
        current_cycle += 1
    return score


def eval_homopolymers(sequence, n_codons=2, homopolymer_threshold=4):
    """Identify consecutive stretches of the same nucleotides
    using a sliding window of a fixed number of codons.
    Then give a score based on the number of these stretches found.

    TODO: crossbreeding function

    Args:
        n_codons (int, optional): Size of window (in codons) to examine.
            Defaults to 2.
        homopolymer_threshold (int): number of consecutive nucleotide
            repeats allowed. Defaults to 4.

    Returns:
        The score for the fitness function
    """
    mutable_seq = sequence.tomutable()

    score = 0

    # look at each (n_codons * 3)-mer
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
            nt_count = 0
            nt_count_largest = 0
            last_letter = 'M'
            letter_largest = 'M'
            for letter in seq:
                if letter is last_letter:
                    nt_count += 1
                if nt_count > nt_count_largest:
                    nt_count_largest = nt_count
                    letter_largest = letter
                else:
                    # TODO: make loop stop if no larger sequence possible
                    last_letter = letter
                    nt_count = 1

                if nt_count_largest <= homopolymer_threshold:
                    keep_looping = False
                    continue

                score += nt_count // n_codons
                keep_looping = True
        return score


def eval_splice_sites(sequence):
    """Identify RNA splice sites within a DNA sequence.
    Then give a score based on number of splice sites in the sequence

    TODO: crossbreeding and change algorithm so it doesn't change

    Args:

    Returns:
        Bio.Seq.Seq: A read-only representation of the new DNA sequence.
    """

    def _pass_back_matches(list_of_sites, curr_dna):
        dna = str(curr_dna)
        sites = set(m for expr in list_of_sites for m in re.finditer(expr, dna))
        try:
            sites.remove(None)
        except KeyError:
            pass
        # remove redundancy
        sites = set((site.span(), site[0]) for site in sites)
        codon_bounds = [
            (s[0][0] // 3, -(-s[0][1] // 3)) for s in sorted(sites, key=lambda x: x[0])
        ]
        return codon_bounds

    def _get_splice_sites(curr_dna):
        donor_sites = _pass_back_matches(_splice_donors, curr_dna)
        acceptor_sites = _pass_back_matches(_splice_acceptors, curr_dna)
        return set(donor_sites + acceptor_sites)

    return len(_get_splice_sites(sequence.tomutable))
    """
    this code will be used to identify problem childs later
    mutable_seq = sequence.tomutable()
    score = 0
    keep_looping = True
    n_times_through_unchanged = 0
    prev_seq = ""
    while keep_looping:
        # look for donor and acceptor seqs
        splice_sites = _get_splice_sites(mutable_seq)

        logger.info("Removing RNA splice site donors and acceptors.")
        for site_bounds in splice_sites:
            site_removed = False

            codon_list = list(range(site_bounds[0], site_bounds[-1] + 1))
            random.shuffle(codon_list)

            indices = [slice(cdn * 3, (cdn + 1) * 3) for cdn in codon_list]
            init_codons = [mutable_seq[idx] for idx in indices if str(mutable_seq[idx])]

            amino_acids = [
                seq3(CodonTable.standard_dna_table.forward_table[str(init)]).upper()
                for init in init_codons
            ]
            synonymous_codons = [codon_use_table[amino_acid][0] for amino_acid in amino_acids]
            substitutions = list(product(*synonymous_codons))
            random.shuffle(substitutions)

            for substitution in substitutions:
                for idx, codon in zip(indices, substitution):
                    mutable_seq[idx] = codon
                    score += 5
                curr_sites = _get_splice_sites(mutable_seq)
                if site_bounds in curr_sites or len(curr_sites) >= len(splice_sites):
                    # put the starting sequence back
                    score -= 5
                    for idx, init_codon in zip(indices, init_codons):
                        mutable_seq[idx] = init_codon
                else:
                    site_removed = True
                    logger.detail("Removed site ({})!".format(site_bounds))
                    logger.detail(
                        "Remaining sites:\n" + "\n".join([str(s) for s in curr_sites])
                    )
                    break

            if site_removed:
                break

        remaining_sites = _get_splice_sites(mutable_seq)
        n_times_through_unchanged += int(str(mutable_seq) == prev_seq)
        prev_seq = str(mutable_seq)

        if not len(remaining_sites) or n_times_through_unchanged == 5:
            keep_looping = False
        """


def eval_gc_content(sequence, gc):
    """Scan across a sequence and determine a score based on how far it is
    from desired GC content

    Todo: crossbreeding

    Args:
        gc (GCParams): A `namedtuple` with fields for name, window_size,
            minimum and maximum GC content.

    Note:
        The following fields of the `GCParams` type are used in this
        function:

        * **window_size** (`int`) – Size of sliding window (in nucelotides) to
          examine for GC content. Window sizes can also be expressed as
          factors of the length of `dna_sequence` by passing a string
          that begins with "x" (e.g. "x0.5").

        * **low** (`float`) – Minimum GC content in window.

        * **high** (`float`) – Maximum GC content in window.
    """

    window_size = gc.window_size  # tuples are immutable
    # some windows may be expressed as function of the sequence length
    if isinstance(window_size, str) and window_size.startswith("x"):
        window_size = int(float(window_size[1:]) * len(sequence))

    # iterate across overlapping chunks of complete codons
    codon_window = window_size // 3 + 1
    mutable_seq = sequence.tomutable()
    score = 0
    # iterate by codon, but map back to sequence-based indices
    for i in range(len(mutable_seq) // 3):
        window = slice(
            i * 3,
            (i + codon_window) * 3
            if (i + codon_window) * 3 < len(mutable_seq)
            else len(mutable_seq),
            )

        gc_percent = GC(mutable_seq[window]) / 100

        if gc_percent > gc.high:
            score += gc_percent - gc.high
        elif gc_percent < gc.low:
            score += gc.low - gc_percent
    return score


def eval_hairpins(sequence, stem_length=10):
    """Identify possible hairpins in population.
    Score based on number of possible hairpins
        Todo: crossbreeding
        Args:
            stem_length (int, optional): Length of hairpin stem to detect.
                Defaults to 10.

        Returns:
            Score
        """
    mutable_seq = sequence.tomutable()
    score = 0
    for i in range(0, len(mutable_seq), stem_length):
        stem_seq = mutable_seq[i: i + stem_length].toseq()
        # include wobble base pairing for G-[CT]
        hairpin_pattern = "".join(
            [nt if nt != "C" else "[CT]" for nt in stem_seq.reverse_complement()]
        )
        for hairpin in re.finditer(hairpin_pattern, str(mutable_seq)):
            score += 1
    return score