import math
import sys
import random
from Bio import Seq
from Bio.Data import CodonTable
from Bio.SeqUtils import seq3
import re
from Bio.Restriction import Analysis
from Bio.SeqUtils import GC

max_generations = 2


class SequenceContainer:

    def __init__(self, sequence):
        setattr(self, "start_fitness", sys.maxsize)
        setattr(self, "hairpins_fitness", sys.maxsize)
        assert sequence is Seq
        setattr(self, "sequence", sequence)
        setattr(self, "gc_fitness", sys.maxsize)
        setattr(self, "homo_fitness", sys.maxsize)
        setattr(self, "host_fitness", sys.maxsize)
        setattr(self, "repeats_fitness", sys.maxsize)
        setattr(self, "restriction_fitness", sys.maxsize)
        setattr(self, "splice_fitness", sys.maxsize)
        setattr(self, "dominated", False)
        setattr(self, "density", sys.maxsize)


def n_dimensional_euclidian_distance(sequence1, sequence2):
    assert sequence2 is SequenceContainer and sequence1 is SequenceContainer
    gc = (getattr(sequence2, "gc_fitness") - getattr(sequence1, "gc_fitness")) ** 2
    homo = (getattr(sequence2, "homo_fitness") - getattr(sequence1, "homo_fitness")) ** 2
    host = (getattr(sequence2, "host_fitness") - getattr(sequence1, "host_fitness")) ** 2
    repeat = (getattr(sequence2, "repeats_fitness") - getattr(sequence1, "repeats_fitness")) ** 2
    restriction = (getattr(sequence2, "restriction_fitness") - getattr(sequence1, "restriction_fitness")) ** 2
    splice = (getattr(sequence2, "splice_fitness") - getattr(sequence1, "splice_fitness")) ** 2
    start = (getattr(sequence2, "start_fitness") - getattr(sequence1, "start_fitness")) ** 2
    hairpin = (getattr(sequence2, "hairpins_fitness") - getattr(sequence1, "hairpins_fitness")) ** 2
    return math.sqrt(gc + homo + host + repeat + restriction + splice + start + hairpin)


def crossover(parents, probability_crossover):
    """number_of_crossovers = len(parents) * probability_crossover
    if number_of_crossovers is 0 or len(parents) < 2 or parents[0] is not SequenceContainer:
        return
    for current_num_crossovers in range(number_of_crossovers):
        parent2 = parent1 = random.randint(1, len(parents))
        while parent1 is parent2:
            parent2 = random.randint(1, len(parents))
        crossover_point = random.randint(1, len(parents-1))
        parents.append(getattr(parents[parent1], "sequence").tomutable()[:crossover_point].append(getattr(parents[parent2], "sequence").tomutable()[crossover_point-1:]))
    """
    return parents


def crossover_and_mutate(parents, probability_crossover, probability_mutation, codon_use_table):
    population = []
    for sequence in parents:
        population.append(mutate_sequence(sequence, codon_use_table, probability_mutation))
    population = crossover(population, probability_crossover)
    return population


def get_parents(archive):
    selected = archive
    selected.sort(
        key=lambda sequence:
        getattr(sequence, "gc_fitness") +
        getattr(sequence, "homo_fitness") +
        getattr(sequence, "host_fitness") +
        getattr(sequence, "repeats_fitness") +
        getattr(sequence, "restriction_fitness") +
        getattr(sequence, "splice_fitness") +
        getattr(sequence, "start_fitness") +
        getattr(sequence, "hairpins_fitness")
    )
    return selected[:0.2 * len(archive)]


def get_density(elem):
    return getattr(elem, "density")


def truncate_similar_individuals(archive, archive_size):
    sorted_archive = sorted(archive, key=get_density, reverse=False)
    return sorted_archive[archive_size:]


def populate_archive_with_remaining_best(archive, archive_size, dominated):
    num_individuals_to_add = archive_size - len(archive)
    dominated_set = sorted(dominated, key=lambda sequence:
                           getattr(sequence, "gc_fitness") +
                           getattr(sequence, "homo_fitness") +
                           getattr(sequence, "host_fitness") +
                           getattr(sequence, "repeats_fitness") +
                           getattr(sequence, "restriction_fitness") +
                           getattr(sequence, "splice_fitness") +
                           getattr(sequence, "start_fitness") +
                           getattr(sequence, "hairpins_fitness"), reverse=False)
    archive.append(dominated_set[num_individuals_to_add:])
    return archive


def get_non_dominated_solutions(union):
    """
    Determine pareto dominance, that is
    dominant: solution that has the best scores without betraying the other scores
    non-dominant: solution that betrays the scores, compared to another one
    its a two way relationship
    :param union:
    :return:
    """
    non_dominated_set = []
    dominated_set = []
    pareto_dominance_graph = [[True for x in range(len(union))] for y in range(len(union))]
    # if any vector is larger it is dominated, else it is not dominated
    for idx in union:
        for jdx in union:
            if getattr(union[idx], "gc_fitness") > getattr(union[jdx], "gc_fitness") or \
                    getattr(union[idx], "homo_fitness") > getattr(union[jdx], "homo_fitness") or \
                    getattr(union[idx], "host_fitness") > getattr(union[jdx], "host_fitness") or \
                    getattr(union[idx], "repeats_fitness") > getattr(union[jdx], "repeats_fitness") or \
                    getattr(union[idx], "restriction_fitness") > getattr(union[jdx], "restriction_fitness") or \
                    getattr(union[idx], "splice_fitness") > getattr(union[jdx], "splice_fitness") or \
                    getattr(union[idx], "start_fitness") > getattr(union[jdx], "start_fitness") or \
                    getattr(union[idx], "hairpins_fitness") > getattr(union[jdx], "hairpins_fitness"):
                pareto_dominance_graph[idx][jdx] = False
    for idx in union:
        add = True
        for jdx in union:
            if pareto_dominance_graph[idx][jdx] is False:
                add = False
                break
        if add is True:
            non_dominated_set.append(union[idx])
        else:
            dominated_set.append(union[idx])
    return non_dominated_set, dominated_set


def calculate_solution_density(union, nearest_neighbours):
    """
    Evaluate distance between kth nearest neighbours, where k is sqrt(len(union))
    :param nearest_neighbours:
    :param union:
    :return:
    """
    kth = math.sqrt(len(union))
    for idx in range(len(union)):
        assert union[idx] is SequenceContainer
        setattr(union[idx], "density", sorted(nearest_neighbours[idx])[kth])


def calculate_raw_fitness(individual):
    """
    Determine the fitness relative to the others based on how dominant they are
    :param individual:
    :return:
    """
    pass


def calculate_fitness(population, gc_parameters, ancestor_sequence, restriction_sites):
    for individual in population:
        assert(individual is SequenceContainer)
        sequence = getattr(individual, "sequence")
        population[sequence][0] = eval_gc_content(sequence, gc_parameters)
        population[sequence][1] = eval_homopolymers(sequence)
        population[sequence][2] = eval_host(sequence, ancestor_sequence)
        population[sequence][3] = eval_repeats(sequence)
        population[sequence][4] = eval_restriction_sites(sequence, restriction_sites)
        population[sequence][5] = eval_splice_sites(sequence)
        population[sequence][6] = eval_start_sites(sequence)
        population[sequence][7] = eval_hairpins(sequence)


def initialize_population(population_size, problem_size, ancestor_sequence, probability_mutation, codon_use_table):
    """
    Initializes a population of sequences based on a (hopefully) codon optimized
    sequence
    :param codon_use_table:
    :param probability_mutation: float determines probability of mutation
    :param population_size: int determines minimum population size
    :param problem_size: int helps determines minimum population size
    :param ancestor_sequence: Bio.seq.seq codon optimized sequence
    :return:
    """
    population = []
    mutable_seq = ancestor_sequence.tomutable()
    if problem_size < population_size:
        size = population_size
    else:
        size = problem_size
    for idx in range(size):
        population.append(
            SequenceContainer(mutate_sequence(mutable_seq.to_seq(), codon_use_table, probability_mutation))
        )
    return population


def mutate_sequence(individual, codon_use_table, mutation_probability=0.05, offset=0):
    """
    Takes a single sequence and gives it a random number of mutations.
    :param individual:
    :param codon_use_table:
    :param mutation_probability:
    :param offset:
    :return:
    """
    assert(individual is SequenceContainer)
    sequence = getattr(individual, "sequence")
    mutable_seq = sequence.tomutable()
    num_codons = len(mutable_seq) // 3
    num_mutations = math.ceil(num_codons * mutation_probability)
    for _ in range(num_mutations):
        position = 3 * random.randrange(0, len(mutable_seq) // 3)
        codon_idx = slice(offset + position, (offset + 3) + position)
        new_codon = mutate_codon(mutable_seq[codon_idx], codon_use_table)
        mutable_seq[codon_idx] = new_codon
    return mutable_seq.toseq()


def mutate_codon(codon_in, codon_use_table):
    """Select a synonymous codon in accordance with the frequency of use
    in the host organism.

    Args:
    codon_in (Bio.Seq.Seq): A single codon.

    Returns:
        Bio.Seq.Seq: A new codon.
    """
    amino_acid = seq3(CodonTable.standard_dna_table.forward_table[str(codon_in)]).upper()
    synonymous_codons, codon_use_freq = codon_use_table[amino_acid]
    if len(synonymous_codons) == 1:
        return codon_in

    # pick new codon
    codon_out = codon_in
    while codon_in == codon_out:
        codon_out = random.choices(synonymous_codons, codon_use_freq).pop()

    return codon_out


def determine_neighbours(union):
    neighbours = [[0 for x in range(len(union))] for y in range(len(union))]
    for idx in range(len(union)):
        for jdx in range(len(union)):

            if idx != jdx:
                distance = n_dimensional_euclidian_distance(union[idx], union[jdx])
                neighbours[idx][jdx] = neighbours[jdx][idx] = distance
    return neighbours


def optimize_with_strength_pareto_evolutionary_algorithm(population_size, archive_size, problem_size,
                                                         probability_crossover, probability_mutation,
                                                         ancestor_sequence, codon_use_table):
    population = initialize_population(population_size, problem_size, ancestor_sequence, probability_mutation, codon_use_table)
    archive = []
    union = []
    for generation in range(0, max_generations):
        calculate_fitness(population, problem_size)
        union = population.append(archive)
        for individual in union:
            calculate_raw_fitness(individual)
        nearest_neighbours = determine_neighbours(union)
        calculate_solution_density(union, nearest_neighbours)
        archive, dominated = get_non_dominated_solutions(union)
        if len(archive) < archive_size:
            archive = populate_archive_with_remaining_best(archive, archive_size, dominated)
        elif len(archive) > archive_size:
            archive = truncate_similar_individuals(archive, archive_size)
        parents = get_parents(archive)
        population = crossover_and_mutate(parents, probability_crossover, probability_mutation, codon_use_table)
    return archive


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


def eval_host(individual, ancestor_sequence):
    assert (individual is SequenceContainer)
    sequence = getattr(individual, "sequence")
    return sum(['host'] for a, b in zip(sequence, ancestor_sequence) if a != b)


def eval_restriction_sites(individual, restrict_sites):
    """
    TODO: Make it remove rest sites
    """
    assert (individual is SequenceContainer)
    sequence = getattr(individual, "sequence")
    # check unwanted restriction sites
    analysis = Analysis(restrictionbatch=restrict_sites, sequence=sequence)
    result = analysis.full()
    # score the sequence based on the number of restriction sites
    score = 0
    for enz, cuts in result.items():
        for cut in cuts:
            score += 1
    return score


def eval_start_sites(individual, ribosome_binding_sites=RibosomeBindingSites, table_name="Standard"):
    """
    TODO: Make it remove start sites
    """
    assert (individual is SequenceContainer)
    sequence = getattr(individual, "sequence")
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


def eval_repeats(individual, window_size=9):
    assert(individual is SequenceContainer)
    sequence = getattr(individual, "sequence")
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


def eval_homopolymers(individual, n_codons=2, homopolymer_threshold=4):
    assert(individual is SequenceContainer)
    sequence = getattr(individual, "sequence")
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


def eval_splice_sites(individual):
    assert(individual is SequenceContainer)
    sequence = getattr(individual, "sequence")

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


def eval_gc_content(individual, gc):
    assert(individual is SequenceContainer)
    sequence = getattr(individual, "sequence")
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


def eval_hairpins(individual, stem_length=10):
    """
    TODO: Make it remove hairpins
    :param individual:
    :param stem_length:
    :return:
    """
    assert(individual is SequenceContainer)
    sequence = getattr(individual, "sequence")
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