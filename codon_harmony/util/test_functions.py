import copy
import math
import random
from Bio import Seq
from Bio.Alphabet import IUPAC


def generate_initial_population(pop_size=10):
    """Create the population of mutated sequences"""
    mutants = {}
    mutable_seq = dna_sequence.tomutable()
    for idx in range(pop_size):
        new_seq = mutable_seq
        new_seq = new_seq.toseq()
        mutants[mutate_sequence(new_seq, idx)]= 0
    return mutants


def mutate_sequence(sequence, mutation_probability=0.0005, offset=0):
    """
    Takes a single sequence and gives it a random number of mutations.
    :param sequence:
    :param mutation_probability:
    :param offset:
    :return:
    """
    mutable_seq = sequence.tomutable()
    num_codons = len(mutable_seq)//3
    num_mutations = math.ceil(num_codons * mutation_probability)
    for _ in range(num_mutations):
        position = random.randrange(0, len(mutable_seq) // 3)
        codon_idx = slice(offset + position, (offset + 3) + position)
        new_codon = mutate_codon(mutable_seq[codon_idx])
        mutable_seq[codon_idx] = new_codon

    return mutable_seq.toseq()


def mutate_codon(codon_in):
    """Select a synonymous codon in accordance with the frequency of use
    in the host organism.

    Args:
        codon_in (Bio.Seq.Seq): A single codon.

    Returns:
        Bio.Seq.Seq: A new codon.
    """
    # pick new codon
    codon_out = codon_in
    while codon_in == codon_out:
        codon_out = "YEE"

    return codon_out


dna_sequence = Seq.Seq("AATTCCGGATCG", IUPAC.ambiguous_dna)
print(generate_initial_population())
