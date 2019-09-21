'''
Initialize population P
Create empty external set E
while generations
    compute fitness of all individuals in population and external set
    copy all individuals evaluating to nondominated vectors in P & E to E
    Remove elements from E if necessary using truncation
    if E not full, fill with dominated individuals from P
    binary tournament selection with replacement to fill mating pool M
    crossover and mutations in M
'''
import logging
import random
import re
import copy
import math
import numpy
import operator
import sys
from itertools import product
import multiprocessing
from Bio import Seq
from Bio.Alphabet import IUPAC
from Bio.Data import CodonTable
from Bio.Restriction import Analysis
from Bio.SeqUtils import CodonUsage, GC, seq3

from test_functions import codon_use_table
from . import codon_use, seq, MutableSeq

max_generations = 2


def crossover_and_mutate(parents, probability_crossover, probability_mutation):
    pass


def get_parents(archive):
    pass


def truncate_similar_individuals(archive):
    pass


def populate_archive_with_remaining_best(archive):
    pass


def get_non_dominated_solutions(union):
    pass


def calculate_solution_density(individual):
    pass


def calculate_raw_fitness(individual):
    pass


def calculate_fitness(individual, problem_size):
    pass


def initialize_population(population_size, problem_size, ancestor_sequence, probability_mutation):
    """
    Initializes a population of sequences based on a (hopefully) codon optimized
    sequence
    :param probability_mutation: float determines probability of mutation
    :param population_size: int determines minimum population size
    :param problem_size: int helps determines minimum population size
    :param ancestor_sequence: Bio.seq.seq codon optimized sequence
    :return:
    """
    population = {}
    mutable_seq = ancestor_sequence.tomutable()
    size = population_size if (problem_size < population_size) else size = problem_size
    for idx in range(size):
        population[mutate_sequence(mutable_seq.to_seq(), probability_mutation)] = sys.maxsize
    return population


def mutate_sequence(sequence, mutation_probability=0.05, offset=0):
    """
    Takes a single sequence and gives it a random number of mutations.
    :param sequence:
    :param mutation_probability:
    :param offset:
    :return:
    """
    mutable_seq = sequence.tomutable()
    num_codons = len(mutable_seq) // 3
    num_mutations = math.ceil(num_codons * mutation_probability)
    for _ in range(num_mutations):
        position = 3 * random.randrange(0, len(mutable_seq) // 3)
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
    amino_acid = seq3(CodonTable.standard_dna_table.forward_table[str(codon_in)]).upper()
    synonymous_codons, codon_use_freq = codon_use_table[amino_acid]
    if len(synonymous_codons) == 1:
        return codon_in

    # pick new codon
    codon_out = codon_in
    while codon_in == codon_out:
        codon_out = random.choices(synonymous_codons, codon_use_freq).pop()

    return codon_out


def optimize_strength_pareto_evolutionary_algorithm(population_size, archive_size, problem_size, probability_crossover, probability_mutation, ancestor_sequence):
    population = initialize_population(population_size, problem_size, ancestor_sequence, probability_mutation)
    archive = []
    for generation in range(0, max_generations):
        for individual in population:
            calculate_fitness(individual, problem_size)
        union = population + archive
        for individual in union:
            calculate_raw_fitness(individual)
            calculate_solution_density(individual)
        archive += get_non_dominated_solutions(union)
        if len(archive) < archive_size:
            populate_archive_with_remaining_best(archive)
        elif len(archive) > archive_size:
            truncate_similar_individuals(archive)
        parents = get_parents(archive)
        population = crossover_and_mutate(parents, probability_crossover, probability_mutation)
    return archive
