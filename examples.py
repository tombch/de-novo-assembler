import assembler as ar
import random
from timeit import default_timer as timer


def main():
    # Example 1
    message = 'a_galaxy_far_far_far_far_away...'
    kmers = ar.get_kmers(message, 3)
    random.shuffle(kmers)
    start = timer()
    assembled_message = ar.assemble(kmers)
    end = timer()
    print(f'assembly time: {end - start}')
    distance = ar.hamming_distance(message, assembled_message)
    print(f'dist(original, assembled) = {distance}')

    # Example 2
    sequence = ar.generate_sequence(100000)
    coverage = 100
    kmers = []
    # Simulate noisy kmers
    for i in range(coverage):
        kmers += ar.mutate(ar.get_kmers(sequence, 100), num_mutations=10000)
    kmer_count = ar.count_kmers(kmers)
    coverage_estimate = ar.estimate_coverage(kmer_count, threshold=10)
    kmers = ar.reduce_kmers(kmer_count, coverage_estimate)
    random.shuffle(kmers)
    start = timer()
    assembled_sequence = ar.assemble(kmers, describe=True)
    end = timer()
    print(f'assembly time: {end - start}')
    distance = ar.hamming_distance(sequence, assembled_sequence)
    print(f'dist(original, assembled) = {distance}')


if __name__ == '__main__':
    main()