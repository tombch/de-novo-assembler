import assembler as ar
import random
from timeit import default_timer as timer


def main():
    # Example 1
    start = timer()
    message = 'a_galaxy_far_far_far_far_away...'
    kmers = ar.get_kmers(message, 3)
    random.shuffle(kmers)
    assembled_message = ar.assemble(kmers)
    distance = ar.hamming_distance(message, assembled_message)
    end = timer()
    print(f'dist(original, assembled) = {distance}')
    print(f'time taken: {end - start}')

    # Example 2
    start = timer()
    sequence = ar.generate_sequence(100000)
    coverage = 100
    kmers = []
    # Simulate kmers drawn from reads
    for i in range(coverage):
        kmers += ar.mutate(ar.get_kmers(sequence, 100), 10000)
    kmer_count = ar.count_kmers(kmers)
    threshold = 10
    coverage = ar.estimate_coverage(kmer_count, threshold)
    kmers = ar.reduce_kmers(kmer_count, coverage)
    random.shuffle(kmers)
    assembled_sequence = ar.assemble(kmers, describe=True)
    distance = ar.hamming_distance(sequence, assembled_sequence)
    end = timer()
    print(f'dist(original, assembled) = {distance}')
    print(f'time taken: {end - start}')


if __name__ == '__main__':
    main()