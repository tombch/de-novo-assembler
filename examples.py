import assembler as ar
import random
from timeit import default_timer as timer


def main():
    # Example 1
    message = 'a_galaxy_far__far___far____far_____away...'
    kmers = ar.get_kmers(message, 10)
    random.shuffle(kmers)
    start = timer()
    assembled_message = ar.assemble(kmers, describe=True)
    end = timer()
    print(f'time taken: {end - start}')
    distance = ar.hamming_distance(message, assembled_message)
    print(f'dist(original, assembled) = {distance}')
    
    print('')

    # Example 2
    sequence = ar.generate_sequence(50000)
    coverage = 100
    kmers = []
    # Simulate noisy kmers
    for i in range(coverage):
        kmers += ar.mutate(ar.get_kmers(sequence, 50), num_mutations=3000)
    random.shuffle(kmers)
    ar.kmer_count_distribution(kmers, threshold=10, plot=True)
    start = timer()
    kmers = ar.estimate_coverage_reduce_kmers(kmers, threshold=10, describe=True)
    assembled_sequence = ar.assemble(kmers, describe=True)
    end = timer()
    print(f'time taken: {end - start}')
    distance = ar.hamming_distance(sequence, assembled_sequence)
    print(f'dist(original, assembled) = {distance}')
    
    print('')

    # Example 3
    with open("genome.txt") as genome_fh:
        genome = genome_fh.readline()
        clean_kmers = ar.get_kmers(genome, 300)
        kmers = []
        coverage = 100
        for i in range(coverage):
            kmers += ar.mutate(clean_kmers, num_mutations=3000)
        random.shuffle(kmers)
        ar.kmer_count_distribution(kmers, threshold=10, plot=True)
        start = timer()
        kmers = ar.estimate_coverage_reduce_kmers(kmers, threshold=10, describe=True)
        assembled_genome = ar.assemble(kmers, describe=True)
        end = timer()
        print(f'time taken: {end - start}')
        distance = ar.hamming_distance(genome, assembled_genome)
        print(f'dist(original, assembled) = {distance}')


if __name__ == '__main__':
    main()