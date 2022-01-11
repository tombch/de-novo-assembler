import assembler as ar
from timeit import default_timer as timer


def main():
    # Example 0
    start = timer()
    sequence = ar.generate_sequence(30000)
    kmers = ar.get_kmers(sequence, k=100)
    graph = ar.de_bruijn_graph(kmers)
    path = ar.eulerian_path(graph)
    assembled_sequence = ar.join(path)
    dist = ar.hamming_distance(sequence, assembled_sequence)
    end = timer()
    print(f'dist(original, assembled) = {dist}')
    print(f'Time taken: {end - start} seconds')

    # Example 1
    # message = 'a_galaxy_far_faar_far_far_away...'
    # kmers = ar.get_kmers(message, 5)
    # kmers = ar.noisify(kmers)
    # print(ar.de_bruijn_graph(kmers))
    # assembled_message = ar.assemble(kmers, describe=True)
    # distance = ar.hamming_distance(message, assembled_message)
    # print(f'dist(original, assembled) = {distance}')
    # print(message)
    # print(assembled_message)

    # Example 2
    # sequence = ar.generate_sequence(50000, bases='atcg')
    # kmers = []
    # # Simulate high coverage
    # kmers = ar.get_kmers(sequence, k=50, coverage=100)
    # # Introduce noise to kmers
    # kmers = ar.noisify(kmers, kmer_mutations=100000, kmer_deletions=100000)
    # # Show distribution of kmer counts to indicate coverage
    # # ar.kmer_count_distribution(kmers, threshold=10, plot=True)
    # # Assemble sequence
    # assembled_sequence = ar.assemble(kmers, describe=True, threshold=10)
    # # Find string distance
    # distance = ar.hamming_distance(sequence, assembled_sequence)
    # print(f'dist(original, assembled) = {distance}')

    # Example 3
    # with open("genome.txt") as genome_fh:
    #     genome = genome_fh.readline()
    #     kmers = []
    #     # Simulate high coverage
    #     kmers = ar.get_kmers(genome, k=300, coverage=100)
    #     # Introduce noise to kmers
    #     kmers = ar.noisify(kmers, kmer_mutations=100000, kmer_deletions=100000)
    #     # Show distribution of kmer counts to indicate coverage
    #     # ar.kmer_count_distribution(kmers, threshold=10, plot=True)
    #     # Assemble sequence
    #     assembled_genome = ar.assemble(kmers, describe=True, threshold=10)
    #     # Find string distance
    #     distance = ar.hamming_distance(genome, assembled_genome)
    #     print(f'dist(original, assembled) = {distance}')


if __name__ == '__main__':
    main()