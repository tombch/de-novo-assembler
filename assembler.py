import copy
import random
import collections
import matplotlib.pyplot as plt
from timeit import default_timer as timer


def generate_sequence(length, bases='ATCG'):
    '''
    Generate and return a random sequence of the given length, from the given bases (characters).
    '''
    bases = list(bases)
    sequence = ''.join(random.choice(bases) for i in range(length))
    return sequence


def hamming_distance(seq1, seq2):
    '''
    Return the Hamming distance between two sequences.
    If seq1 and seq2 have different lengths, returns None.
    '''
    if len(seq1) != len(seq2) or (not isinstance(seq1, str)) or (not isinstance(seq2, str)):
        return None
    else:
        distance = sum(1 for base1, base2 in zip(seq1, seq2) if base1 != base2)
        return distance


def get_kmers(sequence, k, coverage=1):
    '''
    Return all substrings of length k from the given sequence.
    '''
    kmers = []
    for i in range(coverage):
        kmers += [sequence[i : i + k] for i in range(len(sequence) - k + 1)]
    return kmers


def count_kmers(kmers):
    '''
    Transform list of kmers to a kmer_counts dictionary. 
    Unique kmers are the keys, and their number of appearances in the kmers list are the values. 
    '''
    kmer_counts = {}
    for kmer in kmers:
        kmer_counts[kmer] = 0
    for kmer in kmers:
        kmer_counts[kmer] += 1
    return kmer_counts


def get_unique_counts(kmer_counts):
    '''
    Return all unique kmer counts (and the number of kmers with this count) present in the kmer_counts dictionary.
    '''
    unique_counts_dict = {}
    for kmer, count in kmer_counts.items():
        if not (count in unique_counts_dict):
            unique_counts_dict[count] = 1
        else:
            unique_counts_dict[count] += 1
    unique_counts = []
    for count, num_kmers in unique_counts.items():
        unique_counts.append((count, num_kmers))
    return unique_counts


def list_kmers(kmer_counts):
    '''
    Transform kmer_counts dictionary to a list of kmers.
    '''
    kmers = []
    for kmer, count in kmer_counts.items():
        for i in range(count):
            kmers.append(kmer)
    return kmers
    

def kmer_count_distribution(kmers, threshold=0, plot=False):
    '''
    Given a set of kmers, return a list describing the distribution of kmer count values (the number of times a kmer appears) in the set.
    The list can also be optionally plotted by setting plot=True.
    '''
    unique_kmer_counts = get_unique_counts(count_kmers(kmers))
    kmer_counts_dist = sorted([(count, num_kmers) for (count, num_kmers) in unique_kmer_counts if count >= threshold], key=lambda x: x[0])
    if plot:
        plt.scatter(*zip(*kmer_counts_dist))
        plt.xlabel("kmer count")
        plt.ylabel("# kmers")
        plt.show()
    return kmer_counts_dist


def mutate(kmers, num_mutations, shuffle=False):
    '''
    Mutate a random base of a random kmer, with the number of mutations given by num_mutations.
    '''
    bases = ['A', 'T', 'C', 'G']
    # Do not mutate original kmer list, but create a new list instead
    kmers = list(kmers)
    if shuffle:
        random.shuffle(kmers)
    for i in range(num_mutations):
        kmer_index = random.choice(range(len(kmers)))
        base_index = random.choice(range(len(kmers[kmer_index])))
        kmers[kmer_index] = kmers[kmer_index][0 : base_index] + random.choice(bases) + kmers[kmer_index][base_index + 1:]
    return kmers


def delete(kmers, num_deletions, shuffle=False):
    '''
    Pop a given number of kmers off the end of a list of kmers.
    '''
    # Do not delete from original kmer list, but create a new list instead
    kmers = list(kmers)
    if shuffle:
        random.shuffle(kmers)
    for i in range(num_deletions):
        kmers.pop()
    return kmers


def noisify(kmers, kmer_mutations=0, kmer_deletions=0):
    '''
    Introduce noise to a set of kmers by shuffling, mutating and deleting them.
    '''
    kmers = delete(mutate(kmers, num_mutations=kmer_mutations, shuffle=True), num_deletions=kmer_deletions)
    return kmers


def estimate_coverage(kmer_counts, threshold=0):
    '''
    Given a kmer_counts dictionary and a threshold for noise, estimate the number of sequences present.
    Any kmer that appears less than the threshold amount of times, will be assumed to be an error/noise and will be discarded.
    '''
    count_values = {}
    for kmer, count in kmer_counts.items():    
        count_values[count] = 0
    for kmer, count in kmer_counts.items():    
        count_values[count] += 1    
    max_value = 0
    coverage = 0
    for count, value in count_values.items():
        if value > max_value and count > threshold:
            max_value = value
            coverage = count
    return coverage


def reduce_kmers(kmer_counts, coverage):
    '''
    Given a kmer_counts dictionary and a coverage, return a scaled-down set of kmers that represent one instance of the sequence.
    '''
    kmers = []
    for kmer, count in kmer_counts.items():
        remainder = count % coverage
        if remainder < (coverage - remainder):
            count = int((count - remainder) / coverage)
        else:
            count = int((count + (coverage - remainder)) / coverage)
        for i in range(count):
            kmers.append(kmer)
    return kmers


def de_bruijn_graph(kmers):
    '''
    Return a De Bruijn graph from the given kmers.
    '''
    graph = {}
    k_minus_1 = {}
    for kmer in kmers:
        k = len(kmer)
        k_minus_1[kmer] = get_kmers(kmer, k - 1)
        k_1, k_2 = k_minus_1[kmer]
        # Each node records deques of incoming and outgoing nodes
        graph[k_1] = {'in' : collections.deque([]), 'out' : collections.deque([])}
        graph[k_2] = {'in' : collections.deque([]), 'out' : collections.deque([])}
    for kmer in kmers:
        k_1, k_2 = k_minus_1[kmer]
        graph[k_2]['in'].append(k_1)
        graph[k_1]['out'].append(k_2)
    return graph


def eulerian_path(graph):
    '''
    Use Hierholzer's algorithm to find and return an Eulerian path in the given graph.
    If an Eulerian path does not exist, returns None.
    '''
    edge_count = 0
    start_count = 0
    end_count = 0
    equal_count = 0
    current_node = ""
    # Find valid start node, and determine path existence
    for node, edge_nodes in graph.items():
        edge_count += len(edge_nodes['out'])
        if len(edge_nodes['out']) - len(edge_nodes['in']) == 1:
            current_node = node
            start_count += 1
        elif len(edge_nodes['in']) - len(edge_nodes['out']) == 1:
            end_count += 1
        elif len(edge_nodes['in']) == len(edge_nodes['out']):
            equal_count += 1
    # Return None if an Eulerian path does not exist
    correct_degree_count = start_count + end_count + equal_count
    if start_count > 1 or end_count > 1 or correct_degree_count != len(graph):
        return None
    # If no valid start node was determined, we can choose any node with non-zero degree
    # We look for the first node with an outgoing edge and assign this as the start node
    if not current_node:
        for node, edge_nodes in graph.items():
            if len(edge_nodes['out']) > 0:
                current_node = node
                break
    # Now we traverse the graph to find an Eulerian path.
    path = collections.deque([])
    history = collections.deque([])
    graph = copy.deepcopy(graph)
    while len(path) != edge_count + 1:
        if len(graph[current_node]['out']) == 0:
            if len(history) != 0:
                next_node = history[-1]
                history.pop()
            path.appendleft(current_node)
        else:
            next_node = graph[current_node]['out'][0]
            graph[current_node]['out'].popleft()
            history.append(current_node)
        current_node = next_node
    return path


def join(path):
    '''
    Join together the given Eulerian path to form a sequence.
    '''
    path = list(path)
    sequence = path[0]
    for i, node in enumerate(path[1:], start=1):
        # TODO: why did I do this?
        k = len(path[i - 1])
        sequence += node[k - 1:]
    return sequence


def assemble(kmers, describe=False, single_coverage=False, threshold=0):
    '''
    Assemble a sequence from the given kmers.
    If single_coverage=False, then (before constructing a De Bruijn graph), the coverage is estimated, and the set of kmers are reduced to a coverage of 1.
    '''
    if describe:
        start = timer()
        if not single_coverage:
            print('Finding kmer counts...', end=' ', flush=True)
            kmer_counts = count_kmers(kmers)
            print('done.')
            print('Estimating coverage...', end=' ', flush=True)
            coverage_estimate = estimate_coverage(kmer_counts, threshold=threshold)
            print('done.')
            print(f'Estimated coverage: {coverage_estimate}')
            if coverage_estimate != 1:
                print('Reducing coverage of kmers to 1...', end=' ', flush=True)
                kmers = reduce_kmers(kmer_counts, coverage_estimate)
                print('done.')
        print('Generating De Bruijn graph...', end=' ', flush=True)
        graph = de_bruijn_graph(kmers)
        print('done.')
        print('Finding Eulerian path...', end=' ', flush=True)
        path = eulerian_path(graph)
        print('done.')
        print('Joining path to form sequence...', end=' ', flush=True)
        sequence = join(path)
        print('done.')
        end = timer()
        print(f'Time taken: {end - start} seconds')
    else:
        if not single_coverage:
            kmer_counts = count_kmers(kmers)
            kmers = reduce_kmers(kmer_counts, estimate_coverage(kmer_counts, threshold=threshold))
        sequence = join(eulerian_path(de_bruijn_graph(kmers)))
    return sequence