import random
import collections


def generate_sequence(length):
    '''
    Generate and return a random sequence of the given length.
    '''
    bases = ['A', 'T', 'C', 'G']
    sequence = ''.join(random.choice(bases) for i in range(length))
    return sequence


def hamming_distance(seq1, seq2):
    '''
    Return the Hamming distance between two sequences.
    If seq1 and seq2 have different length, returns None.
    '''
    if len(seq1) != len(seq2):
        return None
    else:
        distance = sum(1 for base1, base2 in zip(seq1, seq2) if base1 != base2)
        return distance


def get_kmers(sequence, k):
    '''
    Get all substrings of length k from the given sequence.
    '''
    kmers = [sequence[i : i + k] for i in range(len(sequence) - k + 1)]
    return kmers


def count_kmers(kmers):
    '''
    Transform list of kmers to a kmer_count dictionary. 
    Unique kmers are the keys, and their number of appearances in the kmers list are the values. 
    '''
    kmer_count = {}
    for kmer in kmers:
        kmer_count[kmer] = 0
    for kmer in kmers:
        kmer_count[kmer] += 1
    return kmer_count


def get_unique_counts(kmer_count):
    '''
    Return all unique count values present in the kmer_count dictionary.
    '''
    unique_counts = []
    for kmer, count in kmer_count.items():
        if not (count in unique_counts):
            unique_counts.append(count)
    return unique_counts


def list_kmers(kmer_count):
    '''
    Transform kmer_count dictionary to a list of kmers.
    '''
    kmers = []
    for kmer, count in kmer_count.items():
        for i in range(count):
            kmers.append(kmer)
    return kmers


def mutate(kmers, num_mutations):
    '''
    Mutate a random base of a random kmer, with the number of mutations given by num_mutations.
    '''
    bases = ['A', 'T', 'C', 'G']
    for i in range(num_mutations):
        kmer_index = random.choice(range(len(kmers)))
        base_index = random.choice(range(len(kmers[kmer_index])))
        kmers[kmer_index] = kmers[kmer_index][0 : base_index] + random.choice(bases) + kmers[kmer_index][base_index + 1:]
    return kmers


def estimate_coverage(kmer_count, threshold):
    '''
    Given a kmer_count dictionary and a threshold for noise, estimate the number of sequences present.
    Any kmer that appears less than the threshold amount of times, will be assumed to be an error/noise and will be discarded.
    '''
    count_values = {}
    for kmer, count in kmer_count.items():    
        count_values[count] = 0
    for kmer, count in kmer_count.items():    
        count_values[count] += 1    
    max_value = 0
    coverage = 0
    for count, value in count_values.items():
        if value > max_value and count > threshold:
            max_value = value
            coverage = count
    return coverage


def reduce_kmers(kmer_count, coverage):
    '''
    Given a kmer_count dictionary and a coverage, return a scaled-down set of kmers that represent one instance of the sequence.
    '''
    kmers = []
    for kmer, count in kmer_count.items():
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
    Get a De Bruijn graph from the given kmers.
    '''
    graph = {}
    for kmer in kmers:
        k = len(kmer)
        k_m1, k_m2 = get_kmers(kmer, k - 1)
        graph[k_m1] = {'in' : collections.deque([]), 'out' : collections.deque([])} # Incoming edge nodes, outgoing edge nodes
        graph[k_m2] = {'in' : collections.deque([]), 'out' : collections.deque([])}
    for kmer in kmers:
        k = len(kmer)
        k_m1, k_m2 = get_kmers(kmer, k - 1)
        graph[k_m2]['in'].append(k_m1)
        graph[k_m1]['out'].append(k_m2)
    return graph


def eulerian_path(graph):
    '''
    Use Hierholzer's algorithm to find and return a Eulerian path in the given graph.
    If a Eulerian path does not exist, returns None.
    '''
    path = collections.deque()
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
    # Return None if a Eulerian path does not exist
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
    history = collections.deque()
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
    sequence = path[0]
    for i, node in enumerate(path):
        if i != 0:
            k = len(path[i - 1])
            sequence += node[k - 1:]
    return sequence


def assemble(kmers, describe=False):
    '''
    Assemble a sequence from the given kmers.
    '''
    if describe:
        print('Generating De Bruijn graph...', end=' ', flush=True)
        graph = de_bruijn_graph(kmers)
        print('done.')
        print('Finding Eulerian path...', end=' ', flush=True)
        path = eulerian_path(graph)
        print('done.')
        print('Joining sequence...', end=' ', flush=True)
        sequence = join(path)
        print('done.')
    else:
        sequence = join(eulerian_path(de_bruijn_graph(kmers)))
    return sequence