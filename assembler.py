import random


def get_kmers(sequence, k):
    '''
    Get all substrings of length k from the given sequence.
    '''
    kmers = [sequence[i : i + k] for i in range(len(sequence) - k + 1)]
    return kmers


def de_bruijn_graph(kmers):
    '''
    Get a De Bruijn graph from the given kmers.
    '''
    graph = {}
    for kmer in kmers:
        k = len(kmer)
        k_m1, k_m2 = get_kmers(kmer, k - 1)
        graph[k_m1] = {'in' : [], 'out' : []} # Incoming edge nodes, outgoing edge nodes
        graph[k_m2] = {'in' : [], 'out' : []}
    for kmer in kmers:
        k = len(kmer)
        k_m1, k_m2 = get_kmers(kmer, k - 1)
        graph[k_m2]['in'].append(k_m1)
        graph[k_m1]['out'].append(k_m2)
    return graph


def get_eulerian_path(graph):
    '''
    Use Hierholzer's algorithm to find and return a Eulerian path in the given graph.
    If a Eulerian path does not exist, returns None.
    '''
    path = []
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
    history = []
    while len(path) != edge_count + 1:
        if len(graph[current_node]['out']) == 0:
            if len(history) > 0:
                next_node = history[-1]
                history.pop(-1)
                for i, out_node in enumerate(graph[next_node]['out']):
                    if out_node == current_node:
                        graph[next_node]['out'].pop(i)
                        break
            path = [current_node] + path
        else:
            next_node = graph[current_node]['out'][0]
            graph[current_node]['out'].pop(0)
            history.append(current_node)
        current_node = next_node
    return path


def assemble(path):
    '''
    Assemble a sequence from the given Eulerian path.
    '''
    sequence = path[0]
    for i, node in enumerate(path):
        if i != 0:
            k = len(path[i - 1])
            sequence += node[k - 1:]
    return sequence


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
    seq1_length = len(seq1)
    if seq1_length != len(seq2):
        return None
    else:
        mismatch_generator = (1 for i in range(seq1_length) if seq1[i] != seq2[i])
        distance = sum(x for x in mismatch_generator)    
        return distance


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


def count_kmers(kmers):
    '''
    Transform the kmers from a direct list of kmers to a dictionary representation. 
    Unique kmers are the keys, and their number of appearances in the kmer list are the values. 
    '''
    kmer_count = {}
    for kmer in kmers:
        kmer_count[kmer] = 0
    for kmer in kmers:
        kmer_count[kmer] += 1
    return kmer_count


def display_unique_counts(kmer_count):
    unique_counts = []
    for kmer, count in kmer_count.items():
        if not (count in unique_counts):
            unique_counts.append(count)
    return unique_counts


def estimate_num_sequences(kmer_count, threshold):
    count_values = {}
    for kmer, count in kmer_count.items():    
        count_values[count] = 0
    for kmer, count in kmer_count.items():    
        count_values[count] += 1    
    max_value = 0
    num_sequences = 0
    for count, value in count_values.items():
        if value > max_value and count > threshold:
            max_value = value
            num_sequences = count
    return num_sequences


def reduce_kmers(kmer_count, num_sequences):
    kmers = []
    for kmer, count in kmer_count.items():
        remainder = count % num_sequences
        if remainder < (num_sequences - remainder):
            count = int((count - remainder) / num_sequences)
        else:
            count = int((count + (num_sequences - remainder)) / num_sequences)
        for i in range(count):
            kmers.append(kmer)
    return kmers


message = 'a_galaxy_far_far_far_far_away...'
assembled_message = assemble(get_eulerian_path(de_bruijn_graph(get_kmers(message, 3))))
print(f'dist(original, assembled) = {hamming_distance(message, assembled_message)}')

sequence = generate_sequence(30000)
num_sequences = 100
kmers = []
# Simulate kmers drawn from reads
for i in range(num_sequences):
   kmers += mutate(get_kmers(sequence, 100), 5000)
kmer_count = count_kmers(kmers)
threshold = 10
num_sequences = estimate_num_sequences(kmer_count, threshold)
kmers = reduce_kmers(kmer_count, num_sequences)
random.shuffle(kmers)
assembled_sequence = assemble(get_eulerian_path(de_bruijn_graph(kmers)))
print(f'dist(original, assembled) = {hamming_distance(sequence, assembled_sequence)}')