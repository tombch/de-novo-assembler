#include <iostream>
#include <string>
#include <random>
#include <unordered_map>
#include <deque>
#include <chrono>
#include <limits>
#include <algorithm>


typedef std::unordered_map<std::string, std::pair<std::deque<std::string>, std::deque<std::string>>> directedGraph;
typedef std::chrono::_V2::system_clock::time_point timePoint;


float duration(timePoint start, timePoint end)
{
    return float(std::chrono::duration_cast<std::chrono::microseconds>(end - start).count()) / 1000000;
}


std::string genSequence(int length, bool time = false) 
{
    /* 
    Generate and return a random sequence (of the given length) of bases.
    */
    timePoint start = std::chrono::high_resolution_clock::now();
    std::string bases = "ATCG";
    std::random_device rd;
    std::mt19937 rng(rd());
    std::uniform_int_distribution<int> uniformDist(0, 3);
    std::string sequence = "";
    sequence.reserve(length);
    int randomInteger;
    for (int i = 0; i < length; i++) 
    {
        randomInteger = uniformDist(rng);
        sequence += bases[randomInteger];
    }
    timePoint end = std::chrono::high_resolution_clock::now();
    if (time == true)
        std::cout << "genSequence: " << duration(start, end) << " seconds" << std::endl; 
    return sequence;
}


int hammingDistance(std::string &seq1, std::string &seq2, bool time = false)
{
    /*
    Return the Hamming distance between two sequences.
    */
    timePoint start = std::chrono::high_resolution_clock::now();
    int distance;
    if (seq1.length() != seq2.length()) 
        distance = -1;
    else 
    {
        distance = 0;
        int seqLength = seq1.length();
        for (int i = 0; i < seqLength; i++) 
        {
            if (seq1[i] != seq2[i]) 
                distance += 1;
        }
    }
    timePoint end = std::chrono::high_resolution_clock::now();
    if (time == true)
        std::cout << "hammingDistance: " << duration(start, end) << " seconds" << std::endl; 
    return distance;
}


std::vector<std::string> genKmers(std::string &sequence, int k, bool time = false)
{
    /*
    Return all substrings of length k from the given sequence.
    */
    timePoint start = std::chrono::high_resolution_clock::now();
    std::vector<std::string> kmers;
    if ((sequence.length() == 0) || (k == 0) || (sequence.length() < k))
        kmers.push_back("");
    else
    {
        int numKmers = sequence.length() - k + 1;
        kmers.reserve(numKmers);
        std::string kmer = sequence.substr(0, k);
        kmers.push_back(kmer);
        int kmerEndIndex = k;
        for (int i = 1; i < numKmers; i++)
        {
            kmer = kmers.back();
            kmer.erase(0, 1);
            kmer += sequence[kmerEndIndex];
            kmers.push_back(kmer);
            kmerEndIndex += 1;
        }
    }
    timePoint end = std::chrono::high_resolution_clock::now();
    if (time == true)
        std::cout << "genKmers: " << duration(start, end) << " seconds" << std::endl; 
    return kmers;
}


void shuffleKmers(std::vector<std::string> &kmers, bool time = false)
{
    /*
    Shuffle a set of kmers.
    */
    timePoint start = std::chrono::high_resolution_clock::now();
    auto rng = std::default_random_engine {};
    std::shuffle(std::begin(kmers), std::end(kmers), rng);
    timePoint end = std::chrono::high_resolution_clock::now();
    if (time == true)
        std::cout << "shuffleKmers: " << duration(start, end) << " seconds" << std::endl; 
}


directedGraph deBruijnGraph(std::vector<std::string> &kmers, bool time = false)
{
    /*
    Return a De Bruijn graph from the given kmers.
    */
    timePoint start = std::chrono::high_resolution_clock::now();
    directedGraph graph;
    if (kmers.size() == 0)
    {
        graph[""].first.push_back("");
        graph[""].second.push_back("");
    }
    else
    {
        int kMinus1 = kmers[0].length() - 1;
        int numKmers = kmers.size();
        for (int i = 0; i < numKmers; i++)
        {
            std::string kMinus1Left = kmers[i].substr(0, kMinus1);
            std::string kMinus1Right = kmers[i].substr(1, kMinus1);
            graph[kMinus1Right].first.push_back(kMinus1Left);
            graph[kMinus1Left].second.push_back(kMinus1Right);
        }
    }
    timePoint end = std::chrono::high_resolution_clock::now();
    if (time == true)
        std::cout << "deBruijnGraph: " << duration(start, end) << " seconds" << std::endl; 
    return graph;
}


std::deque<std::string> eulerianPath(directedGraph &graph, bool time = false)
{
    /*
    Use Hierholzer's algorithm to find and return an Eulerian path in the given graph.
    If an Eulerian path does not exist, returns an empty path.
    */
    timePoint start = std::chrono::high_resolution_clock::now();
    int edgeCount = 0;
    int startCount = 0;
    int endCount = 0;
    int nonZeroEqualCount = 0;
    std::string currentNode = "";
    std::deque<std::string> path;
    std::unordered_map<std::string, int> outIndexes;
    for (auto &node : graph)
    {
        auto &nodeVal = node.first;
        int inNodesSize = node.second.first.size();
        int outNodesSize = node.second.second.size();
        edgeCount += outNodesSize;
        if (outNodesSize - inNodesSize == 1)
        {
            if (currentNode == "")
                currentNode = nodeVal;
            startCount += 1;
        }
        else if (inNodesSize - outNodesSize == 1)
            endCount += 1;
        else if ((inNodesSize == outNodesSize) && inNodesSize != 0)
            nonZeroEqualCount += 1;
        outIndexes[nodeVal] = 0;
    }
    if ((startCount > 1) || (endCount > 1) || (startCount + endCount + nonZeroEqualCount != graph.size()))
        path.push_back("");
    else
    {
        if (currentNode == "")
        {
            for (auto &node : graph)
            {
                if (node.second.second.size() > 0)
                {
                    currentNode = node.first;
                    break;
                }
            }
        }
        std::deque<std::string> history;
        std::string nextNode;
        int edgeCountPlusOne = edgeCount + 1;
        while (path.size() != edgeCountPlusOne)
        {
            auto &currentOutNodes = graph[currentNode].second;
            auto &currentOutIndex = outIndexes[currentNode];
            if (currentOutIndex == currentOutNodes.size())
            {
                if (history.size() != 0)
                {
                    nextNode = history.back();
                    history.pop_back();
                }
                path.push_front(currentNode);
            }
            else
            {
                nextNode = currentOutNodes[currentOutIndex];
                currentOutIndex += 1;
                history.push_back(currentNode);
            }
            currentNode = nextNode;
        }
    }
    timePoint end = std::chrono::high_resolution_clock::now();
    if (time == true)
        std::cout << "eulerianPath: " << duration(start, end) << " seconds" << std::endl; 
    return path;
}


std::string join(std::deque<std::string> &path, bool time = false)
{
    /*
    Join together the given Eulerian path to form a sequence.
    */
    timePoint start = std::chrono::high_resolution_clock::now();
    int pathSize = path.size();
    std::string sequence;
    if (pathSize == 0)
        sequence = "";
    else
    {
        sequence = path[0];
        for (int i = 1; i < pathSize; i++)
            sequence += path[i].back();
    }
    timePoint end = std::chrono::high_resolution_clock::now();
    if (time == true)
        std::cout << "join: " << duration(start, end) << " seconds" << std::endl; 
    return sequence;
}


// Work in progress
int estimateCoverage(std::vector<std::string> &kmers, int threshold, bool time = false)
{
    /*
    Estimate the number of sequences present in a set of kmers.
    Any kmer that appears less than the threshold amount of times in the set will be assumed to be noise.
    */
    timePoint start = std::chrono::high_resolution_clock::now();
    int coverage;
    if (kmers.size() == 0)
        coverage = -1;
    else
    {
        std::unordered_map<std::string, int> kmerFrequencies;
        for (auto &kmer : kmers)
        {
            if (kmerFrequencies.find(kmer) == kmerFrequencies.end())
                kmerFrequencies[kmer] = 1;
            else
                kmerFrequencies[kmer] += 1;
        }
        std::unordered_map<int, int> kmerFrequencyCounts;
        for (auto &pair : kmerFrequencies)
        {
            int &frequency = pair.second;
            if (kmerFrequencyCounts.find(frequency) == kmerFrequencyCounts.end())
                kmerFrequencyCounts[frequency] = 1;
            else
                kmerFrequencyCounts[frequency] += 1;
        }
        coverage = std::numeric_limits<int>::min();
        for (auto &pair : kmerFrequencyCounts)
        {
            auto &frequency = pair.first;
            auto &count = pair.second; 
            if ((count > coverage) && (frequency > threshold))
                coverage = frequency;
        }
    }
    timePoint end = std::chrono::high_resolution_clock::now();
    if (time == true)
        std::cout << "estimateCoverage: " << duration(start, end) << " seconds" << std::endl; 
    return coverage;
}


// Work in progress
std::vector<std::string> normaliseKmers(std::vector<std::string> &kmers, bool time = false)
{
    timePoint start = std::chrono::high_resolution_clock::now();
    timePoint end = std::chrono::high_resolution_clock::now();
    if (time == true)
        std::cout << "normaliseKmers: " << duration(start, end) << " seconds" << std::endl;
    return kmers;
}


// Work in progress
std::string assemble(std::vector<std::string> &reads, bool time = false)
{
    /*
    Assemble a sequence from the given reads.
    */
    timePoint start = std::chrono::high_resolution_clock::now();
    int numReads = reads.size();
    if (numReads == 0)
        return "";
    std::vector<std::string> kmers;
    int k = std::numeric_limits<int>::max();
    int readLength;
    for (int i = 0; i < numReads; i++)
    {
        readLength = reads[i].length();
        if (readLength < k)
            k = readLength - 1;
    }
    for (int i = 0; i < numReads; i++)
    {
        std::vector<std::string> readKmers = genKmers(reads[i], k);
        int readKmersSize = readKmers.size();
        for (int j = 0; j < readKmersSize; j++)
            kmers.push_back(readKmers[j]);
    }
    int coverage = estimateCoverage(kmers, 0);
    if (coverage > 1)
        kmers = normaliseKmers(kmers);
    auto graph = deBruijnGraph(kmers);
    auto path = eulerianPath(graph);
    auto sequence = join(path);
    timePoint end = std::chrono::high_resolution_clock::now();
    if (time == true)
        std::cout << "assemble: " << duration(start, end) << " seconds" << std::endl; 
    return sequence;
}


int main()
{
    timePoint start = std::chrono::high_resolution_clock::now();
    std::string seq = genSequence(1000000, true);
    std::vector<std::string> kmers = genKmers(seq, 50, true);
    // kmers.reserve(2 * kmers.size()); // Don't understand why
    // kmers.insert(kmers.end(), kmers.begin(), kmers.end());
    shuffleKmers(kmers);   
    int coverage = estimateCoverage(kmers, 0, true);
    std::cout << "coverage: " << coverage << std::endl;
    // Sequence assembly
    auto graph = deBruijnGraph(kmers, true);
    auto path = eulerianPath(graph, true);
    auto assembledSeq = join(path, true);
    auto dist = hammingDistance(seq, assembledSeq, true);
    std::cout << "dist: " << dist << std::endl;
    timePoint end = std::chrono::high_resolution_clock::now();
    std::cout << "main: " << duration(start, end) << " seconds" << std::endl;
}