"""Find Kmer in sequence"""
import os
import argparse
import statistics as st
import networkx as nx
import random


def read_fastq(fastq_file):
    """Read fastq and return list of sequence
    Takes in arguments a FATSQ file"""
    fq_file = open(fastq_file, 'r')
    for _ in fq_file:
        yield next(fq_file).strip()
        next(fq_file)
        next(fq_file)


def cut_kmer(seq, kmer_size):
    """Return list of Kmer
    Takes in arguments a sequence and a Kmer size"""
    for i in range(len(seq)-kmer_size+1):
        yield seq[i:i+kmer_size]


def build_kmer_dict(fastq_file, kmer_size):
    """Return dictionnary of Kmer and occurency of the Kmer
    Takes in arguments a FATSQ file and a Kmer size"""
    kmer_dict = {}
    for i in read_fastq(fastq_file):
        for kmer in cut_kmer(i, kmer_size):
            if not kmer in kmer_dict:
                kmer_dict[kmer] = 1
            else:
                kmer_dict[kmer] += 1
    return kmer_dict


def build_graph(kmer_dict):
    """Build tree of kmer prefixes and suffixes"""
    graph = nx.DiGraph()
    for kmer, val in kmer_dict.items():
        graph.add_edge(kmer[:-1], kmer[1:], weight=val)
    return graph


def get_starting_nodes(graph):
    """Prend un graphe et retourne une liste de noeuds d'entrés"""
    entres = []
    for node in graph.nodes:
        if len(list(graph.predecessors(node))) == 0:
            entres.append(node)
    return entres


def get_sink_nodes(graph):
    """Prend un graphe et retourne une liste de noeuds de sortis"""
    sortis = []
    for node in graph.nodes:
        if len(list(graph.successors(node))) == 0:
            sortis.append(node)
    return sortis


def get_contigs(graph, entres, sortis):
    """Retourne liste de tuple(contig, taille du contig)"""
    all_paths = []
    for node in entres:
        for node2 in sortis:
            path = list(nx.all_simple_paths(graph, node, node2))
            if len(path) > 0:
                all_paths.append(path)
    tupple = []
    for i in range(len(all_paths)):
        for k in range(len(all_paths[i])):
            contig = str(all_paths[i][k][0])
            for j in range(1, len(all_paths[i][k])):
                string = str(all_paths[i][k][j][-1])
                contig += string[-1]
            tupple.append([contig, len(contig)])
    return tupple


def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))


def save_contigs(tupple, file_name):
    """Sort un fichier avec la liste tupple"""
    file = open(file_name, 'w+')
    for i in range(len(tupple)):
        file.write('>contig_' + str(i) + ' len=' + str(tupple[i][1]) + '\n' + str(fill(tupple[i][0])) + '\n')
    file.close()


def std(list_val):
    """Take list of values and return standard deviation"""
    return st.stdev(list_val)


def path_average_weight(graph, path):
    """Take a graph and a path and return average weigth"""
    new_G = graph.subgraph(path)
    wei = []
    for arretes in new_G.edges(data=True):
        wei.append(arretes[2]['weight'])
    mean_wei = st.mean(wei)
    return mean_wei


def remove_paths(graph, path, delete_entry_node, delete_sink_node):
    """Take graph and path and remove nodes
    delete_entry_node pour indiquer si les noeuds d’entrée seront supprimés
    et delete_sink_node pour indiquer si les noeuds de sortie seront supprimés
    Return clean graph"""
    new_G = graph
    for i in range(len(path)):
        new_G.remove_nodes_from(path[i][1:-1])
        if delete_entry_node == True:
            new_G.remove_node(path[i][0])
        if delete_sink_node == True:
            new_G.remove_node(path[i][-1])
    return new_G


def select_best_path(graph, path, path_len, path_wei,
                     delete_entry_node=False, delete_sink_node=False):
    """Take graph, path, path length, path weigth and delete entry or sink node
    Function return clean graph
    By default delete entry or sink node is set as False"""
    new_G = graph
    random.seed(9001)
    best_wei_indice = []
    for i,e in enumerate(path_wei):
        if e == max(path_wei):
            best_wei_indice.append(i)
    best_len_indice = []
    for i,e in enumerate(path_len):
        if e == max(path_len):
            best_len_indice.append(i)
    if len(best_wei_indice) == 1:
        path.remove(path[best_wei_indice[0]])
        new_G = remove_paths(new_G, path, delete_entry_node, delete_sink_node)
    elif len(best_len_indice) == 1:
        path.remove(path[best_len_indice[0]])
        new_G = remove_paths(new_G, path, delete_entry_node, delete_sink_node)
    else:
        rand = random.randint(0,len(path)-1)
        path.remove(path[rand])
        new_G = remove_paths(new_G, path, delete_entry_node, delete_sink_node)
    return new_G
    

def solve_bubble(graph, ancestor_node, descendant_node):
    """Take graph, ancestor node, descendant node and return clean graph"""
    new_G = graph
    return new_G




graph_1 = nx.DiGraph()
graph_1.add_weighted_edges_from([(1, 2, 10), (3, 2, 10), (2, 4, 15), 
                                     (4, 5, 15), (2, 10,10), (10, 5,10),
                                     (2, 8, 3), (8, 9, 3), (9, 5, 3),
                                     (5, 6, 10), (5, 7, 10)])
graph_1 = solve_bubble(graph_1, 2, 5)
assert (2,8) not in graph_1.edges()
assert (8,9) not in graph_1.edges()
assert (9,5) not in graph_1.edges()
assert (2,10) not in graph_1.edges()
assert (10, 5) not in graph_1.edges()
assert (2,4) in graph_1.edges()
assert (4,5) in graph_1.edges()
assert 8 not in graph_1.nodes()
assert 9 not in graph_1.nodes()
assert 10 not in graph_1.nodes()
assert 2 in graph_1.nodes()
assert 5 in graph_1.nodes()
graph_2 = nx.DiGraph()
graph_2.add_weighted_edges_from([(1, 2, 10), (3, 2, 10), (2, 4, 10), 
                                     (4, 5, 10), (2, 10,10), (10, 5,10),
                                     (2, 8, 10), (8, 9, 10), (9, 5, 10),
                                     (5, 6, 10), (5, 7, 10)])
graph_2 = solve_bubble(graph_2, 2, 5)
assert (2,4) not in graph_2.edges()
assert (4,5) not in graph_2.edges()
assert (2,10) not in graph_1.edges()
assert (10, 5) not in graph_1.edges()
assert (2,8) in graph_2.edges()
assert (8,9) in graph_2.edges()
assert (9,5) in graph_2.edges()


def simplify_bubbles(graph):
    """Take graph and return it without bubble"""


def solve_entry_tips(graph, entry_node):
    """Take graph and entry nodes and return graph whithout entry indesirable path"""


def solve_out_tips(graph, sink_node):
    """Take graph and sink nodes and return graph whithout entry indesirable path"""


def main():
    """Main function"""
    parser = argparse.ArgumentParser(description='Graphe de de Bruijn')
    parser.add_argument('-i', metavar='FASTQ', type=str, help='File FASTQ', required=True)
    parser.add_argument('-k', metavar='Kmer', type=int, help='Kmer size', default='21')
    #parser.add_argument('-o', metavar='Confiq', type=str, help='Config file', required=True)
    args = parser.parse_args()
    fastq_file = args.i
    kmer_size = args.k
    #config_file = args.o
    
    graph = build_graph(build_kmer_dict(fastq_file, kmer_size))
    tupple = get_contigs(graph, get_starting_nodes(graph), get_sink_nodes(graph))
    save_contigs(tupple, 'eva71_hundred_reads.fq.out')


if __name__ == '__main__':
    main()
