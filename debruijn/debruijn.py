"""Find Kmer in sequence"""
import argparse
import networkx as nx
import os

def std():
    pass

def path_average_weight():
    pass


def remove_paths():
    pass


def select_best_path():
    pass


def solve_bubble():
    pass


def simplify_bubbles():
    pass


def solve_entry_tips():
    pass


def solve_out_tips():
    pass


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
    for kmer,val in kmer_dict.items():
        graph.add_edge(kmer[:-1],kmer[1:], weight=val)
    return graph


def get_starting_nodes(graph):
    """Prend un graphe et retourne une liste de noeuds d'entrés"""
    entres = []
    for node in graph.nodes:
        if len(list(graph.predecessors(node)))==False:
            entres.append(node)
    return entres


def get_sink_nodes(graph):
    """Prend un graphe et retourne une liste de noeuds de sortis"""
    sortis = []
    for node in graph.nodes:
        if len(list(graph.successors(node)))==False:
            sortis.append(node)
    return sortis


def get_contigs(graph,entres,sortis):
    """Retourne liste de tuple(contig, taille du contig)"""
    all_paths = []
    for node in entres:
        for node2 in sortis:
            path = list(nx.all_simple_paths(graph,node,node2))
            if len(path) > 0:
                all_paths.append(path)
    tupple = []
    for i in range(len(all_paths)):
        for k in range(len(all_paths[i])):
            contig = str(all_paths[i][k][0])
            for j in range(1,len(all_paths[i][k])):
                string = str(all_paths[i][k][j][-1])
                contig += string[-1]
            tupple.append([contig,len(contig)])
    return tupple


def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))


def save_contigs(tupple,file_name):
    """Sort un fichier avec la liste tupple"""
    file = open("../data/"+file_name,'w+')
    for i in range(len(tupple)):
        file.write('>Contig numero ' + str(i) + ' len=' + str(tupple[i][1]) + '\n' + str(fill(tupple[i][0])) + '\n')


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
    save_contigs(tupple,'eva71_hundred_reads.fq.out')
    

if __name__ == '__main__':
    main()
