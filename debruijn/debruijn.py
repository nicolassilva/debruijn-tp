"""Find Kmer in sequence"""
import argparse
import networkx as nx

def get_starting_nodes():
    pass


def std():
    pass


def get_sink_nodes():
    pass


def path_average_weight():
    pass


def remove_paths():
    pass


def select_best_path():
    pass


def save_contigs():
    pass


def get_contigs():
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
    g = build_graph(build_kmer_dict(fastq_file, kmer_size))
    print(g)

if __name__ == '__main__':
    main()
