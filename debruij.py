import argparse

def main():
    """Main function"""
    parser = argparse.ArgumentParser(description='Graphe de de Bruijn')
    i = parser.add_argument('-i', metavar='FASTQ', type=str, help='File FASTQ')
    k = parser.add_argument('-k', metavar='Kmer', type=int, help='Kmer size', default='21')
    o = parser.add_argument('-o', metavar='Confiq', type=str, help='Config file')
    args = parser.parse_args()

if __name__ == '__main__':
    #fast_file = i
    kmer_size = k
    #config_file = o
    print(k)
    main()
