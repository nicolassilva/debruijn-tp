"""Find Kmer in sequence"""
import os
import argparse
import statistics as st
import random
import networkx as nx


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
    for node in graph.nodes():
        if len(list(graph.predecessors(node))) == 0:
            entres.append(node)
    return entres


def get_sink_nodes(graph):
    """Prend un graphe et retourne une liste de noeuds de sortis"""
    sortis = []
    for node in graph.nodes():
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
        file.write('>contig_' + str(i) + ' len=' + str(tupple[i][1]) + '\n' 
                   + str(fill(tupple[i][0])) + '\n')
    file.close()


def std(list_val):
    """Take list of values and return standard deviation"""
    return st.stdev(list_val)


def path_average_weight(graph, path):
    """Take a graph and a path and return average weigth"""
    new_g = graph.subgraph(path)
    wei = []
    for arretes in new_g.edges(data=True):
        wei.append(arretes[2]['weight'])
    mean_wei = st.mean(wei)
    return mean_wei


def remove_paths(graph, path, delete_entry_node, delete_sink_node):
    """Take graph and path and remove nodes
    delete_entry_node pour indiquer si les noeuds d’entrée seront supprimés
    et delete_sink_node pour indiquer si les noeuds de sortie seront supprimés
    Return clean graph"""
    new_g = graph
    for i in range(len(path)):
        new_g.remove_nodes_from(path[i][1:-1]) #On retire les noeuds entre les entrés et les sortis
        if delete_entry_node == True:
            new_g.remove_node(path[i][0])
        if delete_sink_node == True:
            new_g.remove_node(path[i][-1])
    return new_g


def select_best_path(graph, path, path_len, path_wei,
                     delete_entry_node=False, delete_sink_node=False):
    """Take graph, path, path length, path weigth and delete entry or sink node
    Function return clean graph
    By default delete entry or sink node is set as False"""
    new_g = graph
    random.seed(9001)
    best_wei_indice = []
	#Boucle permet de regarder l'index dans la liste du chemin le plus lourd
    for i, exp in enumerate(path_wei):
        if exp == max(path_wei):
            best_wei_indice.append(i)
    best_len_indice = []
	#Boucle permet de regarder l'index dans la liste du chemin le plus long
    for i, exp in enumerate(path_len):
        if exp == max(path_len):
            best_len_indice.append(i)
	#Boucle à 3 conditions:
	#- Si un chemin est plus lourd alors il est choisi et les autres sont retirés
    if len(best_wei_indice) == 1:
        path.remove(path[best_wei_indice[0]])
        new_g = remove_paths(new_g, path, delete_entry_node, delete_sink_node)
	#- Si les poids sont identiques, on regarde les longueurs des chemins
	# Si un chemin est plus long alors il est choisi et les autres sont retirés
    elif len(best_len_indice) == 1:
        path.remove(path[best_len_indice[0]])
        new_g = remove_paths(new_g, path, delete_entry_node, delete_sink_node)
	#Si aucun chemin n'est plus lourd ou plus long, alors on sélectionne un chemin au hasard
	#Ce chemin est gardé et les autres retirés du graphe
    else:
        rand = random.randint(0, len(path)-1)
        path.remove(path[rand])
        new_g = remove_paths(new_g, path, delete_entry_node, delete_sink_node)
    return new_g


def solve_bubble(graph, ancestor_node, descendant_node):
    """Take graph, ancestor node, descendant node and return clean graph"""
    new_g = graph
	#On regarde tout les chemins possible entre l'ancestor node et le descendant node
    all_path = nx.all_simple_paths(graph, ancestor_node, descendant_node)
    all_path = list(all_path)
    weigth = []
	#On calcul les poids des chemins trouvés
    for i in range(len(list(all_path))):
        weigth.append(path_average_weight(new_g, all_path[i]))
    length = []
	#On calcul les longueurs des chemins trouvés
    for i in range(len(list(all_path))):
        length.append(len(all_path[i]))
	#On selectionne le meilleur chemin
    new_g = select_best_path(new_g, all_path, length, weigth)
    return new_g


def simplify_bubbles(graph):
    """Take graph and return it without bubble"""
    new_g = graph
    bubble = []
	#Pour trouver les bulles on regarde les noeuds qui ont plus que 1 predecesseur
    for node in new_g.nodes():
        pred = list(graph.predecessors(node))
		#Si un noeud à plusieurs prédécesseur alors on regarde si ces prédécesseurs on un ancêtre commun
		#Si c'est le cas, alors on ajoute dans une liste le noeud de début et le noeud de fin de bulle
        if len(pred) > 1:
            anc = nx.lowest_common_ancestor(new_g, pred[0], pred[1])
            bubble.append([anc, node])
	#On utilise la fonction solve_bubble pour éliminer les bulles en envoyant dans la fonction
	#les début et fin de bulles
    for i in range(len(bubble)):
        new_g = solve_bubble(new_g, bubble[i][0], bubble[i][1])
    return new_g


def solve_entry_tips(graph, entry_node):
    """Take graph and entry nodes and return graph whithout entry indesirable path"""
    new_g = graph
    des = []
	#On regarde les prédecesseurs de chaque noeuds
	#Si plus que 1 prédécesseurs alors on ajoute le noeud dans une liste
    for node in new_g.nodes():
        pred = list(new_g.predecessors(node))
        if len(pred) > 1:
            des.append(node)
	#Si la liste de noeud est à 0 alors il n'y a pas plusieurs noeuds d'entrés
	#On return le graph sans modification
    if len(des) == 0:
        return new_g
    path = []
    weight = []
    length = []
	#On calcule les différents chemin d'entrés possible, les poids et longueurs des chemins
    for i in range(len(entry_node)):
        path.append(list(nx.all_simple_paths(new_g, entry_node[i], des[0])))
        path[i] = path[i][0]
        weight.append(path_average_weight(new_g, path[i]))
        length.append(len(path[i]))
	#On sélectionne le meilleur chemin
    new_g = select_best_path(new_g, path, length, weight, True, False)
    return new_g


def solve_out_tips(graph, sink_node):
    """Take graph and sink nodes and return graph whithout entry indesirable path"""
    new_g = graph
    anc = []
	#On récupere les noeuds ancêtres des noeuds de sortis et on les stoke dans une liste
    for i in range(1, len(sink_node)):
        anc.append(nx.lowest_common_ancestor(new_g, sink_node[i-1], sink_node[i]))
    path = []
    weight = []
    length = []
	#On calcule les différents chemin d'entrés possible, les poids et longueurs des chemins
    for i in range(len(sink_node)):
        path.append(list(nx.all_simple_paths(new_g, anc[0], sink_node[i])))
        path[i] = path[i][0]
        weight.append(path_average_weight(new_g, path[i]))
        length.append(len(path[i]))
	#On sélectionne le meilleur chemin
    new_g = select_best_path(new_g, path, length, weight, False, True)
    return new_g


def main():
    """Main function"""
	#On demande les arguments
    parser = argparse.ArgumentParser(description='Graphe de de Bruijn')
    parser.add_argument('-i', metavar='FASTQ', type=str, help='File FASTQ', required=True)
    parser.add_argument('-k', metavar='Kmer', type=int, help='Kmer size', default='21')
    parser.add_argument('-o', metavar='Confiq', type=str, help='Contig file', required=True)
    args = parser.parse_args()
	#On import les arguments dans des variables
    fastq_file = args.i
    kmer_size = args.k
    contig_file = args.o
    #On construit le graphe
    graph = build_graph(build_kmer_dict(fastq_file, kmer_size))
    graph = simplify_bubbles(graph)
    graph = solve_entry_tips(graph, get_starting_nodes(graph))
    graph = solve_out_tips(graph, get_sink_nodes(graph))
	#On crée un contig contenant notre séquence et on le sauvagarde
    tupple = get_contigs(graph, get_starting_nodes(graph), get_sink_nodes(graph))
    save_contigs(tupple, contig_file)
	#Création index: makeblastdb -in eva71.fna -dbtype nucl
	#Blast verification: blastn -query eva71_plus_perfect.fq.out -db eva71.fna


if __name__ == '__main__':
    main()
