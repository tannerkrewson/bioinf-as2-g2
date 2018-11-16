from readfasta import readfasta
import glob, os
import numpy
import multiprocessing

from tree_analysis import build_clade_count_dict, clade_search
from bootstrapping import generate_bootstrap_genes, calculate_confidences, \
progressive_alignment, reorder_alignments
from seq_distance import find_distance, dK2P
from upgma import calculate_upgma

'''
main: read in the fasta files, and do things TODO
'''
def main():
    print( "*****\nBioinformatics - Assignment 2 - Group 2\n*****\n" )

    # scan in all fasta files in the "genes" directory
    os.chdir( os.getcwd() + "/genes/" )
    file = glob.glob( "*.fasta" )[0]

    # read all the genes from the fasta file
    print( file )
    genes = readfasta( file )

    original_tree = generate_tree( genes )

    print( original_tree[0] )

    clade_count_dict = {}

    # count the clades of the original tree and add them as keys
    build_clade_count_dict( original_tree[0], clade_count_dict )

    BOOTSTRAP_TIMES = 20
    
    multi_aligned_sequences = \
    progressive_alignment(genes, original_tree[0], original_tree[1])

    multi_aligned_sequences = reorder_alignments(multi_aligned_sequences)

    for sequence in multi_aligned_sequences:
        print(sequence[0][:120])

    #count each clade in bootstrap trees matching a clade from original tree
    for i in range (0, BOOTSTRAP_TIMES):
        bootstrapped_genes = generate_bootstrap_genes(multi_aligned_sequences)

        this_tree = generate_boots_tree( bootstrapped_genes )

        print("Bootstrap Tree ", i)
        print(this_tree)

        clade_search( this_tree, clade_count_dict )

    # return a dict containing the clades as keys mapped to their confidence
    clade_confidences = \
    calculate_confidences( clade_count_dict, BOOTSTRAP_TIMES )

    for clade, confidence in clade_confidences.items():
        print("clade: ", clade)
        print("confidence: ", confidence)
    

def generate_tree( genes ):
    # find the distance between each gene
    distance_matrix = numpy.zeros((len(genes), len(genes)), dtype=float)

    # if distances already calulated and stored in file, read them in
    read_precalculated_distances(distance_matrix)

    # fill the alignments matrix
    alignments_matrix = []
    for i in range( 0, len( genes ) ):
        new_row = []
        for j in range( 0, len( genes ) ):
            new_row.append([])
        alignments_matrix.append(new_row)

    # if alignments already calulated and stored in file, read them in
    # read_precalculated_alignments(alignments_matrix)

    def store_distance(result):
        distance = result[0]
        i = result[1]
        j = result[2]
        #aligned_genes = result[3]

        distance_matrix[i, j] = distance
        #alignments_matrix[i][j].append(aligned_genes[0])
        #alignments_matrix[i][j].append(aligned_genes[1])
        print(i, j)

    pool = multiprocessing.Pool()

    for i in range( 0, len( genes ) ):
        for j in range( i+1, len( genes ) ):
            # if the distance and alignment don't exist, recalculate them

            if distance_matrix[i, j] == 0:
                # parallel
                pool.apply_async(find_distance, args = \
                    (genes[i], genes[j], i, j), callback = store_distance)
                
                # non-parallel
                # store_distance(find_distance(genes[i], genes[j], i, j))
    
    pool.close()
    pool.join()

    # store the calculated distances in a file
    file = open("distances_database_new.csv", "w")
    for i in range( 0, len( genes ) ):
        for j in range( i+1, len( genes ) ):
            file.write(str(i) + "," + str(j) + "," + \
                str(distance_matrix[i, j]) + "\n")
    file.close()

    #exit()

    # store each alignment in a file
    if not os.path.isfile("alignments_database.csv"):
        file = open("alignments_database.csv", "w")
        for i in range( 0, len( genes ) ):
            for j in range( i+1, len( genes ) ):
                if len(alignments_matrix[i][j]) == 2:
                    file.write(str(i) + "," + str(j) + "," + \
                        alignments_matrix[i][j][0] + "," + \
                        alignments_matrix[i][j][1] + "\n")
        file.close()

    # use upgma to generate a tree from the distance matrix
    return [calculate_upgma( distance_matrix ), alignments_matrix]

def generate_boots_tree( genes ):
    # find the distance between each gene
    distance_matrix = numpy.zeros((len(genes), len(genes)), dtype=float)

    def store_distance(result):
        distance = result[0]
        i = result[1]
        j = result[2]

        distance_matrix[i, j] = distance

    pool = multiprocessing.Pool()

    for i in range( 0, len( genes ) ):
        for j in range( i+1, len( genes ) ):
            pool.apply_async(find_boots_distance, args = \
                (genes[i], genes[j], i, j), callback = store_distance)
    
    pool.close()
    pool.join()

    # use upgma to generate a tree from the distance matrix
    return calculate_upgma( distance_matrix )

def find_boots_distance( gene1, gene2, i, j):
    
    distance = dK2P( gene1, gene2 )

    return [ distance, i, j]

def read_precalculated_distances( distance_matrix ):
    if not os.path.isfile("distances_database.csv"):
        return
    file = open("distances_database.csv", "r") 
    for line in file: 
        things = line.split(",")
        i = int(things[0])
        j = int(things[1])
        if len(things) == 3:
            distance_matrix[i, j] = float(things[2])

def read_precalculated_alignments( alignments_matrix ):
    if not os.path.isfile("alignments_database.csv"):
        return
    file = open("alignments_database.csv", "r") 
    for line in file: 
        things = line.split(",")
        i = int(things[0])
        j = int(things[1])
        if len(things) == 4:
            print("found in file", i , j)
            alignments_matrix[i][j].append(things[2])
            alignments_matrix[i][j].append(things[3])

if __name__ == '__main__':
    main()
