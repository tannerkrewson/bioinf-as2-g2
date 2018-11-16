from readfasta import readfasta
import glob, os
import numpy
import multiprocessing

from tree_analysis import build_clade_count_dict, clade_search
from bootstrapping import generate_bootstrap_genes, calculate_confidences
from seq_distance import find_distance
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

    print( original_tree )

    '''
    clade_count_dict = {}

    # count the clades of the original tree and add them as keys
    build_clade_count_dict( original_tree, clade_count_dict )

    BOOTSTRAP_TIMES = 10

    # count each clade in bootstrapped trees that match a clade from the original tree
    for i in range (0, BOOTSTRAP_TIMES):
        bootstrapped_genes = generate_bootstrap_genes( genes )

        this_tree = generate_tree( bootstrapped_genes )

        clade_search( this_tree, clade_count_dict )

    # return a dictionary containing the clades as keys mapped to their confidence
    clade_confidences = calculate_confidences( clade_count_dict, BOOTSTRAP_TIMES )

    # TODO: add function to append clade confidence values to the original tree
    '''


def generate_tree( genes ):
    # find the distance between each gene
    distance_matrix = numpy.zeros((len(genes), len(genes)), dtype=float)
    read_precalculated_distances(distance_matrix)

    alignments_matrix = []
    for i in range( 0, len( genes ) ):
        new_row = []
        for j in range( 0, len( genes ) ):
            new_row.append([])
        alignments_matrix.append(new_row)

    def store_distance(result):
        distance = result[0]
        i = result[1]
        j = result[2]

        distance_matrix[i, j] = distance
        alignments_matrix[i][j].append(result[3][0])
        alignments_matrix[i][j].append(result[3][1])

    pool = multiprocessing.Pool()

    for i in range( 0, len( genes ) ):
        for j in range( i+1, len( genes ) ):
            if distance_matrix[i, j] == 0:
                pool.apply_async(find_distance, args = (genes[i], genes[j], i, j), callback = store_distance)
    
    pool.close()
    pool.join()

    file = open("distances.csv", "w")
    for i in range( 0, len( genes ) ):
        for j in range( i+1, len( genes ) ):
          file.write(str(i) + "," + str(j) + "," + str(distance_matrix[i, j]) + "\n")
    file.close()

    for i in range( 0, len( genes ) ):
        for j in range( i+1, len( genes ) ):
            file = open(str(i) + "," + str(j) + ".txt", "w")
            file.write(genes[i][0] + "\n")
            file.write(alignments_matrix[i][j][0] + "\n")
            file.write(genes[j][0] + "\n")
            file.write(alignments_matrix[i][j][1] + "\n")
            file.close()

    # use upgma to generate a tree from the distance matrix
    return calculate_upgma( distance_matrix )


def read_precalculated_distances( distance_matrix ):
    if not os.path.isfile("distance_database.csv"):
        return
    file = open("distance_database.csv", "r") 
    for line in file: 
        things = line.split(",")
        i = int(things[0])
        j = int(things[1])
        if len(things) == 5:
            distance_matrix[i, j] = float(things[4])

if __name__ == '__main__':
    main()
