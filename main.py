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
    for file in glob.glob( "*.fasta" ):

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
    distance_matrix = numpy.zeros((len(genes), len(genes)), dtype=int)

    pool = multiprocessing.Pool()

    for i in range( 0, len( genes ) ):
        for j in range( i+1, len( genes ) ):
            pool.apply_async(find_distance, args = (genes[i][1], genes[j][1], i, j, distance_matrix))
    
    pool.close()
    pool.join()

    # use upgma to generate a tree from the distance matrix
    return calculate_upgma( distance_matrix )


if __name__ == '__main__':
    main()
