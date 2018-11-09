from readfasta import readfasta
import glob, os

from tree_analysis import build_clade_count_dict, count_clades
from bootstrapping import get_bootstrapped_genes
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

        original_tree = generate_tree(genes)

        # count the clades of the original tree and add them as keys
        clade_count_dict = build_clade_count_dict( original_tree )

        BOOTSTRAP_TIMES = 10

        # count each clade in bootstrapped trees that match a clade from the original tree
        for i in range (0, BOOTSTRAP_TIMES):
            bootstrapped_genes = get_bootstrapped_genes( genes )

            this_tree = generate_tree( bootstrapped_genes )

            count_clades( this_tree, clade_count_dict )

        # find the confidence of each clade
        for clade in clade_count_dict:
            the_count = clade_count_dict[clade]

            clade_count_dict[clade] = the_count / BOOTSTRAP_TIMES

        # at this point, the clade_count_dict now maps clades to their respective confidence level
        print( clade_count_dict )


def generate_tree(genes):
    # find the distance between each gene
    distance_matrix = []
    for i in range(0, len(genes)):

        # make it 2d by adding new list everytime
        new_row = []

        for j in range(i, len(genes)):
            new_row.append( find_distance(genes[i], genes[j]) )
        
        distance_matrix.append( new_row )

    # use upgma to generate a tree from the distance matrix
    return calculate_upgma( distance_matrix )


main()
