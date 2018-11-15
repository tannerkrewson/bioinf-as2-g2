from readfasta import readfasta
from random import randint
from scipy import stats
import glob, os

'''
main: read in the fasta files, run our algorithm to guess the reading
    frame, and statistically compare it to random
'''
def main():
    print( "*****\nBioinformatics - Assignment 2 - Group 2\n*****\n" )

    # scan in all fasta files in the "genes" directory
    os.chdir( os.getcwd() + "/genes/" )
    for file in glob.glob( "*.fasta" ):
        print( file )

        genes = readfasta( file )
        for gene in genes:
            print( gene[1][:60] )

main()
