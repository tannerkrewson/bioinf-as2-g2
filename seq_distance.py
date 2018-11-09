import math

def find_distance( gene1, gene2 ):

    aligned_genes = align_gene( gene1, gene2 )
    
    return dK2P( aligned_genes[0], aligned_genes[1] )

#x and y are the two sequences being compared 
def dK2P( x, y): 
    count_transition = 0 #the total differ by transition
    count_transversion = 0 #the total that differ by transversion

    transitions = [('G','A'), ('A', 'G'), ('C', 'T'), ('T', 'C')]
    transversions = [('G','T'), ('T', 'G'), ('A', 'T'), ('T', 'A'),('G', 'C'), ('C', 'G'), ('A', 'C'), ('C', 'A')]

    for i in range(len(x)):
        if x[i] != y[i]:
            site = (x[i], y[i])
            for j in range(len(transitions)):
                if site == transitions[j]:
                    count_transition = count_transition + 1
            for j in range(len(transversions)):
                if site == transversions[j]:
                    count_transversion = count_transversion + 1


    S = count_transition / len(x)
    V = count_transversion / len(x)


    distance = -0.5 * math.log(1 - 2*S - V) - 0.25 * math.log(1 - 2*V) #K2P Formula 

    return distance
