def find_distance( gene1, gene2 ):

    aligned_genes = align_gene( gene1, gene2 )
    
    return dK2P( aligned_genes[0], aligned_genes[1] )