import random

def generate_bootstrap_genes( sequences, original_tree, alignments ):
    new_sequences = []

    aligned_sequences = progressive_alignment(sequences, original_tree, alignments)

    for i in range(len(aligned_sequences[0])):
        random_index = random.randint(0, len(aligned_sequences[0]) - 1)
        for j in range(len(aligned_sequences)):
            new_sequences[j] += aligned_sequences[j][random_index]

    return new_sequences

def progressive_alignment( sequences, guide_tree, alignments ):

    if type(guide_tree[0]) == tuple:
        seq_1 = progressive_alignment( sequences, guide_tree[0] )
    else:
        leaf_aligns = alignments[guide_tree[0]][guide_tree[1]]
        return leaf_aligns
        # seq_1 = sequences[guide_tree[0]]

    if type(guide_tree[1]) == tuple:
        seq_2 = progressive_alignment( sequences, guide_tree[1] )
    else:
        leaf_aligns = alignments[guide_tree[0]][guide_tree[1]]
        return leaf_aligns
        # seq_2 = sequences[guide_tree[1]]
    
    return align_alignments(seq_1, seq_2)

def align_alignments( alignment_1, alignment_2 ):
    gp = -2.5
    mb = 1
    mp = -2

    A = numpy.zeros((len(alignment_2), len(alignment_1)), dtype=int)

    #initialize first row
    for j in range(0, len(alignment_1[0])): 
        A[0, j] = j * gp

    #initialize first column
    for i in range(1, len(alignment_2[0])): 
        A[i, 0] = i * gp

    for i in range(1, len(alignment_2[0])):
        if (i % 100) == 0:
            print(i)
        for j in range(1, len(alignment_1[0])):
            x = A[i-1, j-1]
            x += calculate_cell_addition( alignment_1, alignment_2, i, j)

            y = A[i, j-1] + (gp * len(alignment_1)) #value from the left
            z = A[i-1, j] + (gp * len(alignment_2)) #value from above

            A[i, j] = max(x, y, z)

    #figure out how to trace back
    i = len(alignment_2[0])-1
    j = len(alignment_1[0])-1
    
    traceBackDirections = ""
    while i > 0 or j > 0:
        x = A[i-1, j-1]
        x += calculate_cell_addition( alignment_1, alignment_2, i, j)
        y = A[i, j-1] + gp #value from the left
        z = A[i-1, j] + gp #value from above

        if A[i, j] == x and i >= 0 and j >= 0:
            #trace back to the diagonal
            traceBackDirections = traceBackDirections + "d" 
            i = i - 1
            j = j - 1
        elif A[i, j] == y and i >= 0 and j >= 0:
            #trace back to the left
            traceBackDirections = traceBackDirections + "b"
            j = j - 1
        else:
            #trace back to above
            traceBackDirections = traceBackDirections + "u"
            i = i - 1

    #align with the trace back directions 
    seq1spot = len(alignment_1[0]) - 1
    seq2spot = len(alignment_2[0]) - 1
    new_alignment_1 = ["" for i in len(alignment_1)]
    new_alignment_2 = ["" for i in len(alignment_2)]
    for i in range(len(traceBackDirections)):
        if traceBackDirections[i] == "d":
            for j in range(len(alignment_1)):
                new_alignment_1[j] = alignment_1[j][seq1spot] + new_alignment_1[j]
            for j in range(len(alignment_2)):
                new_alignment_2[j] = alignment_2[j][seq2spot] + new_alignment_2[j]
            seq1spot = seq1spot - 1
            seq2spot = seq2spot - 1
        if traceBackDirections[i] == "u":
            for j in range(len(alignment_1)):
                new_alignment_1[j] = '-' + new_alignment_1[j]
            for j in range(len(alignment_2)):
                new_alignment_2[j] = alignment_2[j][seq2spot] + new_alignment_2[j]
            seq2spot = seq2spot - 1
        if traceBackDirections[i] == "b":
            for j in range(len(alignment_1)):
                new_alignment_1[j] = alignment_1[j][seq1spot] + new_alignment_1[j]
            for j in range(len(alignment_2)):
                new_alignment_2[j] = '-' + new_alignment_2[j]
            seq1spot = seq1spot - 1

    new_alignments = new_alignment_1 + new_alignment_2
    return new_alignments

def calculate_cell_addition( alignment_1, alignment_2, i, j ):
    mb = 1
    mp = -2
    cell_addition = 0

    for sequence_1 in alignment_1
        for sequence_2 in alignment_2
            if (sequence_2[i] == sequence_1[j]):
                x += mb #value from diagonal if they match
            else:
                x += mp #value from diagonal if they don't match

    return cell_addition

def calculate_confidences( clade_count_dict, bootstrap_num ):
    confidences = {}

    for i, j in clade_count_dict.items():
        confidences[ i ] = j / bootstrap_num

    return confidences