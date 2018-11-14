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

def align_gene(sequence1, sequence2):
    gp = -2
    mb = 1
    mp = -1
    
    seq1 = "-" + sequence1 #seq1 is the sequence on the horizontal
    seq2 = "-" + sequence2 #seq2 is on the vertical

    A = [] #dynamic programing matrix wihtout letters

    #initialize first row
    row = []
    for j in range(len(seq1)): 
        row.append(j * gp)
    A.append(row)

    #initialize first column
    for i in range(1,len(seq2)): 
        A.append([i * gp])

    for i in range(1,len(seq2)):
        currentRow = A[i]
        for j in range(1,len(seq1)):
            if (seq2[i] == seq1[j]):
                x = A[i-1][j-1] + mb #value from diagonal if they match
            else:
                x = A[i-1][j-1] + mp #value from diagonal if they don't match
            y = A[i][j-1] + gp #value from the left
            z = A[i-1][j] + gp #value from above
            currentRow.append(max(x,y,z))
            
        A[i] = currentRow

    #figure out how to trace back
    i = len(seq2)-1
    j = len(seq1)-1
    
    traceBackDirections = ""
    while i > 0 or j > 0:
        if (seq2[i] == seq1[j]):
            x = A[i-1][j-1] + mb #value from diagonal if they match
        else:
            x = A[i-1][j-1] + mp #value from diagonal if they don't match
        y = A[i][j-1] + gp #value from the left
        z = A[i-1][j] + gp #value from above

        if A[i][j] == x and i >= 0 and j >= 0:
            #trace back to the diagonal
            traceBackDirections = traceBackDirections + "d" 
            i = i - 1
            j = j - 1
        elif A[i][j] == y and i >= 0 and j >= 0:
            #trace back to the left
            traceBackDirections = traceBackDirections + "b"
            j = j - 1
        else:
            #trace back to above
            traceBackDirections = traceBackDirections + "u"
            i = i - 1

    #align with the trace back directions 
    seq1spot = len(sequence1) - 1
    seq2spot = len(sequence2) - 1
    newSeq1 = ""
    newSeq2 = ""
    for i in range(len(traceBackDirections)):
        if traceBackDirections[i] == "d":
            newSeq1 = sequence1[seq1spot] + newSeq1
            newSeq2 = sequence2[seq2spot] + newSeq2
            seq1spot = seq1spot - 1
            seq2spot = seq2spot - 1
        if traceBackDirections[i] == "u":
            newSeq1 = '-' + newSeq1
            newSeq2 = sequence2[seq2spot] + newSeq2
            seq2spot = seq2spot - 1
        if traceBackDirections[i] == "b":
            newSeq2 = '-' + newSeq2
            newSeq1 = sequence2[seq2spot] + newSeq2
            seq2spot = seq2spot - 1

    return [newSeq1, newSeq2]
