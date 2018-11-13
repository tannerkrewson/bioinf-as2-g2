def calculate_upgma( distance_matrix ):
    seq_map = []
    for i in range(0, len(distance_matrix)):
        seq_map.append(i)

    new_matrix = distance_matrix.copy()
    
    while len(new_matrix) > 2:
        print_mat(new_matrix)
        print(seq_map)
        min_index = find_lowest_distance( new_matrix )

        X = min_index[0]
        Y = min_index[1]

        print(min_index)

        seq_map[X] = ( seq_map[X], seq_map[Y] )
        seq_map.pop(Y)

        for i in range(0, len(new_matrix)):
            # average the two values, and store in the left/top most
            # one's original position
            new_matrix[i][X] += new_matrix[i][Y]
            new_matrix[i][X] /= 2

            # matrix is a square so i'll do the other direction
            # in the same loop without consequence
            new_matrix[X][i] += new_matrix[Y][i]
            new_matrix[X][i] /= 2

        # remove the non-averaged row
        new_matrix.pop(Y)

        # remove the non-averaged col
        for i in range(0, len(new_matrix)):
            new_matrix[i].pop(Y)
    
    return ( seq_map[0], seq_map[1] )
    

def refactor_matrix( distance_matrix ):
    # for fixing the matrix after creating a clade
    return False

def find_lowest_distance( distance_matrix ):
    min_val = distance_matrix[0][1]
    min_index = [ 0, 1 ]

    for i in range(0, len(distance_matrix)):
        for j in range(i+1, len(distance_matrix[i])):
            if distance_matrix[i][j] < min_val:
                min_index[0] = i
                min_index[1] = j
                min_val = distance_matrix[i][j]

    return min_index

def generate_matrix( sequences ):
    return False

def test():
    test_dist_matrix = [ 
        [0, 13, 11,  7, 12, 16, 15],
        [0,  0,  2, 11, 14, 13,  5],
        [0,  0,  0,  9, 18, 15,  3],
        [0,  0,  0,  0,  8, 14, 13],
        [0,  0,  0,  0,  0, 18, 13],
        [0,  0,  0,  0,  0,  0, 14],
        [0,  0,  0,  0,  0,  0,  0]
    ]

    calculate_upgma(test_dist_matrix)

def print_mat(mat):
    for row in mat:
        print(row)
    print()

test()