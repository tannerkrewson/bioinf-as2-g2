import copy

def calculate_upgma( distance_matrix ):
    seq_tree = []
    seq_map = []
    for i in range(0, len(distance_matrix)):
        seq_tree.append(i)
        seq_map.append([i])

    new_matrix = copy.deepcopy(distance_matrix)
    
    while len(new_matrix) > 2:

        min_index = find_lowest_distance( new_matrix )

        X = min_index[0]
        Y = min_index[1]

        seq_tree[X] = ( seq_tree[X], seq_tree[Y] )
        seq_tree.pop(Y)

        seq_map[X] += seq_map[Y]
        seq_map.pop(Y)

        # remove the non-averaged row
        new_matrix.pop(Y)

        # remove the non-averaged col
        for i in range(0, len(new_matrix)):
            new_matrix[i].pop(Y)

        for i in range(0, X):
            # average the two values, and store in the left/top most
            # one's original position
            sum = 0
            count = 0
            for m in seq_map[X]:
                for n in seq_map[i]:
                    a = min(m, n)
                    b = max(m, n)
                    sum += distance_matrix[a][b]
                    count += 1
            new_matrix[i][X] = sum / count

        for i in range(X+1, len(new_matrix)):
            sum = 0
            count = 0
            for m in seq_map[X]:
                for n in seq_map[i]:
                    a = min(m, n)
                    b = max(m, n)
                    sum += distance_matrix[a][b]
                    count += 1
            new_matrix[X][i] = sum / count

    return ( seq_tree[0], seq_tree[1] )
    

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