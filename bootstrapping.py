def generate_bootstrap_genes( sequences ):
    new_sequences = []

    return new_sequences

def progressive_alignment( sequences, guide_tree):
    prog_align = []

    if(type(guide_tree[0]) == tuple):
        seq_1 = progressive_alignment( sequences, guide_tree[0] )
    else:
        seq_1 = sequences[guide_tree[0]]

    if(type(guide_tree[1]) == tuple):
        seq_2 = progressive_alignment( sequences, guide_tree[1] )
    else:
        seq_2 = sequences[guide_tree[1]]

    prog_align = align_alignments(seq_1, seq_2)
    
    return prog_align

def align_alignments( sequence_1, sequence_2 ):

    return

def calculate_confidences( clade_count_dict, bootstrap_num ):
    confidences = {}

    for i, j in clade_count_dict.items():
        confidences[ i ] = j / bootstrap_num

    return confidences