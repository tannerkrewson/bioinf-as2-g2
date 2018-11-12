def generate_bootstrap_genes( sequences ):
    new_sequences = []

    align_all_sequences( sequences )

    sequence_length = 1 # TODO: add correct value

    for i in range(0, sequence_length):
        rand_pos = random.randint(0, sequence_length)
        for j in range(0, len(sequences)):
            new_sequences[j] += sequences[j][randpos]

    return new_sequences

def align_all_sequences( sequences ):
    # TODO: align all sequences so they can be picked from for bootstrapping
    return False

def calculate_confidences( clade_count_dict, bootstrap_num ):
    confidences = {}

    for i, j in clade_count_dict.items():
        confidences[i] = j / bootstrap_num

    return confidences