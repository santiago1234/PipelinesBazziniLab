from Bio.Seq import Seq


def get_tripletes(seq_d):
    """
    Iterates over the triples of seq_d on the first
    frame
    args:
        seq_d: str representing a dna sequence
    returns:
        an iterator of the tripletes
    """
    return (seq_d[i: i + 3] for i in range(0, len(seq_d) -2, 3))

def find_premature_stop(seq_d):
    """
    Finds the postion of premature stop codon
    if it exist
    args:
        seq_d: str sequence
    returns:
        position of premature if exist
        else -
    """
    def translate_triplet(tri):
        return str(Seq(tri).translate())
    counter = 0
    for triplet in get_tripletes(seq_d):
        if translate_triplet(triplet) == '*':
            return counter
        counter += 1
    return -1

