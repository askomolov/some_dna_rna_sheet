DNA = ('A', 'T', 'G', 'C', 'a', 't', 'g', 'c')
RNA = ('A', 'U', 'G', 'C', 'a', 'u', 'g', 'c')


def input_check(seqs):
    if len(seqs) == 0:
        print('No sequence in input')
        return 0
    for i in seqs:
        if is_dna(i) is False and is_rna(i) is False:
            print('One sequence is not DNA or RNA')
            return 0
    return 1


def transcribe(seq):
    new_seq = seq.replace('T', 'U')
    new_seq = new_seq.replace('t', 'u')
    return new_seq


def reverse(seq):
    return seq[::-1]


def complement(seq):
    new_seq = ''
    for i in seq:
        if i == 'A':
            new_seq = new_seq + 'T'
        elif i == 'a':
            new_seq = new_seq + 't'
        elif i == 'T':
            new_seq = new_seq + 'A'
        elif i == 't':
            new_seq = new_seq + 'a'
        elif i == 'C':
            new_seq = new_seq + 'G'
        elif i == 'c':
            new_seq = new_seq + 'g'
        elif i == 'G':
            new_seq = new_seq + 'C'
        elif i == 'g':
            new_seq = new_seq + 'c'
        elif i == 'U':
            new_seq = new_seq + 'A'
        elif i == 'u':
            new_seq = new_seq + 'a'
    return new_seq


def reverse_complement(seq):
    new_seq = complement(seq)
    new_seq = reverse(new_seq)
    return new_seq


def is_dna(seq):
    return set(seq).issubset(DNA)


def is_rna(seq):
    return set(seq).issubset(RNA)
