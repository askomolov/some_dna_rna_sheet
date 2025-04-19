import numpy as np
from Bio import SeqIO
from Bio.SeqUtils import GC
import statistics
import os
import argparse


class WrongSequenceTypeError(TypeError):
    """
    For all classes.
    """
    pass


class ProteingEncodingError(TypeError):
    """
    For class Protein.
    """
    pass


class BiologicalSequence():
    dna = ('A', 'T', 'G', 'C')
    rna = ('A', 'U', 'G', 'C')
    aa_list = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

    def __init__(self, seq):
        self.seq = seq.upper()
        self.len_seq = len(seq)
        if self._is_seq() is False:
            raise WrongSequenceTypeError('This is not a biological sequence (DNA, RNA or protein)')

    def __len__(self):
        return self.len_seq

    def __str__(self):
        if self.len_seq > 75:
            seq_for_repr = self.seq[0:35] + '...' + self.seq[self.len_seq-35:self.len_seq]
        else:
            seq_for_repr = self.seq
        return f' {self.type} sequence: {seq_for_repr} with the len of {str(self.len_seq)}\n If you need full sequence refer to .seq attribute'

    def __repr__(self):
        return self.__str__()

    def __getitem__(self, slc):
        if isinstance(slc, int):
            return self.seq[slc]
        elif isinstance(slc, slice):
            return self.seq[slc]

    def _is_seq(self):
        if set(self.seq).issubset(self.dna):
            self.type = 'DNA'
            return True
        elif set(self.seq).issubset(self.rna):
            self.type = 'RNA'
            return True
        elif set(self.seq).issubset(self.aa_list):
            self.type = 'Protein'
            return True
        else:
            return False


class NucleicAcidSequence(BiologicalSequence):

    def __init__(self, seq):
        super().__init__(seq)
        if self.type == 'Protein':
            raise WrongSequenceTypeError('This is not a  nucleic acid sequence')

    def reverse(self):
        "Reverse nucleotide sequence"
        if self.type == "DNA":
            return DNASequence(self.seq[::-1])
        else:
            return RNASequence(self.seq[::-1])

    def complement(self):
        "Make complementary DNA sequence"
        new_seq = ''
        for i in self.seq:
            if i == 'A':
                new_seq = new_seq + 'T'
            elif i == 'T':
                new_seq = new_seq + 'A'
            elif i == 'C':
                new_seq = new_seq + 'G'
            elif i == 'G':
                new_seq = new_seq + 'C'
            elif i == 'U':
                new_seq = new_seq + 'A'
        return DNASequence(new_seq)

    def reverse_complement(self):
        "Make reverse complementary DNA sequence"
        rev = self.reverse()
        return rev.complement()


class DNASequence(NucleicAcidSequence):
    def __init__(self, seq):
        super().__init__(seq)
        if self.type == 'RNA':
            raise WrongSequenceTypeError('This is not a DNA')

    def transcribe(self):
        "Transcibe DNA sequence"
        return RNASequence(self.seq.replace('T', 'U'))


class RNASequence(NucleicAcidSequence):
    def __init__(self, seq):
        super().__init__(seq)
        if set(self.seq).issubset(self.rna):  # This part written because by default we consider sequence as dna
            self.type = 'RNA'
        else:
            raise WrongSequenceTypeError('This is not a RNA')


class AminoAcidSequence(BiologicalSequence):
    def __init__(self, seq):
        super().__init__(seq)
        if set(self.seq).issubset(self.aa_list):  # This part written because by default we consider sequence as dna
            self.type = 'Protein'
        else:
            raise WrongSequenceTypeError('This is not a Protein')

    def one_hot_code(self):
        """
        Return 2d np.array(np.float32): peptide in one-hot representation.
        """
        self.prot_oh_encoded = np.zeros((len(self.aa_list), self.len_seq), dtype=np.float32)

        for i in range(self.len_seq):
            for j in range(len(self.aa_list)):
                if self.aa_list[j] == self.seq[i]:
                    self.prot_oh_encoded[j][i] = 1
        return self.prot_oh_encoded

    def one_hot_decode(self):
        """
        Only working with 2d np.array(np.float32)
        """
        ans = ''
        if self.prot_oh_encoded is None:
            raise ProteingEncodingError('Protein have not coded yet. Processed it with Protein(sequrnce).one_hot_code()')
        for i in range(self.prot_oh_encoded.shape[1]):
            for j in range(len(self.aa_list)):
                if self.prot_oh_encoded[j][i] == 1:
                    ans += self.aa_list[j]
        return ans


def filter_fastq(input_fastq: str, output_fastq: str, **kwargs) -> None:
    """
    Filter input fastq files.
    Parameters for filtration could be specified (default): gc_bounds (0,100), length_bounds(0, 2**32), quality_threshold (0).
    Work with fullpath to the input_fastq or if its in working directory.
    Write filtrated reads to ./filtered/output_fastq
    """
    params = process_params(**kwargs)
    good_reads = (
        rec
        for rec in SeqIO.parse(input_fastq, "fastq")
        if (
            statistics.mean(rec.letter_annotations["phred_quality"]) >= params['quality_threshold'] and
            len(rec) > params['length_bounds'][0] and len(rec) < params['length_bounds'][1] and
            GC(rec.seq) > params['gc_bounds'][0] and GC(rec.seq) < params['gc_bounds'][1]
        )
    )
    output_path = os.path.join('filtered', output_fastq)
    if os.path.isdir('filtered'):
        SeqIO.write(good_reads, output_path, "fastq")
    else:
        os.mkdir('filtered')
        SeqIO.write(good_reads, output_path, "fastq")


def process_params(**kwargs) -> dict:
    "Set parameters for selection (default): gc_bounds (0,100), length_bounds(0, 2**32), quality_threshold (0)"
    gc_b = kwargs.get('gc_bounds')
    params = {}
    if gc_b is None:
        gc_b_low = 0
        gc_b_high = 100
    elif isinstance(gc_b, tuple):
        gc_b_low = gc_b[0]
        gc_b_high = gc_b[1]
    else:
        gc_b_low = 0
        gc_b_high = gc_b
    lenght_b = kwargs.get('length_bounds')
    if lenght_b is None:
        lenght_b_low = 0
        lenght_b_high = 2**32
    elif isinstance(lenght_b, tuple):
        lenght_b_low = lenght_b[0]
        lenght_b_high = lenght_b[1]
    else:
        lenght_b_low = 0
        lenght_b_high = lenght_b
    quality_t = kwargs.get('quality_threshold')
    if quality_t is None:
        quality_t = 0
    params = {'gc_bounds': (gc_b_low, gc_b_high), 'length_bounds': (lenght_b_low, lenght_b_high), 'quality_threshold': quality_t}
    return params




parser = argparse.ArgumentParser(
                    prog='Some_dna_rna_sheet fastq filtrator',
                    description='This tool filter input fastq file.\nParameters for filtration could be specified by default:\ngc_bounds 0 100,\nlength_bounds 0 2**32,\nquality_threshold 0. Work with fullpath to the input_fastq or if its in working directory. Write filtrated reads to ./filtered/output_fastq',
                    epilog='Have a good time!')

parser.add_argument('input_fastq', type=str, help='Input fastq  file.')
parser.add_argument('output_fastq', type=str, help='Name of output file')
parser.add_argument('-l', '--lenght', type=int, nargs='+', help='Length bounds for keeping sequence')
parser.add_argument('-g', '--gcbounds', type=int, nargs='+', help='GC content bounds for keeping sequence')
parser.add_argument('-q', '--quality',type=int, help='Quality threshold for keeping sequence')


if __name__ == '__main__':
    # Parse arguments
    args = parser.parse_args()

    # Process lenght bound argument
    if args.lenght is not None:
        if len(args.lenght) == 1:
            args_length_bounds = args.lenght[0]
        elif len(args.lenght) == 2:
            args_length_bounds = tuple(args.lenght)
        else:
            raise ValueError('Too many arguments for length bounds')
    else:
        args_length_bounds = None

    # Process GC content bound argument
    if args.gcbounds is not None:
        if len(args.gcbounds) == 1:
            args_gcbounds = args.gcbounds[0]
        elif len(args.gcbounds) == 2:
            args_gcbounds = tuple(args.gcbounds)
        else:
            raise ValueError('Too many arguments for GC content bounds')
    else:
        args_gcbounds = None

    filter_fastq(input_fastq=args.input_fastq, output_fastq=args.output_fastq, gc_bounds=args_gcbounds, length_bounds=args_length_bounds, quality_threshold=args.quality)
