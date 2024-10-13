import modules.dna_rna_tools_func as drt
import modules.filter_fastq_func as fifa
import os


def run_dna_rna_tools(*args):
    """
    Script wth some tools for DNA and RNA sequences
    Input: One or more DNA or RNA sequences and one action from the list: transcribe, reverse, complement, reverse_complement
    If no DNA or RNA sequences or not nucleotide sequence in input return None.
    If action is:
    1) transcribe - function transcibe all DNA sequence. If RNA in input return None.
    2) reverse - function reverse nucleotide sequence.
    3) complement - function make complementary DNA sequence.
    4) reverse_complement - function make reverse complementary DNA sequence.
    5) else - function return None.
    """
    *seqs, proc = args
    new_seqs = []
    if not drt.input_check(seqs):
        return None
    for i in range(0, len(seqs)):
        if proc == 'transcribe':
            if not drt.is_dna(seqs[i]):
                return print(' Couldn\'t transcribe RNA')
            new_seqs.append(drt.transcribe(seqs[i]))
        elif proc == 'reverse':
            new_seqs.append(drt.reverse(seqs[i]))
        elif proc == 'complement':
            new_seqs.append(drt.complement(seqs[i]))
        elif proc == 'reverse_complement':
            new_seqs.append(drt.reverse_complement(seqs[i]))
        else:
            return print('Unknown procedure')
    if len(seqs) == 1:
        return new_seqs[0]
    else:
        return new_seqs


def filter_fastq(input_fastq: str, output_fastq: str, **kwargs) -> None:
    """
    Filter input fastq files.
    Parameters for filtration could be specified (default): gc_bounds (0,100), length_bounds(0, 2**32), quality_threshold (0).
    Work with fullpath to the input_fastq or if its in working directory.
    Write filtrated reads to ./filtered/output_fastq
    """
    params = fifa.process_params(**kwargs)
    seqs = fifa.fetch_seqs_from_file(input_fastq)
    act_seqs = fifa.seqs_choose(seqs, params)
    output_path = os.path.join('filtered', output_fastq)
    if os.path.isdir('filtered'):
        fifa.write_seqs_to_file(act_seqs, output_path)
    else:
        os.mkdir('filtered')
        fifa.write_seqs_to_file(act_seqs, output_path)
    return None
