import modules.dna_rna_tools_func as drt
import modules.filter_fastq_func as fifa


def run_dna_rna_tools(*args):
    *seqs, proc = args
    new_seqs = []
    if drt.input_check(seqs) == 0:
        return None
    for i in range(0, len(seqs)):
        if proc == 'transcribe':
            if drt.is_dna(seqs[i]) is False:
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


def filter_fastq(seqs, **kwargs):
    params = fifa.check_params(**kwargs)  # Записывыем все параметры фильтраций (введенные или по умолчанию) в один словарь
    act_seqs = fifa.seqs_choose(seqs, params)  # Выбор последовательностей, подходящих по параметрам
    return act_seqs
