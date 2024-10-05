def check_params(**kwargs) -> dict:   # Функция проверяет введены ли параметры фильтрации. Если нет, то задаёт значения по умолчанию.
    "Set parameters for selection (default): gc_bounds (0,100), length_bounds(0, 2**32), quality_threshold (0)"
    gc_b = kwargs.get('gc_bounds')
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
    params1 = {'gc_bounds': (gc_b_low, gc_b_high), 'length_bounds': (lenght_b_low, lenght_b_high), 'quality_threshold': quality_t}
    return params1


def seqs_choose(seqs: dict, params: dict) -> dict:
    "Sequence choose by parameters"
    new_seqs = {}
    for name, seq in seqs.items():
        gc_comp = gc_comp_calc(seq[0])
        avg_quality = quality_calc(seq[1])
        len_seq = len(seq[0])
        if gc_comp <= params['gc_bounds'][1] and gc_comp >= params['gc_bounds'][0] and len_seq <= params['length_bounds'][1] and len_seq >= params['length_bounds'][0] and avg_quality >= params['quality_threshold']:
            new_seqs[name] = seq
    return new_seqs


def gc_comp_calc(seq: str) -> float:
    "Calculate GC composition of sequence (in%)"
    return (seq.count('G')+seq.count('C'))/len(seq)*100


def quality_calc(seq_q: str) -> float:
    "Calculate average quality"
    quality = 0
    for i in range(0, len(seq_q)):
        quality += ord(seq_q[i])-33
    quality = quality/len(seq_q)
    return quality
