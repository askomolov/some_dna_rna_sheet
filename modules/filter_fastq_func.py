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


def fetch_seqs_from_file(input_file: str) -> dict:
    "Read FASTQ file and fetch reads"
    seqs = {}
    with open(input_file) as file:
        file_cont = file.read()
    pre_seqs = file_cont.split('\n')
    for i in range(0, len(pre_seqs) - 1, 4):
        seqs[pre_seqs[i]] = (pre_seqs[i + 1], pre_seqs[i + 3])
    return seqs


def seqs_choose(seqs: dict, params: dict) -> dict:
    "Sequence choose by parameters"
    new_seqs = {}
    for name, seq in seqs.items():
        if check_gc(seq[0], params['gc_bounds']) and check_quality(seq[1], params['quality_threshold']) and check_lenght(seq[0], params['length_bounds']):
            new_seqs[name] = seq
    return new_seqs


def check_gc(seq: str, gc_comp: tuple) -> bool:
    "Calculate GC composition of sequence (in%) and check if it in set parameters"
    gc_content = (seq.count('G') + seq.count('C')) / len(seq) * 100
    if gc_content <= gc_comp[1] and gc_content >= gc_comp[0]:
        return True
    return False


def check_quality(seq_q: str, quality_t: float) -> bool:
    "Calculate average quality and check if it more then set quality_threshold"
    quality = 0
    for i in range(0, len(seq_q)):
        quality += ord(seq_q[i]) - 33
    quality = quality / len(seq_q)
    if quality >= quality_t:
        return True
    return False


def check_lenght(seq: str, length_b: tuple) -> bool:
    "Check if it in lenght of sequence in set parameters"
    if len(seq) >= length_b[0] and len(seq) <= length_b[1]:
        return True
    return False


def write_seqs_to_file(seqs: dict, output_file: str) -> None:
    "Write chosen reads to output FASTQ file"
    with open(output_file, "a") as file:
        for name, seq in seqs.items():
            name_q = name.replace('@', '+')
            file.write(name + '\n')
            file.write(seq[0] + '\n')
            file.write(name_q + '\n')
            file.write(seq[1] + '\n')
    return None
