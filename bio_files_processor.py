def convert_multiline_fasta_to_oneline(input_fasta: str, output_fasta: str = '') -> None:
    """
    Convert sequences in input_fasta from multilines format to oneline format and write the results to output_fasta.
    input_fasta - name of the input file. Should in *.fasta format
    output_fasta - name of the output file. Default: input_fasta with '_output' at the end.
    """
    file_in = open(input_fasta)
    if len(output_fasta) == 0:
        output_fasta = input_fasta.replace('.fasta', '_output.fasta')
    file_out = open(output_fasta, "a")
    file_cont = file_in.readline()
    file_out.write(file_cont)
    file_cont = file_in.readline()
    while (len(file_cont) > 0):
        if file_cont[0] == '>':
            file_out.write('\n' + file_cont)
        else:
            file_out.write(file_cont[:-1])
        file_cont = file_in.readline()
    file_in.close()
    file_out.close()
    return None


def parse_blast_output(input_file: str, output_file: str = '') -> None:
    """
    Parse blast output returns the descriptions of proteins for each Query.
    input_file - name of the input file. Should in *.txt format
    output_fasta - name of the output file. Default: input_fasta with '_output' at the end.
    """
    with open(input_file) as file_in:
        output = []
        for line in file_in:
            if line.startswith('Description'):
                descr_line = file_in.readline()
                descr_line_sp = descr_line.split('[')  # Almost always works to split line in 2 part
                if len(descr_line_sp[0]) == len(descr_line):  # If '[' didn't split the line, then use '...' as it in example file
                    descr_line_sp = descr_line.split('...')
                output.append(descr_line_sp[0])
    output.sort(key=str.casefold)
    if len(output_file) == 0:
        output_file = input_file.replace('.txt', '_output.txt')
    with open(output_file, "a") as file_out:
        for i in output:
            file_out.write(i + '\n')
    return None
