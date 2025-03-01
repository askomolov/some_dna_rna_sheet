# Some DNA RNA sheet 
**Some DNA RNA sheet** - python library, developed to work with DNA and RNA sequences, READs from **.fastq**, sequences from **.fasta** and parse Blast output.

## Library installation
To install this library on your computer you can: 

1. Download this repository to your computer using the following command `git clone <link>`.

2. Download main scripts manually **some_dna_rna_sheet.py**, **bio_files_processor.py** and the folder **modules**.

In both cases, the repository or scripts must be in the working directory.


## Main functionality

This library can:

1. Work with DNA and RNA (DNA transcription, make complementary DNA sequence, make reverse sequence, make reverse complementary DNA sequence). Script <u>some_dna_rna_sheet.py</u>

Usage:
```python
import some_dna_rna_sheet

some_dna_rna_sheet.dna_rna_tools(<Sequences>, <action>)
```
Sequences are specified in the format of a string or a list of strings. The case is not important. An action is specified with one of the strings: 'transcribe', 'reverse','complement','reverse_complement'.

2. Filter READs from a file in the following format **.fastq** based on input parameters or set by default. Script <u>some_dna_rna_sheet.py</u>

Usage:
```python
import some_dna_rna_sheet

some_dna_rna_sheet.filter_fastq(input_fastq = <name_of_input_file.fastq>, output_fastq = <name_of_output_file.fastq>,  <Options>)
```
The input is a file in the format **<name_of_input_file>.fastq**. You can also specify the name of the output file. By default: **<name_of_input_file>_output.fastq**. The output file is saved in the folder **./filtered**

Options: `gc_bounds` - GC filtration range (by default (0, 100)); `length_bounds` - length filtration range (by default (0, 2**32)); `quality_threshold` - threshold value of average READ quality for filtration, by default 0 (scale - phred33).

3. Convert  multiline **.fasta** into oneline. Script <u>bio_files_processor.py</u>

Usage:
```python
import bio_files_processor

bio_files_processor.convert_multiline_fasta_to_oneline(input_fasta = <name_of_input_file.fasta>, output_fasta = <name_of_output_file.fasta>)
```

4. Parse blast output and returns the descriptions of proteins for each Query. Script <u>bio_files_processor.py</u>

Usage:
```python
import bio_files_processor

bio_files_processor.parse_blast_output(input_file = <name_of_input_file.txt>, output_file = <name_of_output_file.txt>)
```

More detailed information about the functions can be obtained by calling the help with the command `help()`

```python
import some_dna_rna_sheet

help(some_dna_rna_sheet)
```


## Contacts
This library was developed as part of the [homework №4](https://github.com/Python-BI-2024-25/course_materials/tree/main/Homeworks/HW4_Modules) and [homework №5](https://github.com/Python-BI-2024-25/course_materials/tree/main/Homeworks/HW5_Files)  by student of [Bioinformatics Institute](https://bioinf.me) 2024/25 year - Komolov Aleksandr.

If you have any questions, please contact us at: askomolov@mail.ru
