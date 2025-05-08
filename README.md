# Some DNA RNA sheet 
**Some DNA RNA sheet** - python OOP-based library, developed to work with DNA, RNA and protein sequences, READs in **.fastq** format, sequences from **.fasta** and it could also parse **Blast** output.

## Library installation
To install this library on your computer you can: 

1. Download this repository to your computer using the following command `git clone <link>`.

2. Download main scripts manually **some_dna_rna_sheet.py**, **bio_files_processor.py**.

In both cases, the repository or scripts must be in the working directory.


## Main functionality

This library have 5 type of classes:
- `BiologicalSequence()` with following attributes: `.seq`, `.len_seq`, `.type`
- `NucleicAcidSequence(BiologicalSequence)` with following methods: `reverse()`, `complement()`, `reverse_complement()`
- `DNASequence(NucleicAcidSequence)` with following methods: `transcribe()`
- `RNASequence(NucleicAcidSequence)`
- `AminoAcidSequence(BiologicalSequence)` with following methods: `one_hot_code()`, `one_hot_decode()`


This  library can:

1. Filter READs from a file in the following format **.fastq** based on input parameters or set by default. Script <u>some_dna_rna_sheet.py</u>

Usage:
```python
import some_dna_rna_sheet

some_dna_rna_sheet.filter_fastq(input_fastq, output_fastq,  **kwargs)
```
The input is a file in the format **<name_of_input_file>.fastq**. You should also specify the name of the output file. The output file is saved in the folder **./filtered** It also creates `log.log` file with useful information.

\**kwargs: `gc_bounds` - GC filtration range (by default (0, 100)); `length_bounds` - length filtration range (by default (0, 2**32)); `quality_threshold` - threshold value of average READ quality for filtration, by default 0 (scale - phred33).


`some_dna_rna_sheet.filter_fastq` also works from command line. Just do the following:
```bash
python some_dna_rna_sheet.py -h
```

2. Convert  multiline **.fasta** into oneline. Script <u>bio_files_processor.py</u>

Usage:
```python
import bio_files_processor

bio_files_processor.convert_multiline_fasta_to_oneline(input_fasta = <name_of_input_file.fasta>, output_fasta = <name_of_output_file.fasta>)
```

3. Parse blast output and returns the descriptions of proteins for each Query. Script <u>bio_files_processor.py</u>

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


## Requirements

To run python scripts some other packages should be installed. Complete list of dependencies presented in `requirementes.txt`. 

It's reccomended to use new environment to work with library. 

Just do the folowwing commands:

`conda create --name test_env --file requirements.txt`

`conda activate test_env`

List of packages contains some to work with jupyter lab. 
To make new kernel use:
`ipython kernel install --user --name=<any_name_for_kernel>`

There is a script for `pytest` to test the correct performance.
Run it with:
```bash
pytest
```

*This library was tested on MacBook Pro (M1).*


## Contacts
This library was developed as part of the [homework №4](https://github.com/Python-BI-2024-25/course_materials/tree/main/Homeworks/HW4_Modules), [homework №5](https://github.com/Python-BI-2024-25/course_materials/tree/main/Homeworks/HW5_Files) and [homework №15](https://github.com/Python-BI-2024-25/course_materials/blob/main/Homeworks/HW15_BioOOP.md)  by student of [Bioinformatics Institute](https://bioinf.me) 2024/25 year - Komolov Aleksandr.

If you have any questions, please contact via mail: askomolov@mail.ru
