import some_dna_rna_sheet as sdrs
import pytest
import os
import numpy as np


# Test DNA and RNA basic functions
def test_reverse():
    inp = 'ATgC'
    target = 'CGTA'
    assert target == sdrs.DNASequence(inp).reverse().seq


def test_transcribe():
    inp = 'GGGCAT'
    target = 'GGGCAU'
    assert target == sdrs.DNASequence(inp).transcribe().seq


def test_complement_RNA():
    inp = 'UAGCUU'
    target = 'ATCGAA'
    assert target == sdrs.RNASequence(inp).complement().seq


def test_complement_DNA():
    inp = 'TAGCTT'
    target = 'ATCGAA'
    assert target == sdrs.DNASequence(inp).complement().seq


def test_reverse_complement_RNA():
    inp = 'UAGCUU'
    target = 'AAGCTA'
    assert target == sdrs.RNASequence(inp).reverse_complement().seq


def test_reverse_complement_DNA():
    inp = 'TAGCTT'
    target = 'AAGCTA'
    assert target == sdrs.DNASequence(inp).reverse_complement().seq


# For protein sequence one-hot encoding
def test_one_hot_code():
    inp = 'LLLLLLKD'
    target = [
       [0., 0., 0., 0., 0., 0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 1.],
       [0., 0., 0., 0., 0., 0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 1., 0.],
       [1., 1., 1., 1., 1., 1., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0.]
    ]
    assert np.array_equal(sdrs.AminoAcidSequence(inp).one_hot_code(), target)


# Check output file exists
@pytest.fixture
def input_data():
    return './example_data/example_fastq.fastq'


@pytest.fixture
def tmp_file():
    file_path = 'tmp.fasta'
    output_path = os.path.join('filtered', file_path)
    yield file_path
    if os.path.exists(output_path):
        os.remove(output_path)


def test_write_fasta_exists(input_data, tmp_file):
    sdrs.filter_fastq(input_data, tmp_file)
    output_path = os.path.join('filtered', tmp_file)
    assert os.path.exists(output_path)
