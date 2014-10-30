#!/usr/bin/env python

import os.path
from StringIO import StringIO
import pytest
import compute_stats

def create_temp_file(tmpdir, filename, data):
    outfile = tmpdir.join(filename).open(mode='w')
    outfile.write(data)
    outfile.close()
    outfile = tmpdir.join(filename).open(mode='r')
    return outfile

@pytest.fixture
def create_target_file():
    data = """hsa-miR-16-5p corresponds with CCNE1 |ENSG00000105173,ENST00000262643,CCNE1,HGNC:26442
hsa-miR-24-3p corresponds with CDK4 |ENSG00000135446,ENST00000257904,CDK4,HGNC:1773
hsa-miR-27a-3p corresponds with FBW7 \n"""
    return StringIO(data)

@pytest.fixture
def create_miranda_file():
    data = """mirna_id         utr_id          score   energy          mirna_hit_start         mirna_hit_end   utr_hit_start   utr_hit_end     aln_length      identity        similarity      aln_mirna       aln_map         aln_utr         seed_mer
hsa-miR-16-5p   ENSG00000105173 100.000000      -15.310000      3       8       409     430     5       100.000000      100.000000      GCGGUUAUAAAUGCACGACGAU                 |||||    CCAGCTGGGCAGGGGGCTGCCC  no_seed_mer
mirna_id         utr_id          score   energy          mirna_hit_start         mirna_hit_end   utr_hit_start   utr_hit_end     aln_length      identity        similarity      aln_mirna       aln_map         aln_utr         seed_mer
hsa-miR-16-5p   ENSG00000159216 151.000000      -15.690000      2       20      1453    1474    18      72.222222       83.333333       GCGGUUAUAAAUGCACGACGAU     ||| || | :|||:||||   CAACAACATCTTTGTGTTGCTT  no_seed_mer
mirna_id         utr_id          score   energy          mirna_hit_start         mirna_hit_end   utr_hit_start   utr_hit_end     aln_length      identity        similarity      aln_mirna       aln_map         aln_utr         seed_mer
hsa-miR-24-3p   ENSG00000001617 123.000000      -17.340000      2       18      1       18      16      75.000000       75.000000       GACAAGGACGACUUGACUCGGU       || ||||  ||| |||   ---GGCCAGCTG-CCTGTGCCT  no_seed_mer
"""
    return StringIO(data)

@pytest.fixture
def create_rnahybrid_file():
    data = """ENSG00000079432|ENST00000160740:605:hsa-miR-24-3p:22:-25.5:0.312022:281:C    GAGGG  GCAA      GGGUGGUGCAGGCC        C: CUGU     CU    GCUGGA              UUGGGCCA : GACA     GG    CGACUU              GACUCGGU :     A      A
ENSG00000082458|ENST00000194900:3221:hsa-miR-24-3p:22:-29.2:0.410775:2714:    A    AAAUGUGUG            C:     CCUG         UGGGCUGAGCCA :     GGAC         ACUUGACUCGGU :GACAA    G
ENSG00000105971|ENST00000222693:2420:hsa-miR-24-3p:22:-25.8:0.805509:1239:A     G        AAAUCUCAUUAUUU        G: CUGUU CU UUGAA              CUGGGCCG : GACAA GA GACUU              GACUCGGU :      G  C
ENSG00000001617|ENST00000002829:960:hsa-miR-16-5p:22:-21.2:0.982341:234:A          UG   GG  C : GCCAG    U  GUG  GC  : CGGUU    A  CAC  CG  :G     AUAA UG   GA  AU
ENSG00000166851|ENST00000300093:304:hsa-miR-24-3p:22:-26.8:0.068572:105:    C  GGUG       A     G:     CC    GCUGGGC GAGCU :     GG    CGACUUG CUCGG :GACAA  A          A     U
"""
    return StringIO(data)

@pytest.fixture
def create_microtar_file():
    data = ""

def test_find_known_targets_file(create_target_file):
    # 
    known_targets = compute_stats.find_known_targets(create_target_file)
    assert len(known_targets) == 2
    assert 'hsa-miR-16-5p' in known_targets and 'hsa-miR-24-3p' in known_targets
    assert not 'hsa-miR-27a-3p' in known_targets

def test_find_miranda_targets(create_miranda_file):
    fileobj = create_miranda_file
    (predicted_targets, predicted_non_targets) = compute_stats.find_miranda_targets(fileobj)
    assert len(predicted_targets) == 1
    assert len(predicted_non_targets) == 2
    assert 'hsa-miR-16-5p' in predicted_targets
    assert 'ENSG00000159216' in predicted_targets['hsa-miR-16-5p']
    assert 'hsa-miR-16-5p' in predicted_non_targets and 'hsa-miR-24-3p' in predicted_non_targets

def test_find_rnahybrid_targets(create_rnahybrid_file):
    (predicted_targets, predicted_non_targets) = compute_stats.find_rnahybrid_targets(create_rnahybrid_file)
    assert len(predicted_targets) == 1
    assert len(predicted_non_targets) == 2
    assert 'hsa-miR-24-3p' in predicted_targets
    assert 'ENSG00000166851' in predicted_targets['hsa-miR-24-3p']
    assert 'hsa-miR-16-5p' in predicted_non_targets and 'hsa-miR-24-3p' in predicted_non_targets

# def test_find_microtar_targets(create_microtar_file):
#     (predicted_targets, predicted_non_targets) = compute_stats.find_rnahybrid_targets(create_microtar_file)
#     assert len(predicted_targets) == 1
#     assert len(predicted_non_targets) == 2
#     assert 'hsa-miR-16-5p' in predicted_targets
#     assert 'ENSG00000171791' in predicted_targets['hsa-miR-16-5p']
#     assert 'hsa-miR-16-5p' in predicted_non_targets and 'hsa-miR-24-3p' in predicted_non_targets

def test_compute_stats(create_target_file, create_miranda_file):
    known_targets = compute_stats.find_known_targets(create_target_file)
    (predicted_targets, predicted_non_targets) = compute_stats.find_miranda_targets(create_miranda_file)
    (TP, FP, TN, FN, sensitivity, specificity) = compute_stats.compute_stats(known_targets, predicted_targets, predicted_non_targets)
    assert TP == 0
    assert FP == 1
    assert TN == 1
    assert FN == 2
    assert sensitivity == 0
    assert specificity == 50.0