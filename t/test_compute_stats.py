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
def create_rnahybrid_file(tmpdir):
    data = """utr_id   utr_length      mirna_id        mirna_length    energy          p-value         hit_position
ENSG00000105173 531     hsa-miR-16-5p   22      -24.8   0.019500        459
utr_id   utr_length      mirna_id        mirna_length    energy          p-value         hit_position
target too long  ENSG00000159216
utr_id   utr_length      mirna_id        mirna_length    energy          p-value         hit_position
ENSG00000001617 960     hsa-miR-24-3p   22      -29.1   0.015504        10
"""
    return StringIO(data)

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
    assert 'hsa-miR-16-5p' in predicted_targets
    assert 'ENSG00000159216' in predicted_targets['hsa-miR-16-5p']
    assert 'hsa-miR-16-5p' in predicted_non_targets and 'hsa-miR-24-3p' in predicted_non_targets

def test_find_microtar_targets(create_microtar_file):
    (predicted_targets, predicted_non_targets) = compute_stats.find_rnahybrid_targets(create_microtar_file)
    assert len(predicted_targets) == 1
    assert len(predicted_non_targets) == 2
    assert 'hsa-miR-16-5p' in predicted_targets
    assert 'ENSG00000159216' in predicted_targets['hsa-miR-16-5p']
    assert 'hsa-miR-16-5p' in predicted_non_targets and 'hsa-miR-24-3p' in predicted_non_targets

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