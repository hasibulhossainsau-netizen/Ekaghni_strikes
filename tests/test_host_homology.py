import pandas as pd
import pytest

from filter_host_homology import BLAST_COLUMNS, filter_blast_hits


THRESHOLDS = {"identity": 80, "coverage": 70, "evalue": 1e-3}


def test_host_filter_applies_all_thresholds_and_deduplicates_regions():
    results = pd.DataFrame([
        ["protein_1", "dog_A", 95.0, 10, 1, 10, 1e-8, 90.0],
        ["protein_1", "dog_B", 91.0, 10, 1, 10, 1e-7, 85.0],
        ["protein_2", "dog_C", 79.0, 10, 1, 10, 1e-8, 90.0],
        ["protein_3", "dog_D", 95.0, 10, 1, 10, 1e-8, 60.0],
        ["protein_4", "dog_E", 95.0, 10, 1, 10, 1e-2, 90.0],
    ], columns=BLAST_COLUMNS)

    filtered = filter_blast_hits(results, THRESHOLDS)

    assert len(filtered) == 1
    assert filtered.iloc[0]["protein_id"] == "protein_1"
    assert filtered.iloc[0]["host_match"] == "dog_A"


def test_host_filter_requires_blast_schema():
    with pytest.raises(ValueError, match="qseqid"):
        filter_blast_hits(pd.DataFrame({"query": ["protein_1"]}), THRESHOLDS)
