import pandas as pd
import pytest
from assemble_vaccine_chimera import get_peptide_list

from final_weighted_prioritization import rank_candidates


WEIGHTS = {
    "mhc_binding": 0.25,
    "mhc2_binding": 0.20,
    "processing_efficiency": 0.10,
    "conservation": 0.15,
    "immunogenicity": 0.20,
    "safety": 0.10,
}


def inputs():
    mhc = pd.DataFrame([
        {"protein_id": "protein_1", "peptide": "AAAA", "start": 1, "end": 4, "percentile_rank": 5.0},
        {"protein_id": "protein_2", "peptide": "BBBB", "start": 5, "end": 8, "percentile_rank": 25.0},
    ])
    topology = pd.DataFrame([
        {"peptide": "AAAA", "topology_score": 0.8},
        {"peptide": "BBBB", "topology_score": 0.6},
    ])
    conservation = pd.DataFrame({
        "position": range(1, 9),
        "window_3mer_avg": [0.9, 0.9, 0.9, 0.9, 0.7, 0.7, 0.7, 0.7],
    })
    blast = pd.DataFrame([
        {"protein_id": "None", "host_match": "None"},
    ])
    ensemble = pd.DataFrame([
        {"peptide": "AAAA", "mhc2_score": 0.8, "proc_score": 0.6, "ensemble_immunogenicity_score": 0.7},
        {"peptide": "BBBB", "mhc2_score": 0.6, "proc_score": 0.4, "ensemble_immunogenicity_score": 0.5},
    ])
    toxicity = pd.DataFrame([
        {"peptide": "AAAA", "toxicity_prediction": "Non-Toxic"},
        {"peptide": "BBBB", "toxicity_prediction": "Non-Toxic"},
    ])
    allergenicity = pd.DataFrame([
        {"peptide": "AAAA", "allergenicity_prediction": "Non-Allergen"},
        {"peptide": "BBBB", "allergenicity_prediction": "Non-Allergen"},
    ])
    return mhc, topology, conservation, blast, ensemble, toxicity, allergenicity


def test_ranking_calculates_weighted_scores_and_orders_candidates():
    result = rank_candidates(*inputs(), WEIGHTS, "window_3mer_avg")

    assert result["peptide"].tolist() == ["AAAA", "BBBB"]
    assert result["prioritization_rank"].tolist() == [1, 2]
    assert result.loc[0, "mhc1_score"] == pytest.approx(0.95)
    assert result.loc[0, "final_vaccine_score"] == pytest.approx(0.8325)


def test_ranking_rejects_missing_core_columns():
    frames = list(inputs())
    frames[0] = frames[0].drop(columns=["percentile_rank"])

    with pytest.raises(ValueError, match="percentile_rank"):
        rank_candidates(*frames, WEIGHTS, "window_3mer_avg")


def test_construct_input_requires_the_canonical_peptide_column():
    with pytest.raises(ValueError, match="peptide"):
        get_peptide_list(pd.DataFrame({"Sequence": ["AAAA"]}), "peptide", 5, "MHC-I")
