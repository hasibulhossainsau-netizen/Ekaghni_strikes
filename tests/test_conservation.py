import pytest

from calculate_conservation import calculate_conservation


def test_conservation_scores_gaps_and_rolling_window():
    result = calculate_conservation(["A-C", "ATC"], window_size=2)

    assert result["consensus_aa"].tolist() == ["A", "T", "C"]
    assert result["raw_score"].tolist() == [1.0, 0.5, 1.0]
    assert result["window_2mer_avg"].tolist() == [1.0, 0.75, 0.75]


def test_all_gap_column_has_zero_score_and_unknown_consensus():
    result = calculate_conservation(["A--", "A--"], window_size=1)

    assert result.loc[1, "consensus_aa"] == "X"
    assert result.loc[1, "raw_score"] == 0.0


@pytest.mark.parametrize(
    ("sequences", "window"),
    [([], 3), (["AA", "A"], 2), (["AA"], 0)],
)
def test_conservation_rejects_invalid_alignment(sequences, window):
    with pytest.raises(ValueError):
        calculate_conservation(sequences, window)
