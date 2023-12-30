from rbc2.reaction_evaluation.complexity import get_complexity


def test_can_get_complexity():
    smi = 'CCC=O'
    score = get_complexity(smi)
    assert score == 1.2394

