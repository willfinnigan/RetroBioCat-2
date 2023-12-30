from rbc2.mcts.mcts_loop.selection import selection, ucb_calc, find_most_promising_child
from rbc2.mcts.tree_node import create_root, MCTS_Node


def test_first_selection_returns_root():
    root = create_root('CCCC=O')
    node = selection(root, 2)
    assert root == node


def test_ucb_calc_returns_1_for_first_new_node():
    score = ucb_calc(node_score=1, node_visits=1, exploration=2, parent_visits=1)
    assert score == 1


def test_find_most_promising_child_returns_node_with_highest_value():
    root = create_root('CCCC=O')
    child_1 = MCTS_Node(parent=root, value=1)
    child_2 = MCTS_Node(parent=root, value=2)
    child_3 = MCTS_Node(parent=root, value=0.5)
    root.children = [child_1, child_2, child_3]

    promising_child = find_most_promising_child(root, 2)
    assert promising_child == child_2

def test_if_child_2_is_terminal_then_child_1_is_returned():
    root = create_root('CCCC=O')
    child_1 = MCTS_Node(parent=root, value=1)
    child_2 = MCTS_Node(parent=root, value=2, terminal=True)
    child_3 = MCTS_Node(parent=root, value=0.5)
    root.children = [child_1, child_2, child_3]

    promising_child = find_most_promising_child(root, 2)
    assert promising_child == child_1

def test_if_all_children_are_terminal_or_searched_then_root_becomes_fully_searched():
    root = create_root('CCCC=O')
    child_1 = MCTS_Node(parent=root, value=1, terminal=True)
    child_2 = MCTS_Node(parent=root, value=2, terminal=True)
    child_3 = MCTS_Node(parent=root, value=0.5, terminal=True)
    root.children = [child_1, child_2, child_3]

    promising_child = find_most_promising_child(root, 2)
    assert promising_child == None
    assert root.fully_searched == True

def test_first_selection_with_children_returns_child_2():
    root = create_root('CCCC=O')
    child_1 = MCTS_Node(parent=root, value=1, depth=1)
    child_2 = MCTS_Node(parent=root, value=2, depth=1)
    child_3 = MCTS_Node(parent=root, value=0.5, depth=1)
    root.children = [child_1, child_2, child_3]

    node = selection(root, 2)
    assert node == child_2

def test_second_selection_with_children_returns_child_2():
    root = create_root('CCCC=O')
    child_1 = MCTS_Node(parent=root, value=1, depth=1)
    child_2 = MCTS_Node(parent=root, value=2, depth=1)
    child_3 = MCTS_Node(parent=root, value=0.5, depth=1)
    root.children = [child_1, child_2, child_3]
    child_4 = MCTS_Node(parent=child_2, value=1, depth=2)
    child_2.children = [child_4]

    node = selection(root, 2)
    assert node == child_4