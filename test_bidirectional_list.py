import pytest

from bidirectional_list import BidirectionalList as BiList, BidirectionalNode as BiNode


def test_empty_list_is_empty():
    l = BiList()
    assert len(l) == 0
    assert list(l) == []

def test_list_initialization():
    l = BiList([42, 43, 44])
    assert len(l) == 3
    assert [el.value for el in l] == [42, 43, 44]

def test_empty_out_of_bounds_indexing():
    l = BiList()
    with pytest.raises(KeyError):
        l[0]

def test_out_of_bounds_indexing():
    l = BiList([42])
    assert l[0].value == 42
    with pytest.raises(KeyError):
        l[1]
        l[-1]

def test_single_element_reverse_indexing():
    l = BiList([42])
    assert l[0].value == l[-1].value

def test_reverse_indexing():
    l = BiList([42, 43, 44])
    assert l[0].value == l[-3].value
    assert l[1].value == l[-2].value
    assert l[2].value == l[-1].value

def test_append_empty():
    l = BiList()
    l.append(42)
    assert len(l) == 1
    assert [el.value for el in l] == [42]

def test_prepend_empty():
    l = BiList()
    l.prepend(42)
    assert len(l) == 1
    assert [el.value for el in l] == [42]

def test_append():
    l = BiList()
    l.append(42)
    l.append(43)
    l.append(44)
    l.append(45)
    assert len(l) == 4
    assert [el.value for el in l] == [42, 43, 44, 45]

def test_prepend():
    l = BiList()
    l.prepend(42)
    l.prepend(43)
    l.prepend(44)
    l.prepend(45)
    assert len(l) == 4
    assert [el.value for el in l] == [45, 44, 43, 42]

def test_insert_after():
    l = BiList()
    node = l.append(42)
    l.append(44)
    l.append(45)
    l.insert_after(node, 43)
    assert [el.value for el in l] == [42, 43, 44, 45]

def test_insert_before():
    l = BiList()
    l.append(42)
    node = l.append(44)
    l.append(45)
    l.insert_before(node, 43)
    assert [el.value for el in l] == [42, 43, 44, 45]

def test_remove():
    l = BiList()
    l.append(42)
    node = l.append(43)
    l.append(44)
    l.append(45)
    l.remove(node)
    assert len(l) == 3
    assert [el.value for el in l] == [42, 44, 45]

def test_remove_first():
    l = BiList()
    node = l.append(42)
    l.append(43)
    l.append(44)
    l.append(45)
    l.remove(node)
    assert len(l) == 3
    assert [el.value for el in l] == [43, 44, 45]

def test_remove_last():
    l = BiList()
    l.append(42)
    l.append(43)
    l.append(44)
    node = l.append(45)
    l.remove(node)
    assert len(l) == 3
    assert [el.value for el in l] == [42, 43, 44]

def test_empty_after_remove():
    l = BiList()
    node = l.append(42)
    l.remove(node)
    assert len(l) == 0
    assert list(l) == []
