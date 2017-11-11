import collections

Element = collections.namedtuple("Element", ["atom_number", "symbol"])
elems = [
    Element(1, "H"),
    Element(2, "He"),
    #
    Element(3, "Li"),
    Element(4, "Be"),
    Element(5, "B"),
    Element(6, "C"),
    Element(7, "N"),
    Element(8, "O"),
    Element(9, "F"),
    Element(10, "Ne"),
    #
    Element(11, "Na"),
    Element(12, "Mg"),
    Element(13, "Al"),
    Element(14, "Si"),
    Element(15, "P"),
    Element(16, "S"),
    Element(17, "Cl"),
    Element(18, "Ar"),
    #
    Element(19, "K"),
    Element(20, "Ca"),
]


def by_atomic_number(num):
    return [e for e in elems if e.atom_number == num][0]


def find(atom):
    try:
        if isinstance(atom, int):
            return [e for e in elems if e.atom_number == atom][0]
        else:
            return [e for e in elems if e.symbol == atom][0]
    except IndexError:
        raise KeyError(atom)
