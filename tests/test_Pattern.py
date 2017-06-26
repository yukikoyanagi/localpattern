from .. import Pattern


def test_equality():
    # test equality of Segment objects
    a = Pattern.Segment()
    b = Pattern.Segment()
    for i in range(6):
        a.append(Pattern.Atom(i, chr(69+i)))
        b.append(Pattern.Atom(100+i, chr(75+i)))
        c = b[:]
    c[1] = Pattern.Atom(101, chr(66))

    s = Pattern.Bond('h', 3, 104, False)
    t = Pattern.Bond('h', 4, 104, False)
    u = Pattern.Bond('t', 3, 104, False)
    v = Pattern.Bond('g', 3, 104, True)

    P1 = Pattern.Pattern(segments=[a,b], bonds=[s], rotation=123)
    P2 = Pattern.Pattern(segments=[a,c], bonds=[s], rotation=123)
    P3 = Pattern.Pattern(segments=[a,b], bonds=[t], rotation=123)
    P4 = Pattern.Pattern(segments=[a,b], bonds=[u], rotation=123)
    P5 = Pattern.Pattern(segments=[a,b], bonds=[v], rotation=123)
    P6 = Pattern.Pattern(segments=[a,b], bonds=[s], rotation=12)

    assert P1 != P2
    assert P1 != P3
    assert P1 != P4
    assert P1 != P5
    assert P1 == P6
