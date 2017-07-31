from .. import Pattern


def test_equality():
    # test equality of Segment objects
    a = Pattern.Segment()
    b = Pattern.Segment()
    for i in range(6):
        a.append(i)
        b.append(100+i)
        c = b[:]
    c[-1] = 110

    s = Pattern.Bond('h', 3, 104, False, None)
    t = Pattern.Bond('h', 4, 104, False, None)
    u = Pattern.Bond('t', 3, 104, False, 4)
    v = Pattern.Bond('g', 3, 104, True, None)

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

    P1.residue = 'RSPB'
    P6.residue = 'RSPP'
    assert P1 != P6
    P6.residue = 'RSPB'
    assert P1 == P6


def test_lieswithin():
    p = Pattern.Pattern()
    p.segments = [Pattern.Segment([i]) for i in xrange(0, 101, 10)]
    p.bonds = [Pattern.Bond('T', i, i+10, True, 3) for i in xrange(0, 101, 10)]

    assert p.lieswithin(p.segments[-1], p.bonds[0], 0, 2) is False
    assert p.lieswithin(p.segments[-1], p.bonds[0], 0, 10)

    q = Pattern.Pattern()
    q.segments = [Pattern.Segment(range(i, i+6)) for i in xrange(0, 31, 10)]
    a = Pattern.Bond('H', 2, 12, False, None)
    b = Pattern.Bond('T', 14, 24, False, 2)
    c = Pattern.Bond('H', 21, 31, True, None)
    q.bonds = [a, b, c]

    assert q.lieswithin(q.segments[-1], a, 1, 1)
    assert not q.lieswithin(q.segments[-1], a, 1, 0)

    q.bonds = [a, c]
    assert not q.lieswithin(q.segments[-1], a, 1, 1)


def test_dist():
    p = Pattern.Pattern()
    p.segments = [Pattern.Segment(range(i, i+8)) for i in xrange(0, 21, 10)]
    a = Pattern.Bond('H', 2, 12, False, None)
    b = Pattern.Bond('T', 14, 24, False, 3)
    c = Pattern.Bond('H', 7, 15, True, None)
    p.bonds = [a, b]

    assert p.dist(7, 23) == 10
    assert p.dist(7, 23, hlimit=0, tlimit=1) == float('inf')
    assert p.dist(7, 23, hlimit=1, tlimit=0) == float('inf')

    p.bonds.append(c)
    assert p.dist(7, 23) == 4

    p.segments.append(Pattern.Segment(range(30, 38)))
    d = Pattern.Bond('T', 13, 33, True, 5)
    p.bonds.append(d)
    assert p.dist(7, 23) == 4
    assert p.dist(25, 34) == 5
    assert p.dist(25, 34, tlimit=2) == 5
    assert p.dist(25, 34, tlimit=1) == float('inf')
    p.bonds.remove(b)
    assert p.dist(7, 23) > 999

    e = Pattern.Bond('H', 25, 34, False, None)
    p.bonds.append(b)
    p.bonds.append(e)
    assert p.dist(25, 34, hlimit=1, tlimit=2) == 1
    assert p.dist(25, 34, hlimit=0, tlimit=2) == 5


def test_localise():
    p = Pattern.Pattern()
    s1 = Pattern.Segment(range(110, 120))
    s2 = Pattern.Segment(range(50, 60))
    s3 = Pattern.Segment(range(250, 270))
    p.segments = [s1, s2, s3]
    a = Pattern.Bond('H', 113, 260, False, None)
    b = Pattern.Bond('T', 54, 114, False, 3)
    p.bonds = [a, b]

    exp = Pattern.Pattern()
    s4 = Pattern.Segment(range(10))
    s5 = Pattern.Segment(range(100, 120))
    s6 = Pattern.Segment(range(200, 210))
    exp.segments = [s4, s5, s6]
    c = Pattern.Bond('H', 3, 110, False, None)
    d = Pattern.Bond('T', 204, 4, False, 3)
    exp.bonds = [c, d]

    assert p.localise(a) == exp


def test_handleresidue():
    p = Pattern.Pattern()
    r = 'IGKP'
    p.residue = r
    p.handleresidue(2)
    assert p.residue == 'LAEP'
    p.residue = r
    p.handleresidue(4)
    assert p.residue == 'LLEE'
    p.residue = '??CX'
    p.handleresidue(3)
    assert p.residue == 'XXAX'