from .. import Pattern


def test_equality():
    # test equality of Segment objects
    A = Pattern.Segment()
    B = Pattern.Segment()
    C = Pattern.Segment()
    for i in range(6):
        A.append(Pattern.Atom(i, chr(69+i)))
        B.append(Pattern.Atom(i+2, chr(69+i)))
        C.append(Pattern.Atom(i, chr(75+i)))
    assert A == B
    assert A != C

    # test equality of Pattern objects
    P1 = Pattern.Pattern()
    P2 = Pattern.Pattern()
    P3 = Pattern.Pattern()
