import os
from .. import Protein
from .. import Pattern
from .. import Option
from .. import conv


def test_fromfiles():
    prt = Protein.Protein('test')
    cd = os.path.dirname(__file__)
    prt.fromfiles(os.path.join(cd, 'data', 'tst001.txt'),
                  os.path.join(cd, 'data', 'tert001.txt'))
    assert len(prt.hbonds) == 5
    assert prt.hbonds[2].donor == 1323
    assert prt.hbonds[2].acceptor == 1313
    er = [0.491790, 0.011741, 0.870635,
          -0.303610, 0.939465, 0.158829,
          -0.816066, -0.342444, 0.465584]
    ev = conv.mat2ev(conv.lst2mat(er))
    el = [f for f in ev[1]] + [ev[0]]
    assert max([i - j for (i, j) in zip(prt.hbonds[2].rotation, el)]) < 1**(-6)
    assert len(prt.tbonds) == 1
    assert prt.tbonds[0].left == 1267
    assert prt.tbonds[0].right == 1333
    er = [-0.3885338861613514,
          0.35503192946219037,
          -0.8502903906117677,
          0.24198545160937596,
          0.9297146738108265,
          0.2776214446152077,
          0.8890919302939176,
          -0.09789256541322623,
          -0.44713821702361]
    ev = conv.mat2ev(conv.lst2mat(er))
    el = [f for f in ev[1]] + [ev[0]]
    assert max([i - j for (i, j) in zip(prt.tbonds[0].rotation, el)]) < 1**(-6)


def test_istwisted():
    prt = Protein.Protein('test')
    cd = os.path.dirname(__file__)
    prt.fromfiles(os.path.join(cd, 'data', 'tst001.txt'),
                  os.path.join(cd, 'data', 'tert001.txt'))
    assert not prt.istwisted(prt.hbonds[2].rotation)
    assert prt.istwisted(prt.tbonds[0].rotation)


def test_getbondsat():
    prt = Protein.Protein('test')
    cd = os.path.dirname(__file__)
    prt.fromfiles(os.path.join(cd, 'data', 'tst001.txt'),
                  os.path.join(cd, 'data', 'tert001.txt'))
    assert len(prt.getbondsat(1323)) == 1
    assert isinstance(prt.getbondsat(1267)[0], Protein.Tbond)
    assert prt.getbondsat(1267)[0].right == 1333


def test_grow():
    h1 = Protein.Hbond(30, 61, 'ABCD', range(9), 1)
    t1 = Protein.Tbond(31, 46, 'ABCD', range(9), 3.14)
    prt = Protein.Protein('test')
    prt.hbonds = [h1]
    prt.tbonds = [t1]
    pat = Pattern.Pattern(segments=[[30], [61]],
                          bonds=[Pattern.Bond('H', 30, 61, False, None)],
                          residue='ABCD',
                          rotation=h1.rotation)
    expected = Pattern.Pattern(segments=[[29,30,31], [60,61,62]],
                               bonds=[Pattern.Bond('H', 30, 61, False, None)],
                               residue='ABCD',
                               rotation=h1.rotation)
    pat = prt.grow(pat)
    assert pat.bonds == expected.bonds
    assert pat.segments == expected.segments
    assert pat == expected
    expected.segments = [range(28,33), [46], range(59,64)]
    expected.bonds.append(Pattern.Bond('T', 31, 46, False, 3.14))
    assert prt.grow(pat) == expected

    prt.tbonds = []
    h2 = Protein.Hbond(31, 62, 'ABCD', range(9), 2)
    prt.hbonds.append(h2)
    pat = Pattern.Pattern(segments=[[30], [61]],
                          bonds=[Pattern.Bond('H', 30, 61, False, None)],
                          residue='ABCD',
                          rotation=h1.rotation)
    pat = prt.grow(pat)
    assert len(pat.bonds) == 2


def test_findpattern2():
    cd = os.path.dirname(__file__)
    opt = Option.Option(os.path.join(cd, 'data', 'step105_opts'))
    h1 = Protein.Hbond(15, 25, 'ABCD', range(9), 1)
    t1 = Protein.Tbond(24, 36, 'ABCD', range(9), 3.14)
    prt = Protein.Protein('test')
    prt.hbonds = [h1]
    prt.tbonds = [t1]
    center = Pattern.Bond('H', 15, 25, prt.istwisted(h1.rotation), None)
    pat = prt.findpattern2(h1, opt)
    cbond = Pattern.Bond('H', 3, 103, center.twisted, center.vdw)
    exp = Pattern.Pattern(segments=[range(7),
                                    range(100, 107),
                                    range(200, 203)],
                          bonds=[cbond,
                                 Pattern.Bond('T', 102, 201, False, 3.14)],
                          residue='LXLE',
                          rotation=h1.rotation,
                          length=-2)

    assert pat == exp

    h2 = Protein.Hbond(243, 251, 'LXLE', range(9), 2)
    t2 = Protein.Tbond(250, 340, 'ABCD', range(9), 0.45)
    prt2 = Protein.Protein('test2')
    prt2.hbonds = [h2]
    prt2.tbonds = [t2]
    pat2 = prt2.findpattern2(h2, opt)
    assert pat == pat2

    opt = Option.Option(os.path.join(cd, 'data', 'step106_opts'))
    h1 = Protein.Hbond(2, 34, 'GAX?', range(9), 1)
    t1 = Protein.Tbond(23, 33, 'ABCD', range(0, -9, -1), 3.1)
    prt.hbonds = [h1]
    prt.tbonds = [t1]
    pat = prt.findpattern2(h1, opt)
    exp = Pattern.Pattern(segments=[range(15),
                                    range(100, 124)],
                          bonds=[Pattern.Bond('H', 7, 116, False, None),
                                 Pattern.Bond('T', 105, 115, False, 3.1)],
                          residue='AAXX',
                          rotation=h1.rotation,
                          length='L')
    assert pat == exp

    opt = Option.Option(os.path.join(cd, 'data', 'step105_opts'))
    p = Protein.Protein('test3')
    hf = os.path.join(cd, 'data', 'tst001.txt')
    tf = os.path.join(cd, 'data', 'tert001.txt')
    p.fromfiles(hf, tf)
    hbnd = p.hbonds[-2]
    pat = p.findpattern2(hbnd, opt)
    assert len(pat.bonds) == 4
    tbonds = [b for b in pat.bonds if b.type == 'T']
    assert len(tbonds) == 1

    opt = Option.Option(os.path.join(cd, 'data', 'step148_opts'))
    p = Protein.Protein('test4')
    hf = os.path.join(cd, 'data', 'tst002.txt')
    tf = os.path.join(cd, 'data', 'tert002.txt')
    p.fromfiles(hf, tf)
    pat = p.findpattern2(p.hbonds[0], opt)
    assert isinstance(pat, Pattern.Pattern)

    opt = Option.Option(os.path.join(cd, 'data', 'step116_opts'))
    p = Protein.Protein('test5')
    hf = os.path.join(cd, 'data', 'tst003.txt')
    p.fromfiles(hf, None)
    assert str(p.findpattern2(p.hbonds[0], opt)) == ('0-6:100-106_'
                                                     'H3:103-'
                                                     'H6:R-'
                                                     'HR:2-'
                                                     'HR:5-_'
                                                     'XXXX_L')


def test_getogtbonds():
    p = Protein.Protein('testprot')
    p.hbonds = [Protein.Hbond(3, 13, 'ABCD', range(4), 1),
                Protein.Hbond(24, 34, 'ABCD', range(4), 2)]
    p.tbonds = [Protein.Tbond(11, 21, 'ABCD', range(9), 0.12),
                Protein.Tbond(14, 35, 'ABCD', range(9), 0.12)]
    pat = Pattern.Pattern('testpat')
    pat.segments = [Pattern.Segment(range(i, i+7)) for i in range(0, 21, 10)]
    pat.bonds = [Pattern.Bond('H', 3, 13, False, None),
                 Pattern.Bond('T', 11, 21, False, 0.12)]
    l = p.getogtbonds(pat)
    assert len(l) == 1
    assert l[0] == (Pattern.Bond('T', 14, 35, False, 0.13), 14, 35)


def test_addt():
    def startpat():
        sp = Pattern.Pattern('pat')
        sp.segments = [Pattern.Segment(range(10))]
        return sp

    p = Protein.Protein('test')
    r = 'ABCD'
    t = range(9)
    v = 0.12
    p.tbonds = [Protein.Tbond(3, 13, r, t, v),
                Protein.Tbond(12, 24, r, t, v),
                Protein.Tbond(14, 34, r, t, v),
                Protein.Tbond(15, 35, r, t, v),
                Protein.Tbond(33, 45, r, t, v)]
    pat = startpat()
    b = Pattern.Bond('T', 3, 13, False, v)
    newpat = p.addt(pat, b, 13, 0, 0)
    assert newpat.segments == [range(10), [13]]

    pat = startpat()
    newpat = p.addt(pat, b, 13, 1, 1)
    assert newpat.segments == [range(10), range(12, 15)]
    assert len(newpat.bonds) == 1

    pat = startpat()
    newpat = p.addt(pat, b, 13, 2, 1)
    assert newpat.segments == [range(10), range(11, 16), [24], [34]]
    assert len(newpat.bonds) == 3

    pat = startpat()
    newpat = p.addt(pat, b, 13, 7, 2)
    atoms = set([a for seg in newpat.segments for a in seg])
    assert atoms == set(range(40) + range(42, 49))
    assert len(set(newpat.bonds)) == 5



