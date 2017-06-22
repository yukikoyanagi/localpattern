import os
from .. import Protein


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
    assert max([i - j for (i, j) in zip(prt.hbonds[2].rotation, er)]) < 1**(-6)
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
    assert max([i - j for (i, j) in zip(prt.tbonds[0].rotation, er)]) < 1**(-6)
