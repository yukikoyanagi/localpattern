import os
import Protein


def test_fromfiles():
    prt = Protein.Protein('test')
    cd = os.path.dirname(__file__)
    prt.fromfiles(os.path.join(cd, 'data', 'tst001.txt'),
                  os.path.join(cd, 'data', 'tert001.txt'))
    assert len(prt.hbonds) == 6