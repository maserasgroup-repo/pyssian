from pyssian.chemistryutils import is_basis,is_method
import unittest

class ChemistryUtilsTest(unittest.TestCase):
    def setUp(self):
        self.Valid_Basis = '6-311+g(d,p) 6-31g* cc-pVTZ D95V* LanL2DZ SDD28 Def2SVP UGBS2P2++'.split()
        self.Fake_Basis = '6-311g+(d,p) 6-31*g ccpVTZ D96V* LanL2TZ SDD Def2SP UGBS2++P2'.split()
        self.Valid_Methods = 'ub3lyp mp2 casscf ccsd(t) rm062x WB97XD pbepbe'.split()
        self.Fake_Methods = 'bu3lyp pm2 bw97xd m06-2x pbepbe0'.split()
        self.Usual_Keywords = 'opt freq scrf scf #p calcfc empiricaldispersion'.split()
    def test_valid_isbasis(self):
        msg = 'Valid basis not properly recognized'
        for valid in self.Valid_Basis:
            self.assertTrue(is_basis(valid),msg)
    def test_fake_isbasis(self):
        msg = 'Fake basis recognized as valid'
        for fake in self.Fake_Basis:
            self.assertFalse(is_basis(fake),msg)
    def test_valid_ismethod(self):
        msg = 'Valid method not properly recognized'
        for valid in self.Valid_Methods:
            self.assertTrue(is_method(valid),msg)
    def test_fake_ismethod(self):
        msg = 'Fake method recognized as valid'
        for fake in self.Fake_Methods:
            self.assertFalse(is_method(fake),msg)
    def test_usual_keywords(self):
        msg1 = 'Keyword recognized as basis'
        msg2 = 'Keyword recognized as method'
        for keyword in self.Usual_Keywords:
            self.assertFalse(is_basis(keyword),msg1)
            self.assertFalse(is_method(keyword),msg2)
