from pyssian.chemistryutils import is_basis,is_method
import unittest

class ChemistryUtilsTest(unittest.TestCase):
    def setUp(self):
        self.valid_basis = '6-311+g(d,p) 6-31g* cc-pVTZ D95V* LanL2DZ SDD28 Def2SVP UGBS2P2++'.split()
        self.fake_basis = '6-311g+(d,p) 6-31*g ccpVTZ D96V* LanL2TZ SDD Def2SP UGBS2++P2'.split()
        self.valid_methods = 'ub3lyp mp2 casscf ccsd(t) rm062x WB97XD pbepbe'.split()
        self.fake_methods = 'bu3lyp pm2 bw97xd m06-2x pbepbe0'.split()
        self.usual_keywords = 'opt freq scrf scf #p calcfc empiricaldispersion nmr'.split()
    def test_valid_isbasis(self):
        msg = 'Valid basis not properly recognized'
        for valid in self.valid_basis:
            self.assertTrue(is_basis(valid),msg)
    def test_fake_isbasis(self):
        msg = 'Fake basis recognized as valid'
        for fake in self.fake_basis:
            self.assertFalse(is_basis(fake),msg)
    def test_valid_ismethod(self):
        msg = 'Valid method not properly recognized'
        for valid in self.valid_methods:
            self.assertTrue(is_method(valid),msg)
    def test_fake_ismethod(self):
        msg = 'Fake method recognized as valid'
        for fake in self.fake_methods:
            self.assertFalse(is_method(fake),msg)
    def test_usual_keywords(self):
        msg1 = 'Keyword={} recognized as basis'
        msg2 = 'Keyword={} recognized as method'
        for keyword in self.usual_keywords:
            with self.subTest(keyword=keyword):
                self.assertFalse(is_basis(keyword),msg1.format(keyword))
                self.assertFalse(is_method(keyword),msg2.format(keyword))
