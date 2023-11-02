import unittest
from pathlib import Path
from pyssian import GaussianInFile

TEST_FILEDIR = Path(__file__).parent.resolve() / 'test_files'

class GaussianInFileTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.testfiles = TEST_FILEDIR.glob('gaussian_input_*.{gjf,com}')
        cls.files_as_text = []
        for gaufile in cls.testfiles: 
            with open(gaufile,'r') as F: 
                cls.files_as_text.append(F.read())

    def test_read_str_equal(self):
        msg = 'Read file does not match its str representation'
        for gaufile,sol in zip(self.testfiles,self.files_as_text):
            with GaussianInFile(gaufile) as GIF: 
                GIF.read()
            test = str(GIF)
            self.assertEqual(test,sol,msg)

    def test_nprocshared(self): 
        GIF = GaussianInFile()
        GIF.nprocs = 24
        msg = 'incorrect nprocshared keyword creation through property nprocs'
        self.assertTrue('nprocshared' in GIF.preprocessing,msg)
        self.assertEqual(24, GIF.preprocessing['nprocshared'],msg)
        GIF.nprocs = 12
        msg = 'incorrect nprocshared keyword modification through property nprocs'
        self.assertTrue('nprocshared' in GIF.preprocessing,msg)
        self.assertEqual(12, GIF.preprocessing['nprocshared'],msg)

        