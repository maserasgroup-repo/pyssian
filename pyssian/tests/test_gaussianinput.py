import unittest
from pathlib import Path
from pyssian import GaussianInFile

TEST_FILEDIR = Path(__file__).parent.resolve() / 'test_files'

class GaussianInFileTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.testfiles = (list(TEST_FILEDIR.glob('gaussian_input_*.com')) +
                         list(TEST_FILEDIR.glob('gaussian_input_*.gjf')) )
        cls.files_as_text = []
        for gaufile in cls.testfiles:
            with open(gaufile,'r') as F: 
                cls.files_as_text.append(F.read())

    def test_read_without_error(self):
        # Ensure that there is no error during the GaussianInput file 
        # reading. An error in this test indicates a wrong assumption 
        # about the valid structure of a gaussian input file. Note 
        # that if gaussian throws an error, it should throw an error.
        for i,gaufile in enumerate(self.testfiles):
            with self.subTest(testfile=i):
                with GaussianInFile(gaufile) as GIF: 
                    GIF.read()

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

        