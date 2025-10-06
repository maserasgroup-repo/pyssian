import unittest
from pathlib import Path
from itertools import chain, zip_longest

from pyssian.linkjobparsers import *

TEST_FILEDIR = Path(__file__).parent.resolve() / 'test_files'

SMARK = '\n## Split Here ##\n'
MMARK = '\n## Match ##\n'

class TestGeneralLinkJob(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.testfile = TEST_FILEDIR.joinpath('General.txt')
        cls.refile = TEST_FILEDIR.joinpath('General_regex.txt')
        # Read and store the Examples
        with open(cls.testfile,'r') as F:
            txt = F.read()
        cls.objects = [GeneralLinkJob(i) for i in txt.split(SMARK)]

    def test_regex(self):
        msg = 're_enter regex does not match properly'
        regex = GeneralLinkJob.re_enter
        with open(self.refile,'r') as F:
            txt = F.read()
        items = [i.strip() for i in txt.split(SMARK)]
        for obj,solution in zip(self.objects,items):
            txt = obj.text
            test = regex.findall(txt)[0]
            self.assertTrue(test == solution,msg)

    def test_init(self):
        msg = 'Incorrect parsing of GeneralLinkJob'
        for obj in self.objects:
            self.assertTrue(bool(obj.text),msg)

    def test_init_empty(self):
        msg = 'Incorrect empty initialization of GeneralLinkJob'
        obj = GeneralLinkJob('Some Text',asEmpty=True)
        self.assertFalse(bool(obj.text),msg)
        self.assertTrue(obj.number == -1,msg)

class TestLink1(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.testfile = TEST_FILEDIR.joinpath('l1.txt')
        cls.refile = TEST_FILEDIR.joinpath('l1_regex.txt')
        cls.propfile = TEST_FILEDIR.joinpath('l1_properties.txt')
        # Read and store the Examples
        with open(cls.testfile,'r') as F:
            txt = F.read()
        cls.objects = [Link1(i) for i in txt.split(SMARK)]

    def test_init(self):
        msg = 'Incorrect parsing of Link1'
        for obj in self.objects:
            self.assertTrue(bool(obj.text),msg)
            self.assertEqual(obj.number,1)
            self.assertTrue(bool(obj.info),msg)

    def test_regex_link0(self):
        msg = 're_link0 regex does not match properly'
        regex = Link1.re_link0
        with open(self.refile) as F:
            txt = F.read()
        items = [i for i in txt.split(SMARK)]
        solutions = [['%nprocshared=8','%mem=18000MB'],
                     []]
        for txt,solution in zip_longest(items,solutions):
            test = regex.findall(txt)
            self.assertTrue(test == solution,msg)

    def test_regex_commandline(self):
        msg = 're_commandline regex does not match properly'
        regex = Link1.re_commandline
        with open(self.refile) as F:
            txt = F.read()
        items = [i for i in txt.split(SMARK)]
        solutions = [['#p opt=(calcfc,ts,noeigentest) freq b3lyp/6-31+g(d) scrf=(solvent=dich\nloromethane,smd) nosymm empiricaldispersion=gd3\n'],
                     ['#P Geom=AllCheck Guess=TCheck SCRF=Check GenChk RB3LYP/6-31+G(d) Freq\n']]
        for txt,solution in zip_longest(items,solutions):
            test = regex.findall(txt)
            self.assertTrue(test == solution,msg)

    def test_regex_IOps(self):
        msg = 're_IOps regex does not match properly'
        regex = Link1.re_IOps
        with open(self.refile) as F:
            txt = F.read()
        items = [i for i in txt.split(SMARK)]
        solutions = [['1/5=1,10=4,11=1,14=-1,18=20,26=3,38=1/1,3;',
                      '2/9=110,12=2,15=1,17=6,18=5,40=1/2;',
                      '4//1;',
                      '7/10=1,18=20,25=1,30=1/1,2,3,16;',
                      '1/5=1,11=1,14=-1,18=20,26=3/3(-5);'],
                     ['4/5=101/1;',
                      '5/5=2,53=9,98=1/2;',
                      '8/6=4,10=90,11=11/1;',
                      '99//99;']]
        for txt,solution in zip_longest(items,solutions):
            test = regex.findall(txt)
            self.assertTrue(test == solution,msg)

    def test_regex_internaljob(self):
        msg = 're_internaljob regex does not match properly'
        regex = Link1.re_internaljob
        with open(self.refile) as F:
            txt = F.read()
        items = [i for i in txt.split(SMARK)]
        solutions = [[],['2',]]
        for txt,solution in zip_longest(items,solutions):
            test = regex.findall(txt)
            self.assertTrue(test == solution,msg)

    def test_init_empty(self):
        msg = 'Incorrect empty initialization of Link1'
        text = 'Link1:  Proceeding to internal job step number 13.'
        obj = Link1(text,asEmpty=True)
        self.assertFalse(bool(obj.text),msg)
        self.assertEqual(obj.number,1)
        self.assertTrue(bool(obj.info),msg)

    def test_guess_type(self):
        msg = 'Incorrect type guess'
        text = 'Link1:  Proceeding to internal job step number 13.'
        obj = Link1(text,asEmpty=True)
        types = {0:'Constrained Optimization',
                 1:'Optimization',
                 2:'Frequency Calculation',
                 3:'Unidentified'}
        keywords = ['modredundant','opt','modredundant',
                    'opt','freq','opt','freq',
                    'optimization','#p','redundant']
        solutions = [0,0,0,1,1,1,2,3,3,3]
        while keywords:
            keywords.pop(0)
            obj.commandline = ' '.join(keywords)
            sol = types[solutions.pop(0)]
            test = obj._guess_type()
        self.assertTrue(test == sol,msg)

    def test_commandline(self):
        msg = 'commandline does not maintain its integrity'
        with open(self.propfile) as F:
            txt = F.read()
        items = [i.rstrip('\n').split('\n') for i in txt.split(SMARK)]
        for item,obj in zip_longest(items,self.objects):
            _ , sol = item
            self.assertTrue(obj.commandline == sol,msg)

    def test_info(self):
        msg = 'InternalJobInfo not properly added'
        InternalJobInfo = Link1.InternalJobInfo
        with open(self.propfile) as F:
            txt = F.read()
        items = [i.rstrip('\n').split('\n') for i in txt.split(SMARK)]
        for item,obj in zip_longest(items,self.objects):
            Aux = item[0].split(',')
            sol = InternalJobInfo(int(Aux[0]),Aux[1],Aux[2]=='1')
            self.assertTrue(obj.info == sol,msg)

class TestLink101(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.testfile = TEST_FILEDIR.joinpath('l101.txt')
        cls.charge_spin_solutions = [(-1,2),
                                     (0,2),
                                     (2,2),
                                     (2,2)]
        # Read and store the Examples
        with open(cls.testfile,'r') as F:
            txt = F.read()
        cls.objects = [Link101(i) for i in txt.split(SMARK)]

    def test_init(self):
        msg = 'Incorrect parsing of Link101'
        for obj in self.objects:
            self.assertTrue(bool(obj.text),msg)
            self.assertEqual(obj.number,101)
            self.assertTrue(obj.charge is not None,msg)
            self.assertTrue(obj.spin is not None,msg)

    def test_regex_charge(self):
        msg = 're_charge regex does not match properly'
        regex = Link101.re_charge
        solutions = [str(i) for i,_ in self.charge_spin_solutions]
        for obj,solution in zip_longest(self.objects,solutions):
            test = regex.findall(obj.text)[0]
            self.assertTrue(test == solution,msg)

    def test_regex_spin(self):
        msg = 're_spin regex does not match properly'
        regex = Link101.re_spin
        solutions = [str(i) for _,i in self.charge_spin_solutions]
        for obj,solution in zip_longest(self.objects,solutions):
            test = regex.findall(obj.text)[0]
            self.assertTrue(test == solution,msg)

    def test_init_empty(self):
        msg = 'Incorrect empty initialization of Link101'
        obj = Link101('-> Input Text',asEmpty=True)
        self.assertFalse(bool(obj.text),msg)
        self.assertEqual(obj.number,101)

    def test_charge_spin(self):
        msg = '{} is not properly parsed'.format
        for item,obj in zip_longest(self.charge_spin_solutions,self.objects):
            charge, spin = item
            self.assertTrue(obj.charge == charge ,msg('charge'))
            self.assertTrue(obj.spin == spin, msg('spin'))

class TestLink103(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.testfile = TEST_FILEDIR.joinpath('l103.txt')
        cls.parametersfile = TEST_FILEDIR.joinpath('l103_parameters.txt')
        cls.derivativesfile = TEST_FILEDIR.joinpath('l103_derivatives.txt')
        cls.convergencefile = TEST_FILEDIR.joinpath('l103_convergence.txt')
        # Read and store the Examples
        with open(cls.testfile,'r') as F:
            txt = F.read()
        cls.objects = [Link103(i) for i in txt.split(SMARK)]

    def test_init(self):
        msg = 'Incorrect parsing of Link103'
        for obj in self.objects:
            self.assertTrue(bool(obj.text),msg)
            self.assertEqual(obj.number,103,msg)
            self.assertTrue(bool(obj.mode),msg)

    def test_init_empty(self):
        msg = 'Incorrect empty initialization of Link101'
        obj = Link103('-> Input Text',asEmpty=True)
        attrs = ['mode','state','convergence', 'parameters', 'derivatives',
                 'stepnumber', 'scanpoint']
        for attr in attrs:
            test = getattr(obj,attr)
            self.assertTrue(not bool(test) or test is None)

    def test_locate_mode(self):
        msg = 'Mode {} not properly recognized'.format
        solutions = ['Init','Iteration','End','End']
        for i,(obj,solution) in enumerate(zip(self.objects,solutions)):
            with self.subTest(Test_Object=i,mode=solution):
                self.assertTrue(obj.mode == solution,msg(solution))

    def test_regex_parameters(self):
        msg = 're_parameters regex does not match properly'
        regex = Link103.re_parameters
        with open(self.parametersfile,'r') as F:
            solutions = [sample.split('\n') for sample in F.read().split(SMARK)]
        for i,(obj,solution) in enumerate(zip_longest(self.objects,solutions)):
            test = regex.findall(obj.text)
            if not solution:
                with self.subTest(object=i,parameter='None'):
                    self.assertFalse(bool(test))
            for j,(t,s) in enumerate(zip_longest(test,solution,fillvalue='')):
                with self.subTest(object=i,parameter=j):
                    self.assertEqual(t,s,msg)

    def test_regex_derivatives(self):
        msg = 're_derivatives regex does not match properly'
        regex = Link103.re_derivatives
        solutions = []
        with open(self.derivativesfile,'r') as F:
            txt = F.read()
        for sample in txt.split(SMARK):
            Aux = []
            if sample:
                for s in sample.split('\n'):
                    Aux.append(tuple(s.split(',')))
                solutions.append(Aux)
            else:
                solutions.append([])
        for i,(obj,solution) in enumerate(zip_longest(self.objects,solutions)):
            test = regex.findall(obj.text)
            if not solution:
                with self.subTest(object=i,parameter='None'):
                    self.assertFalse(bool(test))
            for j,(t,s) in enumerate(zip_longest(test,solution,fillvalue='')):
                with self.subTest(object=i,parameter=j):
                    self.assertEqual(t,s,msg)

    def test_regex_convergence(self):
        msg = 're_convergence regex does not match properly'
        regex = Link103.re_convergence
        solutions = []
        with open(self.convergencefile,'r') as F:
            txt = F.read()
        for sample in txt.split(SMARK):
            Aux = []
            if sample:
                for s in sample.split('\n'):
                    Aux.append(tuple(s.split(',')))
                solutions.append(Aux)
            else:
                solutions.append([])
        for i,(obj,solution) in enumerate(zip_longest(self.objects,solutions)):
            test = regex.findall(obj.text)
            if not solution:
                with self.subTest(object=i,parameter='None'):
                    self.assertFalse(bool(test))
            for j,(t,s) in enumerate(zip_longest(test,solution,fillvalue='')):
                print(f'subtest {i} parameter={j}')
                print(t,s)
                with self.subTest(object=i,parameter=j):
                    self.assertEqual(t,s,msg)

    def test_locate_numbers(self):
        msg = 'stepnumber {} not properly recognized'.format
        solutions = [0,1,420,66,31]
        for i,(obj,solution) in enumerate(zip_longest(self.objects,solutions)):
            with self.subTest(Test_Object=i,stepnumber=solution):
                self.assertEqual(obj.stepnumber,solution,msg(solution))
        msg = 'scanpoint {} not properly recognized'.format
        solutions = [None,None,None,11,None]
        for i,(obj,solution) in enumerate(zip_longest(self.objects,solutions)):
            with self.subTest(Test_Object=i,scanpoint=solution):
                self.assertEqual(obj.scanpoint,solution,msg(solution))

class TestLink120(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.testfile = TEST_FILEDIR.joinpath('l120.txt')
        cls.re_energy_solutions = ['',
                                   '',
                                   '',
                                   '-2572.457446471110']
        cls.re_energy_partition_solutions = [[],
                                             [],
                                             [],
                                             [['1', 'low','model','0.403686133797'],
                                              ['2','high','model','-2573.050183259390'],
                                              ['3', 'low', 'real','0.996422922077']]]
        # Read and store the Examples
        with open(cls.testfile,'r') as F:
            txt = F.read()
        cls.objects = [Link120(i) for i in txt.split(SMARK)]

    def test_init(self):
        msg = 'Incorrect parsing of Link120'
        for obj in self.objects:
            self.assertTrue(bool(obj.text),msg)
            self.assertEqual(obj.number,120,msg)
            self.assertTrue(obj.energy is None or bool(obj.energy),msg)

    def test_regex_energy(self):
        msg = 're_energy regex does not match properly'
        regex = Link120.re_energy
        solutions = self.re_energy_solutions
        for obj,solution in zip_longest(self.objects,solutions):
            match = regex.findall(obj.text)
            self.assertEqual(bool(match),bool(solution),msg)
            if match: 
                self.assertEqual(match[0],solution,msg)
    
    def test_regex_energy_partitions(self):
        msg = 're_energy regex does not match properly'
        regex = Link120.re_energy_partitions
        solutions = self.re_energy_partition_solutions
        for obj,solution in zip_longest(self.objects,solutions):
            match = regex.findall(obj.text)
            self.assertEqual(bool(match),bool(solution),msg)
            if match:
                for (p,lev,model,energy),sol in zip(match,solution):
                    self.assertEqual([p,lev,model,energy],sol) 
            
    def test_init_empty(self):
        msg = 'Incorrect empty initialization of Link120'
        obj = Link120('Some Text',asEmpty=True)
        self.assertFalse(bool(obj.text),msg)
        self.assertEqual(obj.number,120)

    def test_energy(self):
        msg = 'Energy value not properly read'
        solutions = self.re_energy_solutions
        for obj,solution in zip_longest(self.objects,solutions):
            test = obj.energy
            self.assertEqual(bool(test),bool(solution),msg)
            if solution: 
                self.assertEqual(test,float(solution),msg)
    
    def test_energy_partition(self):
        msg = 'Energy partitions not properly read'
        solutions = self.re_energy_partition_solutions
        EnergyPartition = Link120._EnergyPartition
        for obj,solution in zip_longest(self.objects,solutions):
            match = obj.energy_partitions
            self.assertEqual(bool(match),bool(solution),msg)
            if not match:
                continue
            for (p,level,model,energy),sol in zip_longest(match,solution):
                p,l,m,e = sol
                sol = EnergyPartition(int(p),l,m,float(e))
                test = EnergyPartition(int(p),level,model,float(energy))
                self.assertEqual(test,sol,msg)
    
class TestLink122(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.testfile = TEST_FILEDIR.joinpath('l122.txt')
        
        cls.re_complex_solutions = [
                                 '-622.598617397384',
                                 '-554.437959704946'
        ]
        cls.re_bsse_solutions = [
                                 '0.000785012282',
                                 '0.023214737127'
        ]
        cls.re_fragments_solutions = [
                                 '-622.391294173148',
                                 '-553.667179856454'
        ]
        # Read and store the Examples
        with open(cls.testfile,'r') as F:
            txt = F.read()
        cls.objects = [Link122(i) for i in txt.split(SMARK)]

    def test_init(self):
        msg = 'Incorrect parsing of Link122'
        for obj in self.objects:
            self.assertTrue(bool(obj.text),msg)
            self.assertEqual(obj.number,122,msg)
            self.assertTrue(obj.energy_complex is None or bool(obj.energy_complex),msg)
            self.assertTrue(obj.total_energy_fragments is None or bool(obj.total_energy_fragments),msg)
            self.assertTrue(obj.bsse_correction is None or bool(obj.bsse_correction),msg)

    def test_regex_complex(self):
        msg = 're_energy regex does not match properly'
        regex = Link122.re_complex_energy
        solutions = self.re_complex_solutions
        for obj,solution in zip_longest(self.objects,solutions):
            match = regex.findall(obj.text)
            self.assertEqual(bool(match),bool(solution),msg)
            if match: 
                self.assertEqual(match[0],solution,msg)

    def test_regex_bsse(self):
        msg = 're_energy regex does not match properly'
        regex = Link122.re_bsse
        solutions = self.re_bsse_solutions
        for obj,solution in zip_longest(self.objects,solutions):
            match = regex.findall(obj.text)
            self.assertEqual(bool(match),bool(solution),msg)
            if match: 
                self.assertEqual(match[0],solution,msg)

    def test_regex_fragments(self):
        msg = 're_energy regex does not match properly'
        regex = Link122.re_fragments_energy
        solutions = self.re_fragments_solutions
        for obj,solution in zip_longest(self.objects,solutions):
            match = regex.findall(obj.text)
            self.assertEqual(bool(match),bool(solution),msg)
            if match: 
                self.assertEqual(match[0],solution,msg)

    def test_init_empty(self):
        msg = 'Incorrect empty initialization of Link122'
        obj = Link122('Some Text',asEmpty=True)
        self.assertFalse(bool(obj.text),msg)
        self.assertEqual(obj.number,122)

    def test_energy_complex(self):
        msg = 'Energy value not properly read'
        solutions = self.re_complex_solutions
        for obj,solution in zip_longest(self.objects,solutions):
            test = obj.energy_complex
            self.assertEqual(bool(test),bool(solution),msg)
            if solution: 
                self.assertEqual(test,float(solution),msg)

    def test_total_energy_fragments(self):
        msg = 'Energy value not properly read'
        solutions = self.re_fragments_solutions
        for obj,solution in zip_longest(self.objects,solutions):
            test = obj.total_energy_fragments
            self.assertEqual(bool(test),bool(solution),msg)
            if solution: 
                self.assertEqual(test,float(solution),msg)

    def test_bsse_correction(self):
        msg = 'Energy value not properly read'
        solutions = self.re_bsse_solutions
        for obj,solution in zip_longest(self.objects,solutions):
            test = obj.bsse_correction
            self.assertEqual(bool(test),bool(solution),msg)
            if solution: 
                self.assertEqual(test,float(solution),msg)


class TestLink123(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.testfile = TEST_FILEDIR.joinpath('l123.txt')
        cls.orientationfile = TEST_FILEDIR.joinpath('l123_orientation.txt')
        # Read and store the Examples
        with open(cls.testfile,'r') as F:
            txt = F.read()
        cls.objects = [Link123(i) for i in txt.split(SMARK)]

    def test_init(self):
        msg = 'Incorrect parsing of Link123'
        for obj in self.objects:
            self.assertTrue(bool(obj.text),msg)
            self.assertEqual(obj.number,123,msg)

    def test_init_empty(self):
        msg = 'Incorrect empty initialization of Link101'
        obj = Link123('-> Input Text',asEmpty=True)
        attrs = ['orientation', 'direction', 'step', 'reactioncoord']
        for attr in attrs:
            test = getattr(obj,attr)
            self.assertTrue(not bool(test) or test is None)

    def test_regex_orientation(self):
        msg = 're_orientation regex does not match properly'
        regex = Link123.re_orientation
        solutions = []
        with open(self.orientationfile,'r') as F:
            txt = F.read()
        for sample in txt.split(SMARK):
            Aux = []
            if sample:
                for s in sample.split('\n'):
                    Aux.append(tuple(s.split('\t')))
                solutions.append(Aux)
            else:
                solutions.append([])
        for i,(obj,solution) in enumerate(zip_longest(self.objects,solutions)):
            test = regex.findall(obj.text)
            if not solution:
                with self.subTest(object=i,atom='None'):
                    self.assertFalse(bool(test))
            for j,(t,s) in enumerate(zip_longest(test,solution,fillvalue='')):
                with self.subTest(object=i,atom=j):
                    self.assertEqual(t,s,msg)

    def test_locate_irc_step(self):
        msg0 = 'reactioncoord {} not properly recognized'.format
        msg1 = 'direction {} not properly recognized'.format
        msg2 = 'step {} not properly recognized'.format
        solutions = [(0,'',0),
                     (0,'FORWARD',1),
                     (0.04224,'FORWARD',2),
                     (9.58567,'FORWARD',228)]
        for i,(obj,solution) in enumerate(zip_longest(self.objects,solutions)):
            with self.subTest(Test_Object=i,reactioncoord=solution[0]):
                test = obj.reactioncoord
                sol = solution[0]
                self.assertEqual(test,sol ,msg0(sol))
            with self.subTest(Test_Object=i,reactioncoord=solution[1]):
                test = obj.direction
                sol = solution[1]
                self.assertEqual(test,sol ,msg1(sol))
            with self.subTest(Test_Object=i,reactioncoord=solution[2]):
                test = obj.step
                sol = solution[2]
                self.assertEqual(test,sol ,msg2(sol))

class TestLink202(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.testfile = TEST_FILEDIR.joinpath('l202.txt')
        cls.refile = TEST_FILEDIR.joinpath('l202_regex.txt')
        cls.orifile = TEST_FILEDIR.joinpath('l202_orientations.txt')
        cls.dmatfile = TEST_FILEDIR.joinpath('l202_dmatrices.txt')
        # Read and store the Examples
        with open(cls.testfile,'r') as F:
            txt = F.read()
        cls.objects = [Link202(i) for i in txt.split(SMARK)]
        # Read and store the Solutions to the orientation tests
        with open(cls.orifile,'r') as F:
            txt = F.read()
        cls.orientations = txt.split(SMARK)

    def test_init(self):
        msg = 'Incorrect parsing of Link202'
        for obj in self.objects:
            self.assertTrue(bool(obj.text),msg)
            self.assertEqual(obj.number,202)
            self.assertTrue(bool(obj.orientation),msg)

    def test_regex(self):
        msg = 're_orientation regex does not match properly'
        regex = Link202.re_orientation
        with open(self.refile) as F:
            txt = F.read()
        items = [i.split(MMARK) for i in txt.split(SMARK)]
        for i,(obj,solutions) in enumerate(zip(self.objects,items)):
            text = obj.text
            tests = regex.findall(text)
            for j,(test,sol) in enumerate(zip(tests,solutions)):
                with self.subTest(text=i,match=j):
                    self.assertEqual(test,sol,msg)

    def test_init_empty(self):
        msg = 'Incorrect empty initialization of Link202'
        obj = Link202('Some Text',asEmpty=True)
        self.assertFalse(bool(obj.text),msg)
        self.assertEqual(obj.number,202)
        self.assertFalse(bool(obj.orientation),msg)

    def test_orientation(self):
        msg = 'Orientation not properly parsed'
        AtomCoords = Link202._AtomCoords
        with open(self.orifile) as F:
            txt = F.read()
        items = []
        for sol in txt.split(SMARK):
            lines = sol.split('\n')
            mat = []
            for line in lines:
                if line.strip():
                    values = line.strip().split()
                    a, b, c = int(values[0]),int(values[1]),int(values[2])
                    d, e, f = float(values[3]),float(values[4]),float(values[5])
                    mat.append(AtomCoords(a,b,c,d,e,f))
            items.append(mat)
        for i,(obj,solution) in enumerate(zip(self.objects,items)):
            test = obj.orientation
            with self.subTest(text=i):
                self.assertEqual(test,solution,msg)

    def test_DistanceMatrix(self):
        msg = 'DistanceMatrix not properly parsed'
        obj = self.objects[3]
        dmatrix_sol = []
        with open(self.dmatfile) as F:
            for line in F:
                if line.strip():
                    Aux = line.strip().split()
                    i0 = int(Aux.pop(0))
                    i1 = Aux.pop(0)
                    i2 = list(map(float, Aux))
                    dmatrix_sol.append([i0,i1,i2])
        self.assertTrue(dmatrix_sol == obj.DistanceMatrix,msg)

    def test_AtNum2Sym_error(self):
        msg = 'Obtained an unobtainable Mapping'
        Tests = self.objects
        with self.assertRaises(StopIteration,msg=msg):
            Tests[0].get_atom_mapping()
        with self.assertRaises(StopIteration,msg=msg):
            Tests[1].get_atom_mapping()
        with self.assertRaises(StopIteration,msg=msg):
            Tests[2].get_atom_mapping()

    def test_AtNum2Sym(self):
        msg = 'Incorrect Mapping Obtained'
        SolDict = {1:'H',8:'O',6:'C',7:'N'}
        TestDict = self.objects[3].get_atom_mapping()
        self.assertTrue(TestDict == SolDict,msg)

class TestLink502(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.testfile = TEST_FILEDIR.joinpath('l502.txt')
        cls.re_energies_solutions = ['-2847.29438388','-2847.29924823',
                                     '-2847.29888359','-2847.29999372',
                                     '-2847.30039709','-2847.30030072',
                                     '-2847.30072337','-2847.30086336',
                                     '-2847.30091754','-2847.30086384',
                                     '-2847.30093740']
        cls.re_spin_solutions = [('0.7531', '0.7500'),
                                 ('0.7530', '0.7500'),
                                 ('0.7529', '0.7500'),
                                 ('0.7529', '0.7500'),
                                 ('0.7529', '0.7500'),
                                 ('0.7529', '0.7500'),
                                 ('0.7529', '0.7500'),
                                 ('0.7529', '0.7500'),
                                 ('0.7529', '0.7500'),
                                 ('0.7529', '0.7500'),
                                 ('0.7529', '0.7500')]
        
        # Read and store the Examples
        with open(cls.testfile,'r') as F:
            txt = F.read()
        cls.objects = [Link502(i) for i in txt.split(SMARK)]

    def test_init(self):
        msg = 'Incorrect parsing of Link502'
        for obj in self.objects:
            self.assertTrue(bool(obj.text),msg)
            self.assertEqual(obj.number,502,msg)
            self.assertTrue(obj.energy is not None,msg)

    def test_regex_energy(self):
        msg = 're_EDone regex does not match properly'
        regex = Link502.re_EDone
        solutions = self.re_energies_solutions
        for obj,solution in zip(self.objects,solutions):
            test = regex.findall(obj.text)[0]
            self.assertTrue(test == solution,msg)

    def test_regex_spin(self):
        msg = 're_spin regex does not match properly'
        regex = Link502.re_spin
        solutions = self.re_spin_solutions
        for obj,solution in zip(self.objects,solutions):
            test = regex.findall(obj.text)[0]
            self.assertTrue(test == solution,msg)

    def test_init_empty(self):
        msg = 'Incorrect empty initialization of Link502'
        obj = Link502('Some Text',asEmpty=True)
        self.assertFalse(bool(obj.text),msg)
        self.assertEqual(obj.number,502)

    def test_energy(self):
        msg = 'Energy value not properly read'
        solutions = [float(i) for i in self.re_energies_solutions]
        for obj,solution in zip(self.objects,solutions):
            test = obj.energy
            self.assertTrue(test == solution,msg)
    def test_spin(self):
        msg = 'S**2 values not properly read'
        solutions = [(float(i[0]),float(i[1])) for i in self.re_spin_solutions]
        for obj,solution in zip(self.objects,solutions):
            test = obj.spin
            self.assertTrue(test == solution,msg)

class TestLink601(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.testfile = TEST_FILEDIR.joinpath('l601.txt')
        cls.regexfile = TEST_FILEDIR.joinpath('l601_regex.txt')
        # Read and store the Examples
        with open(cls.testfile,'r') as F:
            txt = F.read()
        cls.objects = [Link601(i) for i in txt.split(SMARK)]

    def test_init(self):
        msg = 'Incorrect parsing of Link601'
        for obj in self.objects:
            self.assertTrue(bool(obj.text),msg)
            self.assertEqual(obj.number,601,msg)
            self.assertTrue(obj.mulliken_heavy,msg)
            self.assertTrue(obj.mulliken,msg)

    def test_regex_MullikenAtoms(self):
        msg = 're_MullikenAtoms does not match the first group\n{}\n!=\n{}'
        msg = msg.format
        msg2 = 're_MullikenAtoms does not match the second group\n{}\n!=\n{}'
        msg2 = msg2.format
        regex = Link601.re_MullikenAtoms
        with open(self.regexfile,'r') as F:
            txt = F.read()
        samples = txt.split(SMARK)
        solutions = []
        for sample in samples:
            Aux = sample.split('\n## re_MullikenHeavy ##\n')[0]
            match0,match1 = Aux.split('\n## Match1 ##\n')
            match0 = match0.split('## Match0 ##\n')[-1]
            solutions.append((match0,match1))
        for obj,solution in zip(self.objects,solutions):
            test = regex.findall(obj.text)[0]
            self.assertTrue(test[0] == solution[0],msg(test[0],solution[0]))
            self.assertTrue(test[1] == solution[1],msg2(test[1],solution[1]))
    def test_regex_MullikenHeavy(self):
        msg = 're_MullikenHeavy does not match the first group\n{}\n!=\n{}'
        msg = msg.format
        msg2 = 're_MullikenHeavy does not match the second group\n{}\n!=\n{}'
        msg2 = msg2.format
        regex = Link601.re_MullikenHeavy
        with open(self.regexfile,'r') as F:
            txt = F.read()
        samples = txt.split(SMARK)
        solutions = []
        for sample in samples:
            Aux = sample.split('\n## re_MullikenHeavy ##\n')[1]
            match0,match1 = Aux.split('\n## Match1 ##\n')
            match0 = match0.split('## Match0 ##\n')[-1]
            solutions.append((match0,match1))
        for obj,solution in zip(self.objects,solutions):
            test = regex.findall(obj.text)[0]
            self.assertTrue(test[0] == solution[0],msg(test[0],solution[0]))
            self.assertTrue(test[1] == solution[1],msg2(test[1],solution[1]))

    def test_init_empty(self):
        msg = 'Incorrect empty initialization of Link502'
        obj = Link601('Some Text',asEmpty=True)
        self.assertFalse(bool(obj.text),msg)
        self.assertEqual(obj.number,601)

class TestLink716(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.testfile        = TEST_FILEDIR.joinpath('l716.txt')
        cls.dipolefile      = TEST_FILEDIR.joinpath('l716_redipole.txt')
        cls.Thermofile      = TEST_FILEDIR.joinpath('l716_reThermo.txt')
        cls.EContribfile    = TEST_FILEDIR.joinpath('l716_reEContrib.txt')
        cls.IRSpectrumfile  = TEST_FILEDIR.joinpath('l716_reIRSpectrum.txt')
        cls.Frequenciesfile = TEST_FILEDIR.joinpath('l716_reFrequencies.txt')
        cls.freqtextfile    = TEST_FILEDIR.joinpath('l716_refreqtxt.txt')
        cls.freqdispfile    = TEST_FILEDIR.joinpath('l716_refreqdisplacements.txt')
        # Read and store the Examples
        with open(cls.testfile,'r') as F:
            txt = F.read()
        cls.objects = [Link716(i) for i in txt.split(SMARK)]

    def test_init(self):
        msg = 'Incorrect parsing of Link716'
        attrs = Link716.__slots__
        for obj in self.objects:
            self.assertTrue(bool(obj.text),msg)
            self.assertEqual(obj.number,716,msg)
        for attr in attrs:
            self.assertTrue(getattr(self.objects[2],attr),msg)

    def test_regex_Thermo(self):
        msg = 're_Thermo does not match. \n{}\n!=\n{}'.format
        msg2 = 're_Thermo matches when it should not {}'.format
        regex = Link716.re_Thermo
        with open(self.Thermofile,'r') as F:
            txt = F.read()
        solutions = txt.split(SMARK)
        for obj,solution in zip(self.objects,solutions):
            test = regex.findall(obj.text)
            if solution != '':
                test = test[0]
                self.assertTrue(test == solution,msg(test,solution))
            else:
                self.assertFalse(bool(test),msg2(test))
    def test_regex_EContrib(self):
        msg = 're_EContrib does not match. \n{}\n!=\n{}'.format
        msg2 = 're_EContrib matches when it should not {}'.format
        regex = Link716.re_EContrib
        with open(self.EContribfile,'r') as F:
            txt = F.read()
        solutions = txt.split(SMARK)
        for obj,solution in zip(self.objects,solutions):
            test = regex.findall(obj.text)
            if solution != '':
                test = test[0]
                self.assertTrue(test == solution,msg(test,solution))
            else:
                self.assertFalse(bool(test),msg2(test))
    def test_regex_IRSpectrum(self):
        msg = 're_IRSpectrum does not match. \n{}\n!=\n{}'.format
        msg2 = 're_IRSpectrum matches when it should not {}'.format
        regex = Link716.re_IRSpectrum
        with open(self.IRSpectrumfile,'r') as F:
            txt = F.read()
        solutions = txt.split(SMARK)
        for obj,solution in zip(self.objects,solutions):
            test = regex.findall(obj.text)
            if solution != '':
                test = test[0]
                self.assertTrue(test == solution,msg(test,solution))
            else:
                self.assertFalse(bool(test),msg2(test))
    def test_regex_Frequencies(self):
        msg = """re_Frequency does not properly match.'.format
        bool(
        {}
        )
        !=bool(
        {}
        )""".format
        msg2 = 'Frequency slice does not match. \n{}\n!=\n{}'.format
        regex = Link716.re_Frequencies
        with open(self.Frequenciesfile,'r') as F:
            txt = F.read()
        samples = txt.split(SMARK)
        solutions = []
        for sample in samples:
            if sample:
                Matches = sample.split('\n## Match ##\n')
                Aux = []
                for match in Matches:
                    Aux.append(tuple(match.split('\n')))
                solutions.append(Aux)
            else:
                solutions.append([])
        for obj,solution in zip(self.objects,solutions):
            test = regex.findall(obj.text)
            self.assertTrue(bool(test) == bool(solution),msg(test,solution))
            for t,s in zip(test,solution):
                self.assertTrue(t == s, msg2(t,s))
    def test_regex_freq_text(self):
        msg = 'the text block does not match'
        regex = Link716.re_freq_text
        with open(self.freqtextfile,'r') as F:
            txt = F.read()
        solutions = txt.split(SMARK)
        for obj,solution in zip(self.objects,solutions):
            test = regex.findall(obj.text)[0]
            self.assertTrue(test == solution,msg)
    def test_regex_freq_text(self):
        msg = 'the displacement: \n{} does not match:\n{}'.format
        regex = Link716.re_freq_displacements
        regex_txt = Link716.re_freq_text
        with open(self.freqdispfile,'r') as F:
            txt = F.read()
        items = txt.split(SMARK)
        solutions = []
        for item in items: 
            if item:
                solutions.append([s+'\n' for s in item.split('\n')])
            else:
                solutions.append([])
        for obj,solution in zip(self.objects,solutions):
            test = regex_txt.findall(obj.text)
            if test:
                subtests = regex.findall(test[0])
                for t,s in zip(subtests,solution): 
                    self.assertEqual(t,s,msg(test,solution))
            else:
                self.assertTrue(test == solution,msg(test,solution))

    def test_regex_dipole(self):
        msg = 're_dipole does not match. \n{}\n!=\n{}'
        msg = msg.format
        regex = Link716.re_dipole
        with open(self.dipolefile,'r') as F:
            txt = F.read()
        samples = txt.split(SMARK)
        solutions = []
        for sample in samples:
            solutions.append(tuple(sample.split('\n')))
        for obj,solution in zip(self.objects,solutions):
            test = regex.findall(obj.text)[0]
            self.assertTrue(test == solution,msg(test,solution))

    def test_init_empty(self):
        msg = 'Incorrect empty initialization of Link716'
        obj = Link716('Some Text',asEmpty=True)
        self.assertFalse(bool(obj.text),msg)
        self.assertEqual(obj.number,716)

class TestLink804(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.testfile = TEST_FILEDIR.joinpath('l804.txt')
        cls.propfile = TEST_FILEDIR.joinpath('l804_properties.txt')
        # Read and store the Examples
        with open(cls.testfile,'r') as F:
            txt = F.read()
        cls.objects = [Link804(i) for i in txt.split(SMARK)]

    def test_init(self):
        msg = 'Incorrect parsing of Link804'
        for obj in self.objects:
            self.assertTrue(bool(obj.text),msg)
            self.assertEqual(obj.number,804)
            self.assertTrue(obj.MP2 is not None,msg)
            self.assertTrue(bool(obj.SpinComponents),msg)

    def test_regex_MP2(self):
        msg = 're_MP2 regex does not match properly'
        regex = Link804.re_MP2
        with open(self.propfile,'r') as F:
            txt = F.read()
        solutions = [lines.split('\n')[0] for lines in txt.split(SMARK)]
        items = [obj.text for obj in self.objects]
        for txt,solution in zip(items,solutions):
            test = regex.findall(txt)[0]
            self.assertTrue(test == solution,msg)

    def test_init_empty(self):
        msg = 'Incorrect empty initialization of Link804'
        text = 'Link1:  Proceeding to internal job step number 13.'
        obj = Link804(text,asEmpty=True)
        self.assertFalse(bool(obj.text),msg)
        self.assertEqual(obj.number,804)
        self.assertTrue(obj.MP2 is None,msg)
        self.assertFalse(bool(obj.SpinComponents),msg)

    def test_SpinComponents(self):
        msg = 'SpinComponents not properly parsed'
        SpinComponent = Link804._SpinComponent
        with open(self.propfile,'r') as F:
            txt = F.read().replace('D','E')
        items = [lines for lines in txt.split(SMARK)]
        for item,obj in zip(items,self.objects):
            sol = []
            lines = item.split('\n')[1:]
            for line in lines:
                a0,a1,a2 = line.split()
                sol.append(SpinComponent(a0,float(a1),float(a2)))
            self.assertTrue(obj.SpinComponents == sol,msg)

    def test_SCSCorr(self):
        msg = 'SCS correction not properly calculated'
        MockObject = Link804('',asEmpty=True)
        SpinComponent = Link804._SpinComponent
        MockObject.SpinComponents = [None,None,None]
        MockDict = dict(Name='',E=0,T=0)
        for item in [(0,0,0),(3,5,3),(1.5,10,1.5)]:
            for i,j in enumerate(item):
                MockDict['E'] = j
                MockObject.SpinComponents[i] = SpinComponent(**MockDict)
            test = MockObject.get_SCScorr()
            aa,ab,bb = item
            sol = (aa+bb)/3.0 +(6.0*ab)/5.0
            self.assertTrue(test == sol,msg)

    def test_MP2(self):
        msg = 'MP2 Energy not properly parsed'
        with open(self.propfile,'r') as F:
            txt = F.read().replace('D','E')
        solutions = [float(lines.split('\n')[0]) for lines in txt.split(SMARK)]
        for sol,obj in zip(solutions,self.objects):
            self.assertTrue(obj.MP2 == sol,msg)

class TestLink913(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.testfile = TEST_FILEDIR.joinpath('l913.txt')
        cls.MP4file = TEST_FILEDIR.joinpath('l913_MP4.txt')
        cls.ccsdtfile = TEST_FILEDIR.joinpath('l913_ccsdt.txt')
        # Read and store the Examples
        with open(cls.testfile,'r') as F:
            txt = F.read()
        cls.objects = [Link913(i) for i in txt.split(SMARK)]

    def test_init(self):
        msg = 'Incorrect parsing of Link913'
        for obj in self.objects:
            self.assertTrue(bool(obj.text),msg)
            self.assertEqual(obj.number,913)

    def test_regex_MP4(self):
        msg = 're_MP4 regex does not match properly'
        regex = Link913.re_MP4
        with open(self.MP4file) as F:
            txt = F.read()
        solutions = [i.split('\n') for i in txt.split(SMARK)]
        items = [obj.text for obj in self.objects]
        for txt,solution in zip(items,solutions):
            test = regex.findall(txt)
            # Test only for the last match as it is the only one stored
            self.assertTrue(test[-1] == solution[-1],msg)

    def test_regex_ccsdt(self):
        msg_match = 're_ccsdt regex does not match properly\n{} != {}'
        msg_nomatch = 're_ccsdt regex does not find a match properly'
        regex = Link913.re_ccsdt
        with open(self.ccsdtfile) as F:
            txt = F.read()
        solutions = [i for i in txt.split(SMARK)]
        items = [obj.text for obj in self.objects]
        for txt,solution in zip(items,solutions):
            test = regex.findall(txt)
            if test: # Test matches
                self.assertTrue(test[0] == solution,msg_match.format(test[0],solution))
            else: # Test non-matches, although the test is within the condition
                self.assertFalse(bool(test),msg_nomatch)

    def test_init_empty(self):
        msg = 'Incorrect empty initialization of Link913'
        obj = Link913('Some Text',asEmpty=True)
        self.assertFalse(bool(obj.text),msg)
        self.assertEqual(obj.number,913)

    def test_MP4(self):
        msg = 'MP4 energy not properly read'
        with open(self.MP4file) as F:
            txt = F.read().replace('D','E')
        solutions = [float(lines.split('\n')[-1].split()[-1])
                                        for lines in txt.split(SMARK)]
        for sol,obj in zip(solutions,self.objects):
            self.assertTrue(obj.MP4 == sol,msg)

    def test_CCSDT(self):
        msg = 'CCSDT energy not properly read'
        with open(self.ccsdtfile,'r') as F:
            txt = F.read().replace('D','E')
        solutions = [float(lines.strip()) if lines.strip() else None
                                        for lines in txt.split(SMARK)]
        for sol,obj in zip(solutions,self.objects):
            self.assertTrue(obj.CCSDT == sol,msg)

class TestLink914(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        ES_mark = '\n## ExcitedState ##\n'
        cls.testfile = TEST_FILEDIR.joinpath('l914.txt')
        cls.propertiesfile = TEST_FILEDIR.joinpath('l914_properties.txt')
        # Read and store the Examples
        with open(cls.testfile,'r') as F:
            txt = F.read()
        cls.objects = [Link914(i) for i in txt.split(SMARK) if i.strip()]
        with open(cls.propertiesfile,'r') as F:
            txt = F.read()
        cls.solutions = [slice.split(ES_mark) for slice in txt.split(SMARK) if slice.strip()]

    def test_init(self):
        msg = 'Incorrect parsing of Link914'
        for obj in self.objects:
            self.assertTrue(bool(obj.text),msg)
            self.assertEqual(obj.number,914)

    def test_regex_ExcitedState(self):
        msg = 're_ExcitedState regex does not match {} with {}'.format
        regex = Link914.re_ExcitedState
        solutions = []
        for sol in self.solutions:
            aux = []
            for ES in sol:
                aux.append(tuple(ES.split('\n')[0].strip().split()))
            solutions.append(aux)
        items = [obj.text for obj in self.objects]
        for txt,solution in zip(items,solutions):
            test = regex.findall(txt)
            # Test only for the last match as it is the only one stored
            self.assertTrue(tuple(test) == tuple(solution),msg(test,solution))

    def test_regex_transition(self):
        msg = 're_transition regex does not match {} with {}'.format
        regex = Link914.re_transition
        # >> obj = solutions[0] = [[('111A','->','222B','-0.00000'),(...)], ...]
        # >> ES  = solutions[0][0] = [('111A','->','222B','-0.00000'),(...)]
        # >> orb = solutions[0][0][0] = ('111A','->','222B','-0.00000')
        solutions = []
        for sol in self.solutions:
            aux = []
            for ES in sol:
                transitions = []
                for transition in ES.split('\n')[1:]:
                    if transition.strip():
                        transitions.append(tuple(transition.strip().split()))
                aux.append(transitions)
            solutions.append(aux)
        items = [obj.text for obj in self.objects]
        for txt,solution in zip(items,solutions):
            tests = regex.findall(txt)
            for test,sol in zip(tests,list(chain.from_iterable(solution))):
                self.assertTrue(test == sol,msg(test,sol))

    def test_init_empty(self):
        msg = 'Incorrect empty initialization of Link913'
        obj = Link914('Some Text',asEmpty=True)
        self.assertFalse(bool(obj.text),msg)
        self.assertEqual(obj.number,914)

    def test_locate_ExcitedStates(self):
        msg = 'ExcitedStates not read properly {} does not match {}'.format
        solutions = []
        for sol in self.solutions:
            aux = []
            for ES in sol:
                a,_,b,c,d,e = ES.split('\n')[0].strip().split()
                a = int(a)
                b,c,d,e = tuple(map(float,(b,c,d,e)))
                aux.append((a,b,c,d,e))
            solutions.append(aux)
        for solution,obj in zip(solutions,self.objects):
            for sol,ES in zip(solution,obj.excitedstates):
                test = tuple(ES[:-1])
                self.assertTrue(sol == test ,msg(sol,test))

    def test_extract_transitions(self):
        msg = 'Wrong Transition Extraction: {1} does not match {0}'.format
        L914 = Link914('',asEmpty=True)
        txt = """
             107B ->166B        0.02496
             165B ->167B        0.02444
             126B <-166B       -0.01650
             Sometexthere
             ## Split Here ##
             165 -> 167        0.02444
             126 <- 166       -0.01650

             """
        solutions = [[('107B','166B',float('0.02496' ),False),
                      ('165B','167B',float('0.02444' ),False),
                      ('166B','126B',float('-0.01650'),True)],
                     [( '165', '167',float('0.02444' ),False),
                      ( '166', '126',float('-0.01650'),True)]]
        for sol,intxt in zip(solutions,txt.split(SMARK)):
            transitions = L914._extract_transitions(intxt)
            for sol,test in zip(sol,transitions):
                self.assertEqual(tuple(test),sol,msg(test,sol))


if __name__ == '__main__':
    unittest.main()
