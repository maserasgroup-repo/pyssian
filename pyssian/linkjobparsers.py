"""
One of the two core modules of pyssian. This module contains the classes and
parsers for the outputs of the different types of links. Only the classes
registered with the appropiate decorator are considered for parsing.
"""

import warnings
from collections import namedtuple
from itertools import cycle
import re

# Auxiliar Documenting Decorators
def Populates(*attributes, defaults=None):
    """
    Helper decorator function to document _locate* methods

    Parameters
    ----------
    *attributes : str
        names of the attributes populated by the method.
    default : list
        list of default values in the same order as the attributes. If none
        provided no "default message" is added.

    Returns
    -------
    subdecorator
        returns a decorator that attaches a message constructed
        with the attributes and defaults specified
    """

    if len(attributes) > 0:
        if defaults is None:
            Out = [f'- {attr}' for attr in attributes]
        elif len(defaults) > 1:
            Out = [f'- {attr} ({default})'
                    for attr, default in zip(attributes,defaults)]
        else:
            Out = [f'- {attr} ({defaults[0]})' for attr in attributes]
        Out.insert(0,'Populates the attributes:')
        Out.append('')
    else:
        Out = [f'Populates the attribute: **{attributes[0]}**',]
    def subdecorator(f):
        match = re.search(r'^[\s]*',f.__doc__)
        lpadding = '\n'
        if match:
            lpadding += match[0]
        f.__doc__ += lpadding + lpadding.join(Out) + '\n'
        return f
    return subdecorator
def SilentFail(f):
    f.__doc__ += '(Fails Silently).'
    return f

#Base Class
class LinkJob(object):
    """LinkJob Base class. Represents the output of a lxxxx.exe.

    Parameters
    ----------
    text : list
        text that corresponds to the output of the lxxxx.exe
    number : int
        Integrer that represents the type of the LinkJob (the default is None).

    Attributes
    ----------
    Register : dict
        A dict where key corresponds to an int, and value corresponds to the
        registered subclass class through @RegisterLinkJob
    number

    text

    """
    __slots__ = ('text','number')
    _token = 'Base'
    Register = {}
    def __init__(self,text,number=None,asEmpty=False):
        self.number = number
        if asEmpty:
            text = ''
        self.text = text
    def __repr__(self):
        cls = type(self).__name__
        return f'<Link {self.number} of class {cls}>'
    def __str__(self):
        return ''.join(self.text)

    def _text_iterbyline(self):
        return [line+'\n' for line in self.text.split('\n')].__iter__()

    @classmethod
    def as_empty(cls,text,number=None):
        instance = cls(text,asEmpty=True)
        if instance.number is None:
            instance.number = cls._token
        return instance

class GeneralLinkJob(LinkJob):
    """
    Subclass of Linkjob to store the information of any Link that does not
    have a specific parser.

    Parameters
    ----------

    text : str
        text that corresponds to the output of the lxxxx.exe
    number : int
        Integrer that represents the type of the LinkJob. (defaults to None)
    asEmpty : bool
        Flag to not parse and store the information of the text.
        (defaults to False)

    Attributes
    ----------

    number

    text

    """

    _token = 'General'

    __slots__ = ('number', 'text')

    re_enter = re.compile(r'(?:Enter.*l)([0-9]{1,4})(?:\.exe)')

    def __init__(self,text,number=None,asEmpty=False):
        cls = self.__class__
        if number is None:
            Match = cls.re_enter.findall(text)
            if Match:
                number = int(Match[0])
            else:
                number = -1
        if asEmpty:
            text = ''
        super().__init__(text,number)

def RegisterLinkJob(cls):
    """Decorator for the LinkJob subclasses, adds the class to the
    LinkJob.Register dictionary with the _token as key.

    Parameters
    ----------
    cls : LinkJob
        LinkJob subclass identified with the cls._token attribute. Which is
        normally an integrer between 1 and 9999.

    """
    key = cls._token
    if key in LinkJob.Register:
        msg = '\nWarning! overwriting :\n{} with {}\n'.format
        warnings.warn(msg(LinkJob.Register[key],cls))
    LinkJob.Register[key] = cls
    return cls

@RegisterLinkJob
class Link1(LinkJob):
    """
    Representation and parser for the output of l1.exe,it also holds
    information about the InternalJob setup.

    Parameters
    ----------
    text : str
        text that corresponds to the output of the l1.exe
    asEmpty : bool
        Flag to not parse and store the information of the text.
        (defaults to False)

    Attributes
    ----------
    info : namedtuple
        tuple that contains information needed to instantiate an InternalJob.
    commandline : str
        contains the information of the gaussian commands
        (defaults to an empty string).
    nprocs : int
        contains the number of the %nprocshared link0 option
    mem : str
        contains the right hand side of the %mem link0 option
    link0 : list
        list of lines containing the link0 specifications.
    IOps : list
        list of lines containing the IOps that appear below the command line
    text
    number
    """

    _token = 1

    re_link0 = re.compile(r'\%.*')
    re_commandline = re.compile(r'(?:\-{10,}\n)(.*?#[\s\S]*?)(?:\-{10,}\n)')
    re_IOps = re.compile(r'.*\;')
    re_internaljob = re.compile(r'Link1\:.*internal.*job.*([0-9]{1,3})\.')

    InternalJobInfo = namedtuple('InternalJobInfo','number type new_InternalJob')

    def __init__(self,text,asEmpty=False):
        super().__init__(text,1)
        self.info = None
        self.commandline = ''
        self.nprocs = None
        self.mem = None
        self.link0 = []
        self.IOps = []
        if asEmpty:
            self._locate_internaljob()
            self.text = '' # clear text
        else:
            self._locate_link0()
            self._locate_commandline()
            self._locate_IOps()
            self._locate_internaljob()

    @Populates('info',defaults=['If none is found defaults to 1'])
    def _locate_internaljob(self):
        """
        Uses regex expressions compiled as class attributes to search the number
        of the "current" internal job at the file.
        """
        cls = self.__class__
        Match = cls.re_internaljob.findall(self.text)
        if Match:
            JobNumber = int(Match[0])
            New_InternalJob = True
        else:
            JobNumber = 1
            New_InternalJob = False
        jobtype = self._guess_type()
        self.info = cls.InternalJobInfo(JobNumber, jobtype, New_InternalJob)

    @Populates('nprocs','mem', 'link0', defaults=['Defaults to None',
    'Defaults to None','If none is found defaults to an empty list'])
    def _locate_link0(self):
        """
        Uses regex expressions compiled as class attributes to search the link0
        options (specified with a % at the start of the option).
        """
        link0_opt = []
        cls = self.__class__
        Match = cls.re_link0.findall(self.text)
        if Match:
            for item in Match:
                if 'nprocshared' in item:
                    self.nprocs = int(item.strip().split('=')[-1])
                elif 'mem' in item:
                    self.mem = item.strip().split('=')[-1]
                else:
                    link0_opt.append(item[1:])
        self.link0 = link0_opt

    @Populates('commandline')
    @SilentFail
    def _locate_commandline(self):
        """
        Uses regex expressions compiled as class attributes to search the
        gaussian commands, which start with a "#".
        """
        cls = self.__class__
        Match = cls.re_commandline.findall(self.text)
        if Match:
            lines = Match[0].split('\n')
            commandline = ''.join([line.lstrip() for line in lines])
            self.commandline = commandline

    @Populates('IOps')
    @SilentFail
    def _locate_IOps(self):
        """
        Uses regex expressions compiled as class attributes to search IOps (see
        Gaussian Manual).
        """
        cls = self.__class__
        Match = cls.re_IOps.findall(self.text)
        if Match:
            self.IOps = [line for line in Match]

    def _guess_type(self):
        """
        Guesses the type of Job depending on the keywords in the
        self.commandline.

        Returns
        -------
        str
            {Constrained Optimization, Optimization,Frequency Calculation, 
            Unidentified, Linked} the Linked type is only assigned externally
        """
        target = self.commandline.lower()
        if 'modredundant' in target:
            jobtype = 'Constrained Optimization'
        elif 'opt' in target:
            jobtype = 'Optimization'
        elif 'freq' in target:
            jobtype = 'Frequency Calculation'
        else:
            jobtype = 'Unidentified'
        return jobtype

@RegisterLinkJob
class Link101(LinkJob):
    """
    Representation and parser for the output of l101.exe

    Parameters
    ----------
    text : str
        text that corresponds to the output of the l101.exe
    asEmpty : bool
        Flag to not parse and store the information of the text.
        (defaults is False)

    Attributes
    ----------
    charge : int
    spin : int
    """
    __slots__ = ('charge','spin')

    _token = 101

    re_charge = re.compile(r'(?:Charge\s=\s*)(-{0,1}[0-9]*)')
    re_spin = re.compile(r'(?:Multiplicity\s=\s*)([0-9]*)')

    def __init__(self,text,asEmpty=False):
        self.charge = None
        self.spin = None
        if asEmpty:
            super().__init__('',101)
        else:
            super().__init__(text,101)
            self._locate_charge_spin()
    @Populates('charge','spin')
    @SilentFail
    def _locate_charge_spin(self):
        """
        Searches the charge and spin.
        """
        cls = self.__class__
        charge = cls.re_charge.findall(self.text)
        spin = cls.re_spin.findall(self.text)
        if charge:
            self.charge = int(charge[0])
        if spin:
            self.spin = int(spin[0])

@RegisterLinkJob
class Link103(LinkJob):
    """
    Representation and parser for the output of l103.exe. Corresponds to a
    single step of a Berny optimization.

    Parameters
    ----------
    text : str
        text that corresponds to the output of the l103.exe
    asEmpty : bool
        Flag to not parse and store the information of the text.
        (defaults to False)

    Attributes
    ----------
    mode : str
        'Init', 'Iteration' or 'End'
    state : str
        Either 'Optimized', 'Non-optimized' or 'Initial'
    conversion : list
        List of namedtuples with fields Item,Value,Threshold,Converged.
    parameters : list
        List of namedtuples with fields Name,Definition,Value,Derivative.
    derivatives : list
        List of namedtuples with fields Var,old,dEdX,dXl,dXq,dXt,new.
    stepnumber
    scanpoint
    text
    number

    """

    __slots__ = ('mode','state','conversion','parameters','derivatives',
                'stepnumber','scanpoint')

    _token = 103

    re_init = re.compile(r'Initialization\spass\.')
    re_end = re.compile(r'(Optimization\sstopped\.)|(Optimization completed\.)')
    re_parameters = re.compile(r'\!.*\!')
    re_derivatives = r'([R|D|X|A][0-9]*)'+r'\s*([-]?[0-9]*\.[0-9]{0,5})'*6
    re_derivatives = re.compile(re_derivatives)
    re_convergence = r'((?:Maximum|RMS)\s*(?:Force|Displacement))'
    re_convergence += r'\s*([0-9]{0,3}\.[0-9]*)\s*([0-9]{0,3}\.[0-9]*)\s*(YES|NO)'
    re_convergence = re.compile(re_convergence)
    re_stepnum = re.compile(r'Step\s*number\s*([0-9]*)\s*out\s*of')
    re_scanpoint = re.compile(r'scan\spoint\s*([0-9]*)\s*out\s*of')

    _Parameter = namedtuple("Parameter", "Name Definition Value Derivative")
    _ConverItem = namedtuple("ConverItem", "Item Value Threshold Converged")
    _Derivative = namedtuple("Derivative", "Var old dEdX dXl dXq dXt new")

    def __init__(self,text,asEmpty=False):
        self.mode = None
        self.state = None
        self.conversion = []
        self.parameters = []
        self.derivatives = []
        self.stepnumber = None
        self.scanpoint = None
        if asEmpty:
            super().__init__('',103)
        else:
            super().__init__(text,103)
            self._locate_mode()
            self._locate_parameters()
            self._locate_derivatives()
            self._locate_convergence()
            self._locate_numbers()

    @Populates('mode',defaults=['Iteration'])
    def _locate_mode(self):
        """
        Uses regex expressions compiled as class attributes to find
        If the Link corresponds to an (Init)ialization pass, (Iteration)
        or (End) converged or not.
        """
        cls = self.__class__ # Alias
        IsInit = cls.re_init.findall(self.text)
        IsEnd  = cls.re_end.findall(self.text)
        if IsInit:
            self.mode = 'Init'
        elif IsEnd:
            self.mode = 'End'
        else:
            self.mode = 'Iteration'

    @Populates('state','parameters')
    @SilentFail
    def _locate_parameters(self):
        """
        Uses regex expressions compiled as class attributes to find
        the parameter table if exists.
        """
        cls = self.__class__
        Match = cls.re_parameters.findall(self.text)
        Parameters = []
        if Match:
            state = Match.pop(0)
            # '! {Initial/Optimized/Non-Optimized} Parameters !'
            self.state = state.replace('!','').strip().split()[0]
            _ = Match.pop(0) # Discard Units
            _ = Match.pop(0) # Discard Table Headings
            for line in Match:
                Columns = line.strip().split()
                Name = Columns[1]
                Definition = Columns[2]
                Value = Columns[3]
                if self.mode == "Init":
                    Derivative = Columns[-2]
                else:
                    Derivative = float(Columns[-2])
                Par = cls._Parameter(Name, Definition,
                                    Value, Derivative)
                Parameters.append(Par)
            self.parameters = Parameters

    @Populates('state','derivatives')
    @SilentFail
    def _locate_derivatives(self):
        """
        Uses regex expressions compiled as class attributes to find
        the derivatives table if exists.
        """
        cls = self.__class__
        Match = cls.re_derivatives.findall(self.text)
        items = []
        if Match:
            for m in Match:
                Var, old, dEdX, dXl, dXq, dX, new = m
                item = cls._Derivative(Var, old, dEdX, dXl, dXq, dX, new)
                items.append(item)
            self.derivatives = items

    @Populates('conversion')
    @SilentFail
    def _locate_convergence(self):
        """
        Uses regex expressions compiled as class attributes to find
        the convergence parameters.
        """
        cls = self.__class__
        Match = cls.re_convergence.findall(self.text)
        items = []
        if Match:
            for m in Match:
                Name = ' '.join(m[0].split())
                val = float(m[1])
                threshold = float(m[2])
                isConverged = m[3] == "YES"
                item = cls._ConverItem(Name,val,threshold,isConverged)
                items.append(item)
            self.conversion = items

    @Populates('stepnumber','scanpoint')
    @SilentFail
    def _locate_numbers(self):
        """
        Uses regex expressions compiled as class attributes to find
        the current step number (of an optimization) and scan number
        (of an scan calculation).
        """
        cls = self.__class__
        Match = cls.re_stepnum.findall(self.text)
        if Match:
            self.stepnumber = int(Match[0])
        elif self.mode == 'Init':
            self.stepnumber = 0
        Match = cls.re_scanpoint.findall(self.text)
        if Match:
            self.scanpoint = int(Match[0])

    def print_convergence(self):
        """ Prints the convergence Table formatted """
        Format = '{0: <22s}\t{1:.6f}\t{2:.6f}\t{3}'.format
        try:
            Out = [Format(*i) for i in self.conversion]
        except ValueError as e:
            print("Non numeric value found in the convergence values")
            Format_str='{0: <22s}\t{1}\t{2}\t{3}'.format
            Out = [Format_str(*i) for i in self.conversion]
        print('\n'.join(Out))

@RegisterLinkJob
class Link120(LinkJob):
    """
    Representation and parser for the output of l120.exe. Corresponds to a
    the energy calculations when using ONIOM.

    Parameters
    ----------
    text : str
        text that corresponds to the output of the l103.exe
    asEmpty : bool
        Flag to not parse and store the information of the text.
        (defaults to False)

    Attributes
    ----------
    energy
    energy_partitions

    """
    __slots__ = ('energy', 'energy_partitions')

    _token = 120

    re_energy = r'extrapolated energy\s=\s*(-?[0-9]*\.[0-9]*)'
    re_energy = re.compile(re_energy)
    re_energy_partitions = r'gridpoint\s*([0-9]*)\smethod:\s*([a-zA-Z]*)\s*system:\s*([a-zA-Z]*)\s*energy:\s*(-?[0-9]*\.[0-9]*)'
    re_energy_partitions = re.compile(re_energy_partitions)
    _EnergyPartition = namedtuple('EnergyPartition','gridpoint level model energy')

    def __init__(self,text,asEmpty=False):
        self.energy = None
        self.energy_partitions = []
        if asEmpty:
            super().__init__('',120)
        else:
            super().__init__(text,120)
            self._locate_energies()

    @Populates('energy','energy_partitions',defaults=['None','[]'])
    @SilentFail
    def _locate_energies(self):
        """
        Uses regex expressions compiled as class attributes to find
        the energy and the different energy_partitions.
        """
        cls = self.__class__
        match = cls.re_energy.findall(self.text)
        if match:
            self.energy = float(match[0])

        EnergyPartition = cls._EnergyPartition
        match = cls.re_energy_partitions.findall(self.text)
        if not match:
            return 
        self.energy_partitions = [EnergyPartition(int(p),level,model,float(energy))
                                  for p,level,model,energy in match]

@RegisterLinkJob
class Link123(LinkJob):
    """
    Representation and parser for the output of l123.exe. Corresponds to irc
    predictor-corrector steps and may hold the orientation (Geometry) at a given
    step.

    Parameters
    ----------
    text : str
        text that corresponds to the output of the l202.exe
    asEmpty : bool
        Flag to not parse and store the information of the text.
        (defaults to False)

    Attributes
    ----------
    orientation : list
        List of namedtuples containing the geometry of the calculated compound.
    step : int
        step number of along the irc calculation.
    direction : str
        FORWARD or REVERSE
    reactioncoord : float
       Value of the reaction coordinate at the current step.

    text
    number

    """
    __slots__ = ('orientation', 'direction', 'step', 'reactioncoord')

    _token = 123

    re_orientation =  r'^\s*([0-9]{1,4})\s*([0-9]{1,3})\s*'
    re_orientation += r'(\-?[0-9]+\.[0-9]*)'
    re_orientation += r'\s*(\-?[0-9]+\.[0-9]*)'*2
    re_orientation = re.compile(re_orientation, re.MULTILINE)

    re_reactioncoord =  r'NET REACTION COORDINATE UP TO THIS POINT\s*='
    re_reactioncoord += r'\s*([0-9]*\.[0-9]*)'
    re_reactioncoord = re.compile(re_reactioncoord)

    re_step = r'Point Number\s*([0-9]*)\s*in\s([A-Z]*)\spath direction'
    re_step = re.compile(re_step)

    re_iscomplete = r'(Calculation of (REVERSE|FORWARD) path complete.)'
    re_iscomplete = re.compile(re_iscomplete)

    re_totalsteps = r'\#\sOF\sPOINTS\sALONG\sTHE\sPATH\s*=\s*([0-9]*)'
    re_totalsteps = re.compile(re_totalsteps)

    _AtomCoords = namedtuple('AtomCoords',
                            'CenterNum AtomicNum AtomType X Y Z')

    def __init__(self,text,asEmpty=False):
        self.orientation = []
        self.direction = ''
        self.step = 0
        self.reactioncoord = 0
        if asEmpty:
            super().__init__('',123)
        else:
            super().__init__(text,123)
            self._locate_orientation()
            self._locate_irc_step()

    @Populates('orientation')
    @SilentFail
    def _locate_orientation(self):
        """
        Uses regex expressions compiled as class attributes to find
        the geometry of the molecule.
        """
        cls = self.__class__
        AtomCoords = cls._AtomCoords
        match = cls.re_orientation.findall(self.text)
        if match:
            for line in match:
                Atom = [int(line[0]),int(line[1]),None]
                Atom.extend([float(i) for i in line[2:]])
                self.orientation.append(AtomCoords(*Atom))

    @Populates('reactioncoord','step','direction')
    @SilentFail
    def _locate_irc_step(self):
        """
        Uses regex expressions compiled as class attributes to find
        the current step info of the irc path.
        """
        cls = self.__class__
        match = cls.re_reactioncoord.findall(self.text)
        if match:
            self.reactioncoord = float(match[0])
        match = cls.re_step.findall(self.text)
        if match:
            self.step = int(match[0][0])
            self.direction = match[0][1]
        else:
            iscomplete = cls.re_iscomplete.findall(self.text)
            if iscomplete:
                self.direction = iscomplete[0][1]
                totalsteps = int(cls.re_totalsteps.findall(self.text)[0])
                self.step = totalsteps + 1

    def print_orientation(self):
        """
        Displays in console the geometry with a format similar to gaussian.
        """
        Format = '{0}\t{1}\t{3: 0.6f}\t{4: 0.6f}\t{5: 0.6f}\n'
        Header = 'Center\tAtNum\t{0:^9s}\t{1:^9s}\t{2:^9s}'.format('X','Y','Z')
        print(Header)
        Out = [Format.format(*i) for i in self.orientation]
        print(''.join(Out))

@RegisterLinkJob
class Link202(LinkJob):
    """
    Representation and parser for the output of l202.exe. Holds the
    orientation (Geometry) at a given step.

    Parameters
    ----------
    text : str
        text that corresponds to the output of the l202.exe
    asEmpty : bool
        Flag to not parse and store the information of the text.
        (defaults to False)

    Attributes
    ----------
    orientation : list
        List of namedtuples containing the geometry of the calculated compound.
    DistanceMatrix : list
        For each item the postitions 0 and 1 correspond to the Atom ID in the
        structure and the Atomic Number/Symbol respectively. Position 1 holds
        a list with the distances with respect to all previous items.
    text
    number

    """
    __slots__ = ('orientation', 'DistanceMatrix')

    _token = 202

    re_orientation = r'(?:.*(?:Standard|Input)\sorientation\:.*\n'
    re_orientation += r'\s?\-+\n.*\n.*\n\s?\-+\n)'
    re_orientation += r'[\s\S]*?'
    re_orientation += r'(?:\s?\-+\n)'
    re_orientation = re.compile(re_orientation)
    re_number = re.compile('[0-9]')

    _AtomCoords = namedtuple('AtomCoords',
                            'CenterNum AtomicNum AtomType X Y Z')

    def __init__(self,text,asEmpty=False):
        self.orientation = []
        self.DistanceMatrix = []
        if asEmpty:
            super().__init__('',202)
        else:
            super().__init__(text,202)
            self._locate_orientation()
            self._locate_DistanceMatrix()

    @Populates('orientation')
    @SilentFail
    def _locate_orientation(self):
        """
        Uses regex expressions compiled as class attributes to find
        the geometry of the molecule.
        """
        cls = self.__class__
        AtomCoords = cls._AtomCoords
        match = cls.re_orientation.findall(self.text)
        if match:
            # If Input and Standard orientation coexist
            # take the last one (Standard, generally)
            lines = [i for i in match[-1].split('\n') if i.strip()]
            _ = lines.pop(0) # Discard "Input orientation" line
            # Discard "Header"
            for i in range(4):
                _ = lines.pop(0)
            # Discard Last Line that
            _ = lines.pop(-1)
            for line in lines:
                Aux = line.split()
                Atom = [int(i) for i in Aux[0:3]]
                Atom.extend([float(i) for i in Aux[3:]])
                self.orientation.append(AtomCoords(*Atom))

    @Populates('DistanceMatrix')
    @SilentFail
    def _locate_DistanceMatrix(self):
        """
        Looks for the keyword 'Distance matrix' and reads the
        Input orientation table that follows
        """
        Iterator = self._text_iterbyline()
        cls = self.__class__
        for line in Iterator:
            if "Distance matrix" in line:
                break
        else: # Fail Silently
            return
        for _ in range(2):
            line = next(Iterator)
        Items = []
        while cls.re_number.match(line.strip()[0]):
            # while the first non blank character is a number
            Atom = line.strip().split()
            Item = []
            if '.' not in line or cls.re_number.match(Atom[1]):
                line = next(Iterator)
                continue
            Item.append(int(Atom[0]))
            Item.append(Atom[1])
            Item.append(list(map(float,Atom[2:])))
            Items.append(Item)
            line = next(Iterator)
        DicMatrix = {i[0]:['X',[]] for i in Items}
        for i in Items:
            Atom = DicMatrix[i[0]]
            if Atom[0] == "X":
                Atom[0] = i[1]
            else:
                if Atom[0] != i[1]:
                    raise RuntimeError('''Two atoms with the same Center number
                                     and different Symbol have been found ''')
            Atom[1].extend(i[2])
        self.DistanceMatrix = []
        for Id, (Sym, Dist) in DicMatrix.items():
            self.DistanceMatrix.append([Id,Sym,Dist])

    def print_orientation(self):
        """
        Displays in console the geometry with a format similar to gaussian.
        """
        Format = '{0}\t{1}\t{3: 0.6f}\t{4: 0.6f}\t{5: 0.6f}\n'.format
        Header = 'Center\tAtNum\t{0:^9s}\t{1:^9s}\t{2:^9s}'.format('X','Y','Z')
        print(Header)
        Out = [Format(*i) for i in self.orientation]
        print(''.join(Out))
    def get_atom_mapping(self):
        """
        Returns a dictionary that relates the atomic number with the symbol
        based on the information contained in the orientation and distance
        matrix.

        Returns
        -------
        Out
            Dictionary with Atomic Number to Symbol mapping
        """
        Out = dict()
        if not self.orientation: # In case it was parched
            raise RuntimeError('Attempting to access a non-existing orientation')
        atoms = [(atom[0],atom[1]) for atom in self.orientation ]
        _iter = self._text_iterbyline()
        for line in _iter:
            if 'Distance matrix' in line:
                break
        try:
            _ = next(_iter)
        except StopIteration:
            raise StopIteration('Attempting to get a mapping from a l202 without Distance matrix')
        for AtId,AtNum in sorted(atoms,key=lambda x:float(x[0])):
            line = next(_iter)
            Aux = line.strip().split()[:2]
            _, AtSym = Aux
            Out[AtNum] = AtSym
        return Out

#@RegisterLinkJob
class Link301(LinkJob):
    """
    Parser for the output of l301.exe, each execution displays
    information about solvation model. Holds information about
    Cavity Surface Area, and Cavity Volume. In development.
    """
    _token = 301

@RegisterLinkJob
class Link502(LinkJob):
    """
    Representation and parser for the output of l502.exe. Corresponds to a
    iterative calculation of the SCF cycles.

    Parameters
    ----------
    text : str
        text that corresponds to the output of the l502.exe
    asEmpty : bool
        Flag to not parse and store the information of the text.
        (defaults to False)

    Attributes
    ----------
    energy : float
        Final Potential Energy of the SCF cycles.
    spin : tuple[float,float]
        S**2 (before, after) annihiliation
    """
    __slots__ = ('energy','spin')

    _token = 502

    re_EDone = re.compile(r'SCF\sDone\:\s*E\(.*\).*=\s*(\-?[0-9]*\.[0-9]+)')
    re_spin = re.compile('S\*\*2 before annihilation\s*(-?[0-9]*\.[0-9]*),\s*after\s*(-?[0-9]*\.[0-9]*)')

    def __init__(self,text,asEmpty=False):
        self.energy = None
        if asEmpty:
            super().__init__('',502)
        else:
            super().__init__(text,502)
            self._locate_energy()
            self._locate_spin()

    @Populates('energy')
    @SilentFail
    def _locate_energy(self):
        """
        Uses regex expressions compiled as class attributes to find the Energy
        rported as 'SCF Done:' and reads the potential energy.
        """
        cls = self.__class__
        Match = cls.re_EDone.findall(self.text)
        if Match:
            self.energy = float(Match[0])

    @Populates('spin')
    @SilentFail
    def _locate_spin(self):
        """
        Uses regex expressions compiled as class attributes to find the Energy
        rported as 'SCF Done:' and reads the potential energy.
        """
        cls = self.__class__
        Match = cls.re_spin.findall(self.text)
        if Match:
            before,after = Match[0]
            self.spin = float(before),float(after)
    

@RegisterLinkJob
class Link508(Link502):
    """
    Representation and parser for the output of l508.exe. Corresponds to a
    iterative calculation of the SCF using a quadratic convergence algorithm.

    Parameters
    ----------
    text : str
        text that corresponds to the output of the l508.exe
    asEmpty : bool
        Flag to not parse and store the information of the text.
        (defaults to False)

    Attributes
    ----------
    energy : float
        Final Potential Energy of the SCF cycles.
    spin : tuple[float,float]
        S**2 (before, after) annihiliation
    """
    __slots__ = ('energy','spin')

    _token = 508

    def __init__(self,text,asEmpty=False):
        super().__init__(text,asEmpty)
        self.number = 508

@RegisterLinkJob
class Link601(LinkJob):
    """
    Representation and parser for the output of l601.exe, 
    each does a population analysis using SCF densities. Holds information 
    about Mulliken Charges, {Dipole, Quadrupole, Traceless Quadrupole, Octapole 
    and Hexadecapole} moments, spin densities... etc. 
    (Currently only Mulliken Charges implemented)

    Parameters
    ----------
    text : str
        text that corresponds to the output of the l508.exe
    asEmpty : bool
        Flag to not parse and store the information of the text.
        (defaults to False)

    Attributes
    ----------
    mulliken_heavy : list
        list of mulliken charges condensed to H atoms. None if not properly parsed.
    mulliken : list
        list of mulliken charges. None if not properly parsed.
    """
    re_MullikenAtoms =  r'(?:Mulliken charges( and spin densities)?.*\n.*\n)'# Start
    re_MullikenAtoms += r'([\s\S]*)' # Body
    re_MullikenAtoms += r'(?:\n.Sum of Mulliken charges.*\n)' # End
    re_MullikenAtoms = re.compile(re_MullikenAtoms)
    re_MullikenHeavy =  r'(?:Mulliken charges( and spin densities)? with hydrogens.*\n.*\n)'
    re_MullikenHeavy += r'([\s\S]*?)'
    re_MullikenHeavy += r'(?:\n.[a-zA-Z].*\n)'
    re_MullikenHeavy = re.compile(re_MullikenHeavy)
    
    _token = 601

    def __init__(self,text,asEmpty=False):
        self._Atom = namedtuple('Atom','number symbol charge spin')
        self.mulliken_heavy = None
        self.mulliken = None
        if asEmpty:
            super().__init__('',601)
        else:
            super().__init__(text,601)
            self._locate_MullikenHeavy()
            self._locate_MullikenAtoms()

    @Populates('mulliken_heavy')
    @SilentFail
    def _locate_MullikenHeavy(self):
        """
        Uses regex expressions compiled as class attributes to find the
        Mulliken charges and spin condensed to heavy atoms.
        """
        cls = self.__class__
        Match = cls.re_MullikenHeavy.findall(self.text)
        if Match:
            lines =  Match[0][1].split('\n')
            if Match[0][0]:
                Atom = lambda x: self._Atom(int(x[0]),x[1],float(x[2]),float(x[3]))
            else:
                Atom = lambda x: self._Atom(int(x[0]),x[1],float(x[2]),spin=None)
            self.mulliken_heavy = [Atom(line.strip().split())
                                    for line in lines if line.strip()]

    @Populates('mulliken')
    @SilentFail
    def _locate_MullikenAtoms(self):
        """
        Uses regex expressions compiled as class attributes to find the
        atomic Mulliken charges.
        """
        cls = self.__class__
        Match = cls.re_MullikenAtoms.findall(self.text)
        if Match:
            lines =  Match[0][1].split('\n')
            if Match[0][0]:
                Atom = lambda x: self._Atom(int(x[0]),x[1],float(x[2]),float(x[3]))
            else:
                Atom = lambda x: self._Atom(int(x[0]),x[1],float(x[2]),spin=None)
            
            self.mulliken = [Atom(line.strip().split())
                             for line in lines if line.strip()]

@RegisterLinkJob
class Link716(LinkJob):
    """
    Representation and parser for the output of l716.exe. Displays
    information about Forces and Frequency Calculations.

    Parameters
    ----------
    text : str
        text that corresponds to the output of the l716.exe
    asEmpty : bool
        Flag to not parse and store the information of the text.
        (defaults to False)

    Attributes
    ----------
    dipole : list
    units : list
    zeropoint : list
    thermal_energy : list
    enthalpy : list
    gibbs : list
    EContribs : list
        List of namedtuples with fields Name,Thermal,CV,S .
    IRSpectrum : str
        Text that corresponds to the IRSpectrum
    mode : str
        Either 'Forces', 'Freq' or 'Other'.
    frequencies : list
    freq_displacements : list
    """
    __slots__ = ('dipole', 'units','frequencies','EContribs', 'IRSpectrum',
                  'zeropoint','thermal_energy', 'enthalpy', 'gibbs', 'mode',
                  'freq_displacements')

    _token = 716

    re_Thermo = r'(Zero-point\scorrection=[\s\S]*)(?:\n\s*E\s\(Thermal\))'
    re_Thermo = re.compile(re_Thermo)
    re_EContrib =  r'(?:.*E\s\(Thermal\)\s*CV\s*S.*\n)' # Header
    re_EContrib += r'([\s\S]*\n)'                    # Body
    re_EContrib += r'(?:.*Q.*Log10\(Q\).*Ln\(Q\).*)' # Header of the next Table
    re_EContrib = re.compile(re_EContrib)
    re_IRSpectrum = re.compile(r'(.*IR\sSpectrum[\s\S]*?\n)(?:\s?\-\-+)')
    re_Frequencies =  r'(?:.*Frequencies\s\-\-)(.*)\n'
    re_Frequencies += r'(?:.*Red\.\smasses\s\-\-)(.*)\n'
    re_Frequencies += r'(?:.*Frc\sconsts\s*\-\-)(.*)\n'
    re_Frequencies += r'(?:.*IR\sInten\s*\-\-)(.*)'
    re_Frequencies = re.compile(re_Frequencies)
    re_freq_text = re.compile(r'  Atom[\s\S]*\n\n')
    re_freq_displacements = r'^ {3,5}'
    re_freq_displacements += r'[0-9]{1,3}\s*'*2
    re_freq_displacements += r'-?[0-9]*\.[0-9]*\s*'*3
    re_freq_displacements += r'[^\n]*\n'
    re_freq_displacements = re.compile(re_freq_displacements,re.MULTILINE)
    re_dipole = r'(?:Dipole.*\=)'+r'(.?[0-9]\.[0-9]*D[\-|\+][0-9]{2,})'*3
    re_dipole = re.compile(re_dipole)

    _Frequency = namedtuple('Frequency','freq redmass forcek IRint')
    _Displacement = namedtuple('FreqDisplacements','atomids atoms xyz')
    _EContrib = namedtuple('EContrib','Name Thermal CV S')

    def __init__(self,text,asEmpty=False):
        self.dipole = []
        self.units = [] # The first one corresponds to the non tabulated
        self.zeropoint = []
        self.thermal_energy = []
        self.enthalpy = []
        self.gibbs = []
        self.EContribs = []
        self.frequencies = []
        self.freq_displacements = []
        self.IRSpectrum = ''
        self.mode = ''
        if asEmpty:
            super().__init__('',716)
        else:
            super().__init__(text,716)
            if 'Forces' in self.text:
                self.mode = 'Forces'
            else:
                self.mode = 'Other'
            self._locate_dipole()
            self._locate_thermochemistry()
            self._locate_IRSpectrum()
            self._locate_frequencies()
            self._locate_freq_displacements()

    @Populates('dipole')
    @SilentFail
    def _locate_dipole(self):
        """
        Uses regex expressions compiled as class attributes to find the Dipole
        Vector.
        """
        cls = self.__class__
        Match = cls.re_dipole.findall(self.text)
        if Match:
            self.dipole = [float(i.replace('D','E')) for i in Match[0]]

    @Populates('frequencies')
    @SilentFail
    def _locate_frequencies(self):
        """
        Uses regex expressions compiled as class attributes to find the
        frequencies table.
        """
        cls = self.__class__
        Frequency = cls._Frequency
        Match = cls.re_Frequencies.findall(self.text)
        frequencies = []
        if Match:
            for i in Match:
                Aux = []
                for line in i:
                    Aux.append(tuple(map(float,line.strip().split())))
                for items in zip(*Aux):
                    frequencies.append(Frequency(*items))
            self.frequencies = frequencies
            self.mode = 'Freq'

    @Populates('freq_displacements')
    @SilentFail
    def _locate_freq_displacements(self): 
        """
        Uses regex expressions compiled as class attributes to find the
        cartesian displacements of the frequencies.
        """
        cls = self.__class__
        Displacement = cls._Displacement
        match = cls.re_freq_text.findall(self.text)
        if match:
            displacements = []
            for text in match[0].split('Atom'):
                if not text.strip(): 
                    continue
                lines = cls.re_freq_displacements.findall(text)
                atomlines = []
                for line in lines:
                    atomline = []
                    if not line.strip():
                        continue
                    atid,atnum,xyz = line.strip().split(maxsplit=2)
                    atomline.append(int(atid))
                    atomline.append(int(atnum))
                    xyz = [float(i) for i in xyz.split()]
                    for x,y,z in zip(xyz[0::3],xyz[1::3],xyz[2::3]): 
                        atomline.append((x,y,z))
                    atomlines.append(atomline)
                atomids,atomnums,*matrices = zip(*atomlines)
                for matrix in matrices: 
                    displacements.append(Displacement(atomids,atomnums,matrix))
            self.freq_displacements = displacements

    @Populates('zeropoint','thermal_energy','enthalpy','gibbs')
    @SilentFail
    def _locate_thermochemistry(self):
        """
        Uses regex expressions compiled as class attributes to find the
        thermochemistry computed properties.
        """
        cls = self.__class__
        EContrib = cls._EContrib
        Match = cls.re_Thermo.findall(self.text)
        if Match:
            lines = [line.strip() for line in Match[0].split('\n') if line.strip()]
            Properties = cycle([self.zeropoint,
                                self.thermal_energy,
                                self.enthalpy,
                                self.gibbs])
            for line,Property in zip(lines,Properties):
                Aux = line.split("=")[-1].strip().split()
                Property.append(float(Aux[0]))
                if len(Aux) == 2:
                    self.units.append(Aux[-1][1:-1])
                Property.append(float(Aux[0]))
        Match = cls.re_EContrib.findall(self.text)
        if Match:
            lines = [line.strip() for line in Match[0].split('\n') if line.strip()]
            Header = lines.pop(0)
            self.units.extend(Header.strip().split())
            for line in lines:
                Name, E, CV, S = line.strip().rsplit(maxsplit=3)
                item = EContrib(Name,float(E),float(CV),float(S))
                self.EContribs.append(item)

    @Populates('IRSpectrum')
    @SilentFail
    def _locate_IRSpectrum(self):
        """
        Uses regex expressions compiled as class attributes to find the IR
        spectrum.
        """
        cls = self.__class__
        Match = cls.re_IRSpectrum.findall(self.text)
        if Match:
            self.IRSpectrum = Match[0]

@RegisterLinkJob
class Link804(LinkJob):
    """
    Representation and parser for the output of l804.exe. Holds the
    information about the MP2 energy calculation

    Parameters
    ----------
    text : str
        text that corresponds to the output of the l804.exe
    asEmpty : bool
        Flag to not parse and store the information of the text.
        (defaults to False)

    Attributes
    ----------
    MP2 : float
        Potential energy with MP2 method
    SpinComponents : list
        List of namedtuples with the fields Name, T, E.
    """
    __slots__ = ('MP2', 'SpinComponents')

    _token = 804

    re_MP2 = re.compile(r'EUMP2\s=\s*(\-?[0-9]+\.[0-9]*D[\-|\+][0-9]{2,3})')

    _SpinComponent = namedtuple('SpinComponent','Name T E')

    def __init__(self,text,asEmpty=False):
        self.MP2 = None
        self.SpinComponents = [] # The first one corresponds to the non tabulated
        if asEmpty:
            super().__init__('',804)
        else:
            super().__init__(text,804)
            self._locate_SpinComponents()
            self._locate_MP2()

    def _locate_MP2(self):
        """ Looks for the keyword 'EUMP2=' and reads the MP2 Potential Energy """
        cls = self.__class__
        Match = cls.re_MP2.findall(self.text)
        if Match:
            self.MP2 = float(Match[0].replace('D','E'))
    def _locate_SpinComponents(self):
        """ Looks for the keywords 'Spin Components' and reads
        the Energies of the spin Components """
        Iterator = self._text_iterbyline()
        cls = self.__class__
        SpinComponent = cls._SpinComponent
        for line in Iterator:
            if 'Spin components' in line:
                break
        else: # fail silently
            return
        for _ in range(3):
            line = next(Iterator)
            # Standarize number of columns
            line = line.replace('=',' ').replace('E2',' ').replace('T2',' ')
            # Change Number Notation
            line = line.replace('D','E')
            # Store spin Component
            Name,T,E = line.split()
            self.SpinComponents.append(SpinComponent(Name,float(T),float(E)))
    def get_SCScorr(self):
        """ Calculates and returns the MP2(SCS) potential energy"""
        aa = self.SpinComponents[0].E
        ab = self.SpinComponents[1].E
        bb = self.SpinComponents[2].E
        return (aa+bb)/3.0 +(6.0*ab)/5.0

@RegisterLinkJob
class Link913(LinkJob):
    """Representation and parser for the output of l913.exe. Holds information
    about CCSD(T) calculations including MP4 Energy.

    Parameters
    ----------
    text : str
        text that corresponds to the output of the l913.exe
    asEmpty : bool
        Flag to not parse and store the information of the text.
        (defaults to False)

    Attributes
    ----------
    MP4 : float
        MP4 potential energy
    CCSDT : float
        CCSD(T) potential energy
    """
    # Current Implementation only stores the MP4(SDQ) Energy, other energies
    # Should be implemented using a different name or reestructure all code

    __slots__ = ('MP4', 'CCSDT')

    _token = 913

    re_MP4 = r'[A-Z]*4\(S?DQ\)=\s*\-?[0-9]+\.[0-9]*D[\+|\-][0-9]{2,3}'
    re_MP4 = re.compile(re_MP4)
    re_ccsdt = re.compile(r'CCSD\(T\)=\s?(\-?[0-9]+\.[0-9]*D[\+|\-][0-9]{2,3})')

    def __init__(self,text,asEmpty=False):
        self.MP4 = None
        self.CCSDT = None
        if asEmpty:
            super().__init__('',913)
        else:
            super().__init__(text,913)
            self._locate_MP4()
            self._locate_CCSDT()
    def _locate_MP4(self):
        """ Looks for the First iteration and then looks for the keyword
        'UMP4(SDQ)' and reads the MP4 Potential Energy """
        cls = self.__class__
        Match = cls.re_MP4.findall(self.text)
        if Match and Match[-1].startswith('UMP4(SDQ)'):
            self.MP4 = float(Match[-1].split("=")[1].replace('D','E'))
    def _locate_CCSDT(self):
        """ Looks for the keywords 'Time for triples' which appears after the
        calculation has converged and reads the CCSD(T) Potential Energy """
        cls = self.__class__
        Match = cls.re_ccsdt.findall(self.text)
        if Match:
            self.CCSDT = float(Match[0].replace('D','E'))

@RegisterLinkJob
class Link914(LinkJob):
    """
    Representation and parser for the output of l914.exe. Holds information
    about TDDFT calculations excited states.

    Parameters
    ----------
    text : str
        text that corresponds to the output of the l913.exe
    asEmpty : bool
        Flag to not parse and store the information of the text.
        (defaults to False)

    Attributes
    ----------
    excitedstates : list
        list of excited states
    """

    __slots__ = ('excitedstates',)

    _token = 914

    _ExcitedState = namedtuple('ExcitedState',
                            'number energy wavelen OStrenght s2 transitions')
    _TransitionData = namedtuple('TransitionData',
                            'donor acceptor contribution isreversed')

    re_ExcitedState =  r'Excited\s*State\s*([0-9]{1,3})\:' # Excited State number
    re_ExcitedState += r'\s*(\S*)'     # Match Pointgroup
    re_ExcitedState += r'\s*(\S*)\seV' # Match Energy
    re_ExcitedState += r'\s*(\S*)\snm' # Match wavelength
    re_ExcitedState += r'\s*f\=(\S*)'  # Match Oscilator Strength
    re_ExcitedState += r'\s*\S*\=(\S*).*' # Match <S**2>
    re_ExcitedState = re.compile(re_ExcitedState)

    re_transition  = r'([0-9]{1,4}[A|B]?)\s*' # match first transition id
    re_transition += r'(\-\>|\<\-)\s*'        # Match "donation" direction
    re_transition += r'([0-9]{1,4}[A|B]?)\s*' # match second transition id
    re_transition += r'(\-?[0-9]*\.[0-9]*)'   # Match contribution
    re_transition = re.compile(re_transition)

    def __init__(self,text,asEmpty=False):
        self.excitedstates = []
        if asEmpty:
            super().__init__('',914)
        else:
            super().__init__(text,914)
            self._locate_ExcitedStates()

    @Populates('excitedstates')
    @SilentFail
    def _locate_ExcitedStates(self):
        """
        Does the logic of finding and slicing the text of the link into the
        different excited states
        """
        ExcitedState = self._ExcitedState
        excitedstates = []
        matches = list(self.re_ExcitedState.finditer(self.text))
        # Fail silently
        if not matches:
            return
        for current,new in zip(matches,matches[1:]):
            start,stop = current.end(),new.start()
            text = self.text[start:stop]
            number,_,energy,wavelen,OStrenght,s2 = current.groups()
            transitions = self._extract_transitions(text)
            excitedstates.append(ExcitedState(int(number),
                                              float(energy),
                                              float(wavelen),
                                              float(OStrenght),
                                              float(s2),
                                              transitions))
        self.excitedstates = excitedstates

    def _extract_transitions(self,text):
        """
        Given the slice of output text that corresponds to the transition
        contributions of a ExcitedState, returns each contribution as a
        TransitionData object.

        Parameters
        ----------
        text : str
            The slice of text that corresponds to the transition
            contributions of an ExcitedState.

        Returns
        -------
        transitions
            list of TransitionData instances.

        """
        TransitionData = self._TransitionData
        transitions_matches = self.re_transition.findall(text)
        transitions = []
        if transitions_matches:
            for donor,direction,acceptor,contribution in transitions_matches:
                isreversed = False
                if direction == '<-':
                    isreversed = True
                    donor,acceptor = acceptor,donor
                transition = TransitionData(donor,
                                            acceptor,
                                            float(contribution),
                                            isreversed)
                transitions.append(transition)
        return transitions

    def print_excitedstates(self,*ESnumbers,show_transitions=False):
        """
        Prints in the console the parameters in of the excited states selected 
        in ascending ordinal order.

        Parameters
        ----------

        *ESnumbers : int
            An undefined number of integrers. Only the excited states that are
            within the provided set will be displayed.
        show_transitions : bool, optional
            Allows control over the display of the transitions, by 
            default False
        """
        for ES in self.excitedstates:
            if ES.number not in ESnumbers:
                continue
            number,energy,wavelen,OStrenght,s2,_ = ES
            print(f'ExcitedState(number={number} energy={energy} ' \
                  f'wavelen={wavelen} OStrenght={OStrenght} ' \
                  f's2={s2} transitions=[...])')
            if show_transitions:
                for transition in ES.transitions:
                    print(f'\t{transition.donor} -> {transition.acceptor}'\
                          f'\t {transition.contribution}')

#@RegisterLinkJob
class Link9999(LinkJob):
    """
    parser for the output of l9999.exe. This Link has 2 termination modes:
    - normal, if so it will display a resume string
    - error, if so it will display or not information about the error
    Depending on the InternalJob, it might display a final structure or not,
    (before printing the resume string). In development.
    """
    _token = 9999
