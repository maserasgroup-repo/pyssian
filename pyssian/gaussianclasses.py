"""
One of the two core libraries of pyssian. Contains the Classes that represent
Gaussian Files (input and output).
"""
import io
import re
from itertools import chain

from .chemistryutils import is_method, is_basis
from .linkjobparsers import LinkJob, GeneralLinkJob

# Pre-Initialized dictionary for the GaussianOutFile.Parse
Available_Linkjobs = {i:GeneralLinkJob for i in range(1,10000)}
for key in LinkJob.Register.keys():
    Available_Linkjobs[key] = LinkJob.Register[key]


class GaussianOutFile(object):
    """Gaussian 09/16 '.log' file parent class, if any special type of calculation
    requires different processing it should be a subclass of this one. Accepts
    a context manager usage similar to 'with open(file) as F:...'
    *For a Gaussian .log file to be parsable it requires that its corresponding 
    input has the additional printout enabled (#p)* 

    Parameters
    ----------
    file : io.TextIOBase or str
        File instance (Result of open(filename,'r')) or valid filename.
    parselist : list
        List of integrers that represent which types of Links to parse
        (the default is None).

    Attributes
    ----------
    InternalJobs
        List of InternalJobs done by gaussian i.e an gaussian calculation with
        the opt freq keywords will run first an InternalJob for the Optimization
        and after an InternalJob for the Frequency calculation.
    """
    _interblock = -1        # interblock token
    _EOF = -9999            # EOF token

    def __init__(self,file,parselist=None):
        cls = self.__class__
        self.InternalJobs = [InternalJob(),]
        if isinstance(file,io.TextIOBase):
            self._file = file
        else:
            self._file = open(file,'r')
        if parselist is None:
            parselist = []
        # Access the dictionary that holds the constructors for each LinkJob
        self._SetParsers(parselist,cls._interblock)
        # Initialize the generators/coroutines
        self._BlockFetcher = self.BlockFetcher(cls._EOF,cls._interblock)
        _ = next(self._BlockFetcher)
        self._BlockHandler = self.BlockHandler()
        _ = next(self._BlockHandler)
    def __repr__(self):
        cls = type(self).__name__
        file = self._file.name.split('/')[-1]
        size = len(self)
        return f'<{cls}({file})> with {size} InternalJobs'
    def __str__(self):
        cls = type(self).__name__
        file = self._file.name.split('/')[-1]
        repr = f'<{cls}({file})>\n'
        indent = '    '
        for InternalJob in self:
            repr += indent + f'{InternalJob} type <{InternalJob.type}>\n'
            for Link in InternalJob:
                 repr += indent*2 + f'{Link}\n'
        return repr
    def __len__(self):
        return len(self.InternalJobs)
    def __getitem__(self,index):
        return self.InternalJobs[index]

    def __enter__(self):
        ''' Wrapper to have similar behaviour to "_io.TextIOWrapper" '''
        return self
    def __exit__(self, exc_type, exc_value, traceback):
        ''' Wrapper to have similar behaviour to "_io.TextIOWrapper" '''
        return self._file.__exit__(exc_type, exc_value, traceback)

    def _SetParsers(self,ParseList,interblock=-1):
        Parsers = Available_Linkjobs.copy()
        assert interblock not in Parsers
        if ParseList:
            if ParseList[0] == -1: # flag to parse all as empty
                for key in Parsers:
                    Parsers[key] = Parsers[key].as_empty
            else: # Set up the appropiate parsers
                for key in Parsers:
                    if key not in ParseList:
                        Parsers[key] = Parsers[key].as_empty
        else: # Parse all normally
            pass
        Parsers[interblock] = GeneralLinkJob.as_empty
        self._Parsers = Parsers
    def print_file_structure(self):
        """Display the structure of links and internal jobs of the file."""
        indent = "  "
        Result = f"{self.__repr__():}\n"
        for intjob in self:
            Result += indent + f"{intjob.__repr__():}\n"
            for link in intjob:
                Result += indent*2 + f"{link.__repr__():}\n"
        print(Result)
    def read(self):
        """Alias of update for consistency with GaussianInFile class"""
        self.update()
    def close(self):
        """Alias to file.close for consistency with the io.TextIOBase class"""
        self._file.close()
    def update(self,clean=True,FinalPrint=False):
        """Tries to fetch new data. If it exists it parses it appropiately
        otherwise it fails silently.

        Parameters
        ----------
        clean : Bool
            If True removes all the EmptyLinkJobs found (the default is True).
        FinalPrint : Bool
            If True after a normal execution has finished it will print in the
            console a message to notify the user (the default is False).
        """
        cls = self.__class__
        BlockFetcher = self._BlockFetcher
        BlockHandler = self._BlockHandler
        BlockType, Block = next(BlockFetcher)
        while BlockType != cls._EOF:
            BlockHandler.send((BlockType, Block))
            BlockType, Block = next(BlockFetcher)
        if self.InternalJobs[0].number is None:
            self.InternalJobs[0].guess_info()
        if clean:
            self.clean()
        if FinalPrint:
            print(f"{self:!r} UPDATED")
    def clean(self):
        """Removes per each InternalJob stored all the EmptyLinkJobs."""
        for InternalJob in self.InternalJobs:
            InternalJob.clean()
    def get_links(self,*LinkIds):
        """Wrapper Method to get a list of Links with certain Ids across
        the different Internal Jobs.

        Parameters
        ----------
        *LinkIds : int
            Integrers that correspond to the type of links to be return.

        Returns
        -------
        list
        """
        LinkLists = [IntJob.get_links(*LinkIds) for IntJob in self.InternalJobs]
        return list(chain(*LinkLists))
    # Generators and Coroutines for File Parsing
    def Reader(self,file):
        """ Generator for line by line reading without stopiteration """
        while True:
            Pos_old = file.tell()
            line = file.readline()
            Pos_new = file.tell()
            ReachedEOF = Pos_old == Pos_new
            yield ReachedEOF,line
    def BlockFetcher(self,EOF=-9999,interblock=-1):
        """
        Generator that yields the text sliced in blocks and their type.

        A block is an iterable of strings and its type refers to a token that
        can be recognized by the ._Parsers variable, something in betweem
        Link Blocks (interblock=-1) or no end was found (EOF=-9999)
        """
        # Regex generation
        re_enter = re.compile(r'(?:Enter.*l)([0-9]{1,4})(?:\.exe)')
        re_exit = re.compile(r'(?:Leave\s*Link\s*)([0-9]{1,4})')
        re_termination = re.compile(r'\s?([a-zA-Z]*)\stermination')

        # Initialize the Reader generator
        Reader = self.Reader(self._file)

        yield 'Initialization done'
        # If more than 1 EndKeyws then BlockType Assesment has to be modified
        while True:
            start = False
            Block = []
            # Ask the Reader until a "start line" is found
            while not start:
                ReachedEOF,line = next(Reader)
                if ReachedEOF:
                    yield EOF, ''
                else:
                    start = re_enter.findall(line)
                    if not start:
                        Block.append(line)
                    else: # Store the number of the Link
                        number = int(start[0])
            # When found yield it as a "InterBlock" and prepare Block
            if Block:
                yield interblock, ''.join(Block)
                Block = [line,]
            else:
                Block.append(line)
            # Now that the start of the Link has been found, accumulate lines
            ## until the end or termination line is found
            end = False
            while not end:
                ReachedEOF,line = next(Reader)
                if ReachedEOF:
                    Target = Block[-10:] + [line,]
                    terminated = re_termination.findall(''.join(Target))
                    if terminated:
                        Block.append(line)
                        #if terminated[0] == 'Normal':
                        #    other = str(BlockTypes[Termination_key])
                        #elif terminated[0] == 'Error':
                        #    other = str(BlockTypes[Error_key])
                        break
                    else:
                        yield EOF, ''
                else:
                    end = re_exit.findall(line)
                Block.append(line)
            # when end found, do return type token and yield the block
            yield number, ''.join(Block)
    def BlockHandler(self):
        """ Coroutine. Receives a block, chooses the parser and parses it """
        # Initialization
        Parsers = self._Parsers
        CurrentJob = self.InternalJobs[-1]
        BlockType, Block = yield 'Initialization done'
        while True:
            #Parser = Parsers.get(BlockType,Parsers[ignore])
            Parser = Parsers[BlockType]
            Link =  Parser(Block)
            if Link.number == 1:
                new_InternalJob = Link.info.new_InternalJob
            else:
                new_InternalJob = False
            if new_InternalJob:
                New = InternalJob()
                self.InternalJobs.append(New)
                CurrentJob = self.InternalJobs[-1]
                CurrentJob.append(Link)
                CurrentJob.guess_info()
            else:
                CurrentJob.append(Link)
            BlockType, Block = yield

class InternalJob(object):
    """Gaussian 09/16 InternalJob parent class, if any special type of Job
    requires different parsing it should be a subclass of this one.

    Parameters
    ----------
    number : int
        ordinal number of the InternalJob (the default is None).

    Attributes
    ----------
    type
        string identifier for the job.
    Links
        List of the different Links that belong to the InternalJob.
    number

    """

    def __init__(self,number=None):
        self.number = number
        self.type = None
        self.Links = []
    def __repr__(self):
        cls = type(self).__name__
        if self.number is None:
            return f'<{cls} Created but Empty>'
        else:
            return f'<{cls} {self.number}>'
    def __str__(self):
        return f'Internal Job {self.number}: {self.type}'
    def __getitem__(self,index):
        return self.Links[index]
    def __len__(self):
        return len(self.Links)

    def append(self,Link):
        # Restrict to Link objects
        if not isinstance(Link, LinkJob):
            raise TypeError(f'{Link:!r} is not of class {LinkJob:!r}')
        self.Links.append(Link)
    def guess_info(self):
        """ Guesses the number and type attributes of itself using the stored
        Links."""
        if self.Links:
            Links = (Link for Link in self.Links if Link.number == 1)
            try:
                StarterLink = next(Links)
                info = StarterLink.info
            except AttributeError:
                pass
            except StopIteration:
                pass
            else:
                self.number = info.number
                self.type = info.type
    def clean(self):
        """Removes all the Empty Link instances within Links."""
        Indices2Remove = []
        for i, Link in enumerate(self.Links):
            if not Link.text:
                Indices2Remove.append(i)
        for index in reversed(Indices2Remove):
            _ = self.Links.pop(index)
    def get_links(self,*LinkIds):
        """Wrapper Method to get a list of Links with certain Ids.

        Parameters
        ----------
        *LinkIds : int
            Integrers that correspond to the type of links to be return.

        Returns
        -------
        list
            List of Link Objects ordered by appearance in the file and filtered
            by Link Number.
        """
        return [Link for Link in self.Links if Link.number in LinkIds]

class GaussianInFile(object):
    """
    Gaussian 09/16 .in file parent class, if any special type of input
    requires different processing it should be a subclass of this one.

    Parameters
    ----------
    file : io.TextIOBase or str
        File instance (Result of open(filename,'r')) or valid filename.

    Attributes
    ----------
    preprocessing : dict
        Dictionary in which each key corresponds to a certain Link0 keyword
    commandline : dict
        Dictionary that contains the information of how the calculation
        will be carried out.
    title : str
        title of the calculation.
    method : str
        If it cannot be recognized in the command line it will be empty.
    basis : str
        If it cannot be recognized in the command line it will be empty.
    spin : int
    charge : int
    geometry : str-ish
        It should be able to write the text block of an input file upon calling
        str(geometry)
    tail : list
        List of str in which each should be separated from the others by a
        single blank line in the input file.
    structure : str
        A string holding the structure of the input file. Used to write new
        Input files.
    nprocs : int
        property to easily access and change the preprocessing['nprocshared'] value
    mem : int
        property to easily access and change the preprocessing['mem'] value
    """
    def __init__(self,file):
        # Do Something
        if isinstance(file,io.TextIOBase):
            self._file = file
        else:
            self._file = open(file,'a+')
            if self._file.tell() != 0:
                self._file.seek(0)
        self._txt = ''
        self.preprocessing = dict() # In the G16 Manual "Link 0 Commands"
        self.commandline = dict() # In the G16 Manual "Route Section"
        self.title = ''
        self._method = ''
        self._basis = ''
        self.spin = 1
        self.charge = 0
        self.geometry = '' # In the G16 Manual "Molecule Specification"
        self.tail = [] # In the G16 Manual "Optional additional sections"
        self.structure = '{preprocessing}\n{commandline}\n\n'
        self.structure += '{title}\n\n'
        self.structure += '{charge} {spin}\n{geometry}\n\n'
        self.structure += '{tail}\n\n'
    def __repr__(self):
        cls = type(self).__name__
        file = self._file.name.split("/")[-1]
        size = len(self)
        return f'<{cls}({file})>'
    def __str__(self):
        # str repr of the preprocessing
        preprocessing = []
        for key,val in self.preprocessing.items():
            if val:
                Aux = f'%{key}={val}'
            else:
                Aux = f'%{key}'
            preprocessing.append(Aux)
        # str repr of the commandline
        commandline = ['#p',]
        for key,val in self.commandline.items():
            if val and (len(val) == 1):
                Aux = f"{key}={','.join(val)}"
            elif val:
                Aux = f"{key}=({','.join(val)})"
            else:
                Aux = f"{key}"
            commandline.append(Aux)
        # Prepare to format as str
        kwargs = dict( preprocessing='\n'.join(preprocessing),
                       commandline=' '.join(commandline),
                       title=self.title,
                       charge=self.charge,
                       spin=self.spin,
                       geometry=str(self.geometry),
                       tail='\n\n'.join(self.tail))
        return self.structure.format(**kwargs)
    def __len__(self):
        return len(str(self))

    def __enter__(self):
        ''' Wrapper to have similar behaviour to "_io.TextIOWrapper" '''
        return self
    def __exit__(self, exc_type, exc_value, traceback):
        ''' Wrapper to have similar behaviour to "_io.TextIOWrapper" '''
        return self._file.__exit__(exc_type, exc_value, traceback)

    @property
    def method(self):
        return self._method
    @method.setter
    def method(self,other):
        self.change_method(other)

    @property
    def basis(self):
        return self._basis
    @basis.setter
    def basis(self,other):
        self.change_basis(other)

    @property
    def nprocs(self):
        return self.preprocessing.get('nprocshared',None)
    @nprocs.setter
    def nprocs(self,other):
        self.preprocessing['nprocshared'] = other

    @property
    def mem(self):
        return self.preprocessing.get('mem',None)
    @mem.setter
    def mem(self,other):
        self.preprocessing['mem'] = other

    def read(self):
        """
        Reads the file and populates the appropiate attributes.
        """
        txt = [line.strip() for line in self._file]
        if not txt: 
            raise EOFError('Attempting to read an empty or non-existent file')
        if txt[-1]: # If the file ends without a blank line add it
            txt.append('')
        if txt[-2]: # If the file ends without two blank lines add one
            txt.append('')
        self._txt = '\n'.join(txt)
        bins = [i for i,line in enumerate(txt) if not line]
        # Ensure that the if the title is empty, the bins are not including it
        bins = [i for i in bins if not set((i-1,i,i+1)).issubset(set(bins))]
        stop = bins[0]
        header = iter(txt[:stop])
        preprocessing = []
        for line in header:
            if line.startswith("%"):
                preprocessing.append(line.lstrip("%"))
            elif line.startswith("#"):
                break
        self.parse_preprocessing(preprocessing)

        # Read the command line assuming that the keywords cannot be in a
        # 'chopped in half' version between lines
        commandline = [line.split(),]
        for line in header:
            commandline.append(line.split())
        self.parse_commandline(commandline)

        # Read the Title Section
        start = bins[0]+1
        stop = bins[1]
        title = [line for line in txt[start:stop]]
        self.title = '\n'.join(title)

        # Read charge and spin
        charge,spin = txt[stop+1].split()
        self.charge,self.spin  = int(charge), int(spin)

        # Now we read the geometry
        start = stop+1
        stop = bins[2]
        geometry = [line for line in txt[start+1:stop]]
        self.parse_geometry(geometry)

        # Now we read the Tail
        tail = []
        if len(txt) > stop+1: # if it exists
            tail = [line for line in txt[stop:]]
            self.parse_tail(tail)
    def close(self):
        """Alias to file.close for consistency with the io.TextIOBase class"""
        self._file.close()
    def write(self,filepath=None):
        """
        Writes the File object to a File. If a filepath is provided it will
        write to that filepath otherwise it will attempt to write to the path
        provided in the initialization.

        Parameters
        ----------
        filepath : str
            A valid filepath.
        """
        self._txt = str(self)
        if filepath is None:
            # Write to self._file
            self._file.write(self._txt)
        else:
            # open the file write and close the file
            with open(filepath,'w') as F:
                F.write(self._txt)

    # Helper Functions for the read function to encapsulate different behaviours
    def parse_preprocessing(self,lines):
        """
        Parses the lines that contain the Link 0 keywords and transforms them
        into a dictionary representation.

        Parameters
        ----------
        lines : list
            list of strings previously stripped. Empty lines will be ignored.
        """
        #The logic of the specification of the preprocessing is below
        ## %kwd=something (most of the keywords follow this shape)
        ## %kwd           (Few keywords follow this shape)
        ## %kwd L{number} [something,or not] (2 kewyords follow this shape)
        # As initial design criteria I'm reducing the third case to the second
        for line in lines:
            Aux = line.split('=')
            if len(Aux) == 1:
                key,val = Aux[0],''
            elif len(Aux) == 2:
                key,val = Aux
            else:
                pass
            self.preprocessing[key] = val
    def parse_commandline(self,lines):
        """
        Parses the lines that contain the calculation commands keywords and
        transforms them into a dictionary representation.

        Parameters
        ----------
        lines : list
            list of strings previously stripped and split. Empty lines will be 
            ignored.
        """
        # the first line contains the "#p", remove it
        start = lines[0][1:]
        others = [i for i in chain(*lines[1:])]
        method_found = False
        basis_found = False
        for item in chain(start,others):
            if is_method(item) and not method_found:
                method_found = True
                self._method = item
                key,val = item, []
            elif is_basis(item) and not basis_found:
                basis_found = True
                self._basis = item
                key,val = item, []
            elif is_basis(item) or is_method(item):
                print('2 Basis or methods found \n')
                print(f'taking {item} as a normal keyword')
                key,val = item, []
            else:
                Aux = item.split('=',1)
                if len(Aux) == 1:
                    key,val = Aux[0],[]
                elif len(Aux) == 2:
                    key,val = Aux
                else:
                    pass # Should only enter with empty items
                if val and val.startswith('('):
                    val = val[1:-1].split(',')
                elif val:
                    val = [val,]
            previously_stored = key in self.commandline
            has_suboptions = bool(self.commandline.get(key,False))
            if not previously_stored:
                self.commandline[key] = val
            elif not has_suboptions:
                self.commandline[key].extend(val)
            else:
                Aux3 = set(self.commandline[key]) + set(val)
                self.commandline[key] = list(Aux3)
    def parse_geometry(self,lines):
        """Parses each line that contains 'Atom x y z' in an appropiate form
        and saves it to self.geometry

        Parameters
        ----------
        lines : list
            list of strings previously stripped. Should not contain empty lines
        """
        # This function is currently set up to only get the geometry as is.
        # In the future it should include the logic to transform between
        # coordinate specifications  (zmatrix,xyz,internal) to enforce a certain
        # geometry or geometry class
        self.geometry = '\n'.join(lines)
    def parse_tail(self,lines):
        """Chops the set of lines into different blocks of text using as
        reference the emptylines/blank lines

        Parameters
        ----------
        lines : list
            list of strings previously stripped.

        Returns
        -------
        type
            Description of returned object.

        Raises
        -------
        ExceptionName
            Why the exception is raised.

        """
        Aux = []
        self.tail = []
        for line in lines:
            if line:
                Aux.append(line)
            elif Aux:
                self.tail.append('\n'.join(Aux))
                Aux = []
            else:
                pass
        else:
            if Aux:
                self.tail.append('\n'.join(Aux))
                print('Parsed an input file withouth blank line ending')


    # Modifying functions
    def pop_chk(self,default=None):
        """
        Removes the chk from the file, returns 'default' if the chk was not
        included already
        """
        if default is not None:
            return self.preprocessing.pop('chk',default)
        else:
            return self.preprocessing.pop('chk')
    def add_chk(self,name=None):
        """
        Adds the chk to the file, with the specified name. If none is provided
        defaults to the file name ended in .chk
        """
        # Adds the chk to the file
        if name is None:
            try:
                name = self._file.name
            except AttributeError:
                name = self._file.split('/')[-1]
            name = name.rsplit('.')[0]+'.chk'
        else:
            if not name.endswith('.chk'):
                name = name +'.chk'
        self.preprocessing['chk'] = name
    def change_method(self,method):
        """Changes appropiately the method of the calculation. Running 
        self.method = method makes a call to this function.

        Parameters
        ----------
        method : str
            A string representation of a valid method

        Raises
        -------
        NotImplementedError
            If the method is not within the registered methods keywords

        """
        if not is_method(method):
            raise NotImplementedError(f'method {method} not implemented')
        key = self._method
        _ = self.commandline.pop(key,None) # Used to ensure deletion of the key
        self._method = method
        self.commandline[method] = ''
    def change_basis(self,basis):
        """Changes appropiately the basis of the calculation. Running 
        self.basis = basis makes a call to this function.

        Parameters
        ----------
        basis : str
            A string representation of a valid method if specified in the
            command line.

        Raises
        -------
        NotImplementedError
            If the basis is not within the registered basis keywords

        """
        if not is_basis(basis):
            raise NotImplementedError(f'basis {basis} not implemented')
        key = self._basis
        _ = self.commandline.pop(key,None) # Used to ensure deletion of the key
        self._basis = basis
        self.commandline[basis] = ''

# TODO: Implement a class to read and manipulate the basis functions in the tail
# class BasisTail(object), whose str function returns things as it should and
# that can have a linked input file object, so that modifying the basis of this
# object will modify the input file basis in certain cases.
