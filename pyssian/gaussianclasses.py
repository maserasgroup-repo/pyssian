"""
One of the two core libraries of pyssian. Contains the Classes that represent
Gaussian Files (input and output).
"""
import io
import re
from itertools import chain
import warnings

from .chemistryutils import is_method, is_basis
from .linkjobparsers import LinkJob, GeneralLinkJob, NPROCSHARED_ALIASES, MEMORY_ALIASES

# Pre-Initialized dictionary for the GaussianOutFile.Parse
available_linkjobs = {i:GeneralLinkJob for i in range(1,10000)}
for key in LinkJob.Register.keys():
    available_linkjobs[key] = LinkJob.Register[key]

EMPTYLINKFLAG = -1 

class GaussianOutFile(object):
    f"""
    Gaussian 09/16 '.log' file parent class, if any special type of calculation
    requires different processing it should be a subclass of this one. Accepts
    a context manager usage similar to 'with open(file) as F:...'
    *For a Gaussian .log file to be parsable it requires that its corresponding 
    input has the additional printout enabled (#p)* 

    Parameters
    ----------
    file : io.TextIOBase or str
        File instance (Result of open(filename,'r')) or valid filename.
    parselist : list[int]
        List of integrers that represent which types of Links to parse.
        If None or an empty list are provided it will attempt to parse
        all LinkJobs. If {EMPTYLINKFLAG} is in the list, all Links will be parsed 
        as_empty (This is usefull when only the structure of the file wants
        to be ascertained). Otherwise, only the specified Links will be parsed. 

    Attributes
    ----------
    jobs : list
        List of InternalJobs done by gaussian i.e an gaussian calculation with
        the opt freq keywords will run first an InternalJob for the Optimization
        and after an InternalJob for the Frequency calculation.
    """
    _interblock = -1        # interblock token
    _EOF = -9999            # EOF token

    def __init__(self,file,parselist:None|list[int]=None):
        cls = self.__class__

        self.jobs = [InternalJob(),]

        if isinstance(file,io.TextIOBase):
            self._file = file
        else:
            self._file = open(file,'r')
        
        if parselist is None:
            parselist = []
        
        # Access the dictionary that holds the constructors for each LinkJob
        self._set_parsers(parselist,cls._interblock)
        
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
        for job in self:
            repr += indent + f'{job} type <{job.type}>\n'
            for link in job:
                repr += indent*2 + f'{link}\n'
        return repr

    def __len__(self):
        return len(self.jobs)

    def __getitem__(self,index):
        return self.jobs[index]

    def __enter__(self):
        """ Wrapper to have similar behaviour to '_io.TextIOWrapper' """
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        """ Wrapper to have similar behaviour to '_io.TextIOWrapper' """
        return self._file.__exit__(exc_type, exc_value, traceback)

    def _set_parsers(self,link_ids:list[int],interblock=-1):
        f"""
        Generates an internal copy of the available linkjob parsers which 
        has the specifications of whether to parse specific Links as empty or
        not. 

        Parameters
        ----------
        link_ids : list[int]
            List of integrers that represent which types of Links to parse.
            If None or an empty list are provided it will attempt to parse
            all LinkJobs. If {EMPTYLINKFLAG} is in the list, all Links will be parsed 
            as_empty (This is usefull when only the structure of the file wants
            to be ascertained). Otherwise, only the specified Links will be parsed.
        interblock : int, optional
            Number used to indicate the text between LinkJob blocks or 
            InternalJobs, by default -1
        """
        parsers = available_linkjobs.copy()
        assert interblock not in parsers

        if not link_ids: # If no specific Links are provided parse everything that is available
            parsers[interblock] = GeneralLinkJob.as_empty

        elif EMPTYLINKFLAG in link_ids: # Signal to 
            for key in parsers:
                parsers[key] = parsers[key].as_empty
        else:
            for key in parsers:
                if key not in link_ids:
                    parsers[key] = parsers[key].as_empty
        
        parsers[interblock] = GeneralLinkJob.as_empty
        self._parsers = parsers

    def print_file_structure(self):
        """Display the structure of links and internal jobs of the file."""
        indent = "  "
        output = f"{self.__repr__():}\n"
        for intjob in self:
            output += indent + f"{intjob.__repr__():}\n"
            for link in intjob:
                output += indent*2 + f"{link.__repr__():}\n"
        print(output)

    def read(self):
        """
        Alias of update for consistency with GaussianInFile class. It 
        automatically removes empty LinkJobs
        """
        self.update()

    def close(self):
        """Alias to file.close for consistency with the io.TextIOBase class"""
        self._file.close()
    
    def update(self,clean=True,final_print=False):
        """
        Tries to fetch new data. If it exists it parses it appropiately
        otherwise it fails silently.

        Parameters
        ----------
        clean : Bool
            If True removes all the EmptyLinkJobs found (the default is True).
        final_print : Bool
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

        if self.jobs[0].number is None:
            self.jobs[0].guess_info()

        if clean:
            self.clean()

        if final_print:
            print(f"{self:!r} UPDATED")
    
    def clean(self):
        """Removes per each InternalJob stored all the EmptyLinkJobs."""
        for internaljob in self.jobs:
            internaljob.clean()

    def get_links(self,*link_ids):
        """
        Wrapper Method to get a list of Links with certain Ids across
        the different Internal Jobs.

        Parameters
        ----------
        *link_ids : int
            Integrers that correspond to the type of links to be return.

        Returns
        -------
        list
        """
        link_lists = [job.get_links(*link_ids) for job in self.jobs]
        return list(chain(*link_lists))
    
    # Generators and Coroutines for File Parsing
    def Reader(self,file):
        """ Generator for line by line reading without stopiteration """
        while True:
            old_pos = file.tell()
            line = file.readline()
            new_pos = file.tell()
            reached_EOF = old_pos == new_pos
            yield reached_EOF,line
    
    def BlockFetcher(self,EOF=-9999,interblock=-1):
        """
        Generator that yields the text sliced in blocks and their type.

        A block is an iterable of strings and its type refers to a token that
        can be recognized by the ._parsers variable, something in betweem
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
            block = []

            # Ask the Reader until a "start line" is found
            while not start:
                reached_EOF,line = next(Reader)
                if reached_EOF:
                    yield EOF, ''
                else:
                    start = re_enter.findall(line)
                    if not start:
                        block.append(line)
                    else: # Store the number of the Link
                        number = int(start[0])
            
            # When found yield it as a "InterBlock" and prepare Block
            if block:
                yield interblock, ''.join(block)
                block = [line,]
            else:
                block.append(line)
            
            # Now that the start of the Link has been found, accumulate lines
            ## until the end or termination line is found
            end = False
            while not end:
                reached_EOF,line = next(Reader)
                if reached_EOF:
                    target = block[-10:] + [line,]
                    terminated = re_termination.findall(''.join(target))
                    if terminated:
                        block.append(line)
                        break
                    else:
                        yield EOF, ''
                else:
                    end = bool(re_exit.findall(line))
                    end = end or bool(re_termination.findall(line))
                block.append(line)
            # when end found, do return type token and yield the block
            yield number, ''.join(block)
    
    def BlockHandler(self):
        """ Coroutine. Receives a block, chooses the parser and parses it """
        # Initialization
        parsers = self._parsers
        current_job = self.jobs[-1]
        first_link = True
        block_type, block = yield 'Initialization done'
        while True:
            parser = parsers[block_type]
            link =  parser(block)
            if link.number != 1:
                current_job.append(link)
            elif first_link:
                current_job.append(link)
                current_job.guess_info()
                first_link = False
            else:
                is_explicit = link.info.new_InternalJob
                if not is_explicit:
                    info = link.InternalJobInfo(self.jobs[-1].number+1,'linked',True)
                    link.info = info
                new_job = InternalJob()
                self.jobs.append(new_job)
                current_job = self.jobs[-1]
                current_job.append(link)
                current_job.guess_info()

            block_type, block = yield

class InternalJob(object):
    """
    Gaussian 09/16 InternalJob parent class, if any special type of Job
    requires different parsing it should be a subclass of this one.

    Parameters
    ----------
    number : int
        ordinal number of the InternalJob (the default is None).

    Attributes
    ----------
    type
        string identifier for the job.
    links
        List of the different Links that belong to the InternalJob.
    number

    """

    def __init__(self,number=None):
        self.number = number
        self.type = None
        self.links = []
    def __repr__(self):
        cls = type(self).__name__
        if self.number is None:
            return f'<{cls} Created but Empty>'
        else:
            return f'<{cls} {self.number}>'
    def __str__(self):
        return f'Internal Job {self.number}: {self.type}'
    def __getitem__(self,index):
        return self.links[index]
    def __len__(self):
        return len(self.links)

    def append(self,link):
        # Restrict to Link objects
        if not isinstance(link, LinkJob):
            raise TypeError(f'{link:!r} is not of class {LinkJob:!r}')
        self.links.append(link)
    def guess_info(self):
        """
        Guesses the number and type attributes of itself using the stored
        Links.
        """
        if self.links:
            links = (link for link in self.links if link.number == 1)
            try:
                start_link = next(links)
                info = start_link.info
            except AttributeError:
                pass
            except StopIteration:
                pass
            else:
                self.number = info.number
                self.type = info.type
    def clean(self):
        """Removes all the Empty Link instances within Links."""
        idx_to_remove = []
        for i, link in enumerate(self.links):
            if not link.text:
                idx_to_remove.append(i)
        for index in reversed(idx_to_remove):
            _ = self.links.pop(index)
    def get_links(self,*link_ids):
        """Wrapper Method to get a list of Links with certain Ids.

        Parameters
        ----------
        *link_ids : int
            Integrers that correspond to the type of links to be return.

        Returns
        -------
        list
            List of Link Objects ordered by appearance in the file and filtered
            by Link Number.
        """
        return [link for link in self.links if link.number in link_ids]

class GaussianInFile(object):
    """
    Gaussian 09/16 input file parent class, if any special type of input
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
    memory : int
        property to easily access and change the preprocessing['mem'] value
    extra_printout : bool
        Controls the behaviour of the command line's "#" or "#p". If True 
        the command line will appear with the "#p". It is True by default.
    solvent : str or None
        When set to None the solvation is still active but no solvent has been
        specified. When set to 'gas' the solvation is completely removed and 
        when set to a different string, the string is assumed to be a valid name
        for the solvent in gaussian.
    solvent_model : str or None
        This property controls the keyword scrf and has three possible values: 
        None -> the scrf keyword is not included
        'pcm' -> the 'pcm' suboption within scrf is included
        'smd' -> the 'smd' suboption within the scrf is included.  
    """
    def __init__(self,file=None):

        if isinstance(file,io.TextIOBase):
            self._file = file
        elif file is None: 
            self._file = None
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
        self.extra_printout = True
        self._nprocs_key = None
        self._mem_key = None

    def __repr__(self):
        cls = type(self).__name__
        name = 'unnamed'
        if hasattr(self._file,'name'): 
            name = self._file.name.split("/")[-1]
        #size = len(self)
        return f'<{cls}({name})>'

    def __str__(self):

        title = self.title
        if not title: 
            title = 'title_placeholder'
        
        structure = [
            self.preprocessing_as_str(),
            self.commandline_as_str(),
            '',
            title,
            '',
            ]

        # The geom=check, geom=read, geom=allcheck may require or not 
        # to include the geometries/charge&spin
        chrg_spin = self.charge_and_spin_as_str()
        if chrg_spin:
            structure.append(chrg_spin)

        geom = self.geometry_as_str()
        if geom: 
            structure.append(geom)
        
        structure.extend(['',
                          '\n\n'.join(self.tail)])
        str_repr = '\n'.join(structure).strip() +'\n\n\n'
        return str_repr

    def __len__(self):
        return len(str(self))

    def __enter__(self):
        """ Wrapper to have similar behaviour to '_io.TextIOWrapper' """
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        """ Wrapper to have similar behaviour to '_io.TextIOWrapper' """
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
        if self._nprocs_key is None: 
            self._search_nprocs_key()
        return self.preprocessing.get(self._nprocs_key,None)
    @nprocs.setter
    def nprocs(self,other):
        if self._nprocs_key is None: 
            self._search_nprocs_key()
        self.preprocessing[self._nprocs_key] = other

    @property
    def memory(self):
        if self._mem_key is None: 
            self._search_mem_key()
        return self.preprocessing.get(self._mem_key,None)
    @memory.setter
    def memory(self,other):
        if self._mem_key is None: 
            self._search_mem_key()
        self.preprocessing[self._mem_key] = other
    
    @property
    def solvent(self):
        if 'scrf' not in self.commandline: 
            return 'gas'
        items = self.commandline.get('scrf',[])
        for item in items: 
            if 'solvent=' in item:
                return item.split('=')[1]
        return None
    
    @solvent.setter
    def solvent(self,other):
        if other is None: 
            self._remove_solvent(remove_solvation=False)
        elif other.lower() == 'gas':
            self._remove_solvent(remove_solvation=True)
        else:
            self._set_solvent(other)
    
    @property
    def solvent_model(self):
        items = self.commandline.get('scrf',[])
        for item in items: 
            if item.lower() in ['pcm','smd']:
                return item
        return None
    
    @solvent_model.setter
    def solvent_model(self,other):
        if other is None: 
            self._remove_solvent_model()
        elif other.lower() in ['pcm','smd']: 
            self._set_solvent_model(other.lower())
        else:
            raise ValueError(f'solvent model "{other}" is not smd/pcm/None')

    def read(self):
        """
        Reads the file and populates the appropiate attributes.
        """
        txt = [line.rstrip() for line in self._file]
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

        # Define the ids of the block where the title starts
        start_idx = 0
        stop_idx = 1

        # Parse preprocessing if there is preprocessing
        preprocessing = []
        for line in header:
            if line.startswith("%"):
                preprocessing.append(line.lstrip("%"))
            elif line.startswith("#"):
                break
        else:
            # if the preprocessing and the command line are separated by newline
            # update accordingly
            header = iter(txt[bins[0]+1:bins[1]])
            line = next(header)
            start_idx = 1
            stop_idx = 2
        
        self.parse_preprocessing(preprocessing)

        commandlines = [line.strip(),]
        for line in header:
            commandlines.append(line.strip())
        self.parse_commandlines(commandlines)

        # Read the Title Section
        start = bins[start_idx]+1
        stop = bins[stop_idx]
        title = [line for line in txt[start:stop]]
        self.title = '\n'.join(title)

        # Read charge and spin
        charge,spin = txt[stop+1].split()
        self.charge,self.spin  = int(charge), int(spin)

        # Now we read the geometry
        start = stop+1
        stop = bins[stop_idx+1]
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
            self._file.write(self._txt)
        else:
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
    def parse_commandlines(self,lines:list[str]):
        """
        Parses the lines that contain the calculation commands keywords and
        transforms them into a dictionary representation.

        Parameters
        ----------
        lines : list[str]
            list of strings. Empty lines will be ignored.
        """
        lines = [line.split() for line in lines if line.strip()]
        # the first line contains the "#p"
        self.extra_printout = lines[0][0] == '#p'
        # parse all the keywords
        start = lines[0][1:]
        others = [i for i in chain(*lines[1:])]
        method_found = False
        basis_found = False
        for item in chain(start,others):
            dummy = item.split('/')
            if '/' in item and is_method(dummy[0]) and is_basis(dummy[1]):
                method_found = True
                basis_found = True
                method,basis = dummy
                self.method = method
                self.basis = basis
                continue
            elif is_method(item) and not method_found:
                method_found = True
                self._method = item
                key,val = item, []
            elif is_basis(item) and not basis_found:
                basis_found = True
                self._basis = item
                key,val = item, []
            elif is_basis(item) or is_method(item):
                warnings.warn('2 Basis or methods found \n')
                warnings.warn(f'taking {item} as a normal keyword')
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
                Aux3 = set(self.commandline[key]).union(set(val))
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

    # Helper functions for writing
    def preprocessing_as_str(self):
        """
        Transforms the preprocessing attribute to a suitable string 
        representation.

        Returns
        -------
        str
            string corresponding to the preprocessing part of the input file. 
        """
        preprocessing = []
        for key,val in self.preprocessing.items():
            if val:
                Aux = f'%{key}={val}'
            else:
                Aux = f'%{key}'
            preprocessing.append(Aux)
        return '\n'.join(preprocessing)
    def commandline_as_str(self): 
        """
        Transforms the commandline attribute to a suitable string 
        representation.

        Returns
        -------
        str
            string corresponding to the commandline of the input file. 
        """
        if self.extra_printout: 
            commandline = ['#p',]
        else:
            commandline = ['#',]
        for key,val in self.commandline.items():
            if val and (len(val) == 1):
                Aux = f"{key}={','.join(val)}"
            elif val:
                Aux = f"{key}=({','.join(val)})"
            else:
                Aux = f"{key}"
            commandline.append(Aux)
        return ' '.join(commandline)
    def charge_and_spin_as_str(self):
        """
        generates a string in the appropriate format for the charge and spin 
        information

        Returns
        -------
        str
            string with the charge and spin data
        """
        keywords = [k.lower() for k in self.commandline.get('geom',[])]
        if 'geom' in self.commandline and 'allcheck' in keywords:
            # the allcheck keyword does not require geometry nor charge and 
            # spin specification
            return ''
        return f'{self.charge} {self.spin}'
    def geometry_as_str(self): 
        """
        generates a string in the appropriate format for the geometry, based 
        on the available geometry attribute and the command line keywords

        Returns
        -------
        str
            string of the geometry
        """
        keywords = [k.lower() for k in self.commandline.get('geom',[])]
        is_geometry_fromfile = any([item in keywords for item in ['check','read','allcheck']])
        if 'geom' in self.commandline and is_geometry_fromfile:
            return ''
        return str(self.geometry).rstrip()

    # Private functions for properties management
    def _search_nprocs_key(self):
        for key in NPROCSHARED_ALIASES: 
            if key in self.preprocessing: 
                self._nprocs_key = key
                break
        else:
            self._nprocs_key = NPROCSHARED_ALIASES[0]
    def _search_mem_key(self):
        for key in MEMORY_ALIASES: 
            if key in self.preprocessing: 
                self._mem_key = key
                break
        else:
            self._mem_key = MEMORY_ALIASES[0]
    def _remove_nprocs(self): 
        key = self._nprocs_key
        self._nprocs_key = None
        return self.preprocessing.pop(key,None)
    def _remove_mem(self): 
        key = self._mem_key
        self._mem_key = None
        return self.preprocessing.pop(key,None)
    def _remove_solvent(self,remove_solvation=False):
        if remove_solvation: 
            _ = self.commandline.pop('scrf')
            return

        items = self.commandline.get('scrf',[])
        for i,item in enumerate(items): 
            if 'solvent=' in item:
                break
        else:
            return
        _ = items.pop(i)
        self.commandline['scrf'] = items
    def _set_solvent(self,other): 
        items = self.commandline.get('scrf',[])
        for i,item in enumerate(items): 
            if 'solvent=' in item:
                items[i] = f'solvent={other}'
                break
        else:
            items.append(f'solvent={other}')
        self.commandline['scrf'] = items
    def _remove_solvent_model(self):
        # If scrf does not exist
        current = self.solvent_model 
        if current is None: 
            return
        
        # Otherwise assume existence
        items = self.commandline.get('scrf',[])
        for i,item in enumerate(items):
            if item == current: 
                break
        _ = items.pop(i)
        self.commandline['scrf'] = items
    def _set_solvent_model(self,other):
        current = self.solvent_model 
        items = self.commandline.get('scrf',[])
        if current is None: 
            items.append(other)
        else:
            for i,item in enumerate(items):
                if item == current: 
                    items[i] = other
        
        self.commandline['scrf'] = items
        
    # Attribute modifying functions
    def geometry_from_chk(self,name=None,with_cs=False): 
        """
        Adds the specific keywords to read the geometry from the chk. 
        If no chk is present in the preprocessing options it will be 
        added. 

        Parameters
        ----------
        name : str, optional
            stem of the chk file. If none provided it will use the default 
            behavior of the add_chk method. 
        with_cs : bool, optional
            If enabled the allchk will be used instead of chk or read, which 
            leads to gaussian reading the charge and spin directly from the 
            chk. 
        """
        
        if 'chk' not in self.preprocessing: 
            self.add_chk(name)
        
        if with_cs: 
            self.add_kwd('allcheck',where='geom')
        else:
            self.add_kwd('check',where='geom')
        self.add_kwd('read',where='guess')
    def geometry_not_from_chk(self):
        """
        Undoes the effects of the geometry_from_chk method, except the creation
        of the chk. 
        """
        keywords = self.commandline.get('geom',[])
        for target in ['allcheck','check','read']: 
            keywords.pop(keywords.index(target))
        
        if keywords: 
            self.commandline['geom'] = keywords
        else: 
            self.pop_kwd('geom')

        keywords = self.commandline.get('guess',[])
        for target in ['allcheck','check','read']: 
            keywords.pop(keywords.index(target))
        
        if keywords: 
            self.commandline['guess'] = keywords
        else: 
            self.pop_kwd('guess')

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
    def disable_extra_printout(self):
        """
        When used, the string representation of the object will not include the 
        #p in the command line.
        """
        self.extra_printout = False
    def enable_extra_printout(self):
        """
        When used, the string representation of the object will include the 
        #p in the command line.
        """
        self.extra_printout = True
    def pop_kwd(self,keyword,where=None):
        """
        Removes a keyword from the command line and returns it. 

        Parameters
        ----------
        keyword : str
            keyword to remove from the command line.
        where : str, optional
            if provided it searches the keyword as a suboption of the "where" 
            keyword i.e. pop_kwd('smd',where='scrf') or 
            pop_kwd('cartesian',where='opt')

        Returns
        -------
        str or None
            Returns the removed keyword (or None when the keyword was not in the
            command line)  
        """
        if where is None:
            return self.commandline.pop(keyword,None)
        
        items = self.commandline.get(where,[])
        if keyword in items:
            kwd = items.pop(items.index(keyword)) 
            self.commandline[where] = items
            return kwd
        return None
    def add_kwd(self,keyword,where=None):
        """
        Adds a keyword to the command line. 

        Parameters
        ----------
        keyword : str
            keyword to add to the command line.
        where : str, optional
            if provided it adds the keyword as a suboption of the "where" 
            keyword i.e. add_kwd('smd',where='scrf') or 
            add_kwd('cartesian',where='opt')

        """
        if where is None:
            self.commandline[keyword] = []
            return
        
        items = self.commandline.get(where,[])
        if keyword in items:
            return
        items.append(keyword)
        self.commandline[where] = items
    def pop_kwds(self,*keywords,where=None): 
        """
        Removes the provided keywords from the command line and returns them. 

        Parameters
        ----------
        *keywords : str
            keywords to remove from the command line.
        where : str, optional
            if provided it searches the keywords as a suboption of the "where" 
            keyword i.e. pop_kwds('smd',where='scrf') or 
            pop_kwds('cartesian',where='opt')

        Returns
        -------
        list[str]
            Returns the removed keywords if they were in the command line
        """
        if where is None:
            output = [self.commandline.pop(k,None) for k in keywords]
            return [o for o in output if o is not None]
        
        items = self.commandline.get(where,[])
        output = []
        indices = [items.index(keyword) for keyword in keywords if keyword in items]
        output = [items.pop(idx) for idx in reversed(indices)]
        self.commandline[where] = items 
        return list(reversed(output))
    def add_kwds(self,*keywords,where=None): 
        """
        Adds keywords to the command line. 

        Parameters
        ----------
        keyword : str
            keyword to add to the command line.
        where : str, optional
            if provided it adds the keywords as suboptions of the "where" 
            keyword i.e. add_kwds('smd','solvent=water',where='scrf') or 
            add_kwds('cartesian','ts','noeigentest',where='opt')

        """
        if where is None:
            for keyword in keywords: 
                self.add_kwd(keyword)
            return
        
        items = self.commandline.get(where,[])
        keywords_to_add = [k for k in keywords if k not in items]
        items.extend(keywords_to_add)
        self.commandline[where] = items

    def pop_l0_kwd(self,keyword,where=None):
        """
        Removes a keyword from the preprocessing and returns it.

        Parameters
        ----------
        keyword : str
            keyword to remove from the preprocessing.
        where : str, optional
            if provided it searches the keyword as a suboption of the "where" 
            keyword i.e. pop_l0_kwd('job1.chk',where='oldchk')

        Returns
        -------
        str or None
            Returns the removed keyword (or None when the keyword was not in the
            preprocessing)  
        """
        if key in NPROCSHARED_ALIASES or where in NPROCSHARED_ALIASES:
            return self._remove_nprocs()

        if key in MEMORY_ALIASES or where in MEMORY_ALIASES:
            return self._remove_mem()

        if where is None:
            return self.preprocessing.pop(keyword,None)

        kwd = self.preprocessing.get(where,None)
        self.preprocessing[where] = ''
        return kwd
    def add_l0_kwd(self,keyword,where=None):
        """
        Adds a keyword to the preprocessing. 

        Parameters
        ----------
        keyword : str
            keyword to add to the preprocessing.
        where : str, optional
            if provided it adds the keyword as a suboption of the "where" 
            keyword i.e. add_l0_kwd('job1.chk',where='chk') or 
            add_kwd('40GB',where='mem')

        """
        if where in NPROCSHARED_ALIASES: 
            self.nprocs = keyword
        elif keyword in NPROCSHARED_ALIASES:
            self.preprocessing[keyword] = ''
        elif where in MEMORY_ALIASES: 
            self.memory = keyword
        elif keyword in MEMORY_ALIASES:
            self.preprocessing[keyword] = ''
        elif where is None:
            self.preprocessing[keyword] = ''
        else:
            self.preprocessing[where] = keyword

    @classmethod
    def from_str(cls,text):
        """
        Creates a GaussianInFile object from a gaussian input file read as a 
        string, parses it and populates the attributes of the class.  

        Parameters
        ----------
        text : str
            Contents of a Gaussian input file as a string. i.e. 
            .. code:: python
               
               with open('somefile.com','r') as F: 
                   text = F.read()
        name : str
            Name to represent the GaussianInFile 
        Returns
        -------
        GaussianInFile
            A populated GaussianInFile object
        """
        new = cls()
        new._file = text.split('\n')
        new.read()
        return new

class MultiGaussianInFile(object): 
    """
    Container class of multiple GausianInFiles linked together in a single 
    gaussian input calculation. 

    Parameters
    ----------
    file : io.TextIOBase or str (the default is None)
        File instance (Result of open(filename,'r')) or valid filename.

    Attributes
    ----------
    jobs : list
        List providing access to each one of the individual GaussianInFile
        objects representing each one of the linked calculations. 
    """
    def __init__(self,file=None):
        # Do Something
        if isinstance(file,io.TextIOBase):
            self._file = file
        elif file is None: 
            self._file = None
        else:
            self._file = open(file,'a+')
            if self._file.tell() != 0:
                self._file.seek(0)
        self.jobs = []

    def __repr__(self):
        cls = type(self).__name__
        file = self._file.name.split("/")[-1]
        #size = len(self)
        return f'<{cls}({file} with njobs={len(self.jobs)})>'
    def __str__(self):
        str_repr = '\n\n--Link1--\n'.join([str(job).rstrip() for job in self.jobs])
        # It is better to enfoce the condition here than for every 
        # possible case of reading a file or adding the tail
        if not str_repr.endswith('\n\n\n'):
            n = 3 - str_repr[-3:].count('\n') 
            str_repr += '\n'*n
        return str_repr
    def __len__(self):
        return len(str(self))

    def __enter__(self):
        ''' Wrapper to have similar behaviour to "_io.TextIOWrapper" '''
        return self
    def __exit__(self, exc_type, exc_value, traceback):
        ''' Wrapper to have similar behaviour to "_io.TextIOWrapper" '''
        return self._file.__exit__(exc_type, exc_value, traceback)

    def read(self):
        """
        Reads the file and populates the appropiate attributes.
        """
        if self._file is None: 
            return
        txt = self._file.read()
        self.jobs = [GaussianInFile.from_str(job.lstrip()) for job in txt.split('--Link1--')] 
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

    def enforce_same_chk(self,chk=None): 
        """
        Enforces the same chk file in all jobs. 

        Parameters
        ----------
        chk : str, optional
            filename of the chk to use, if none is provided it defaults to the 
            one of the first job. 
        """
        if chk is None: 
            chk = self.jobs[0].preprocessing['chk']
        for job in self.jobs: 
            job.add_l0_kwd(chk,where='chk')
    def enforce_continuous_chk(self,basename=None):
        """
        Enforces that the chk of job i-1 is retained and a copy of it is used 
        at the start of job i. For an example with two jobs, the first link0 
        section  will look like:

        ..
        
           %chk=basename_job0.chk
        
        Then the link0 of the second job will look like: 
        
        .. 
           %oldchk=basename_job1.chk
           %chk=basename_job1.chk

        Parameters
        ----------
        basename : str, optional
            basename of the chk to use, if none is provided it defaults to the 
            one of the first job. 
        """
        if basename is None: 
            basename = self.jobs[0].preprocessing['chk'].rsplit('.',maxsplit=1)[0]
        for (i,job) in enumerate(self.jobs):
            job.add_l0_kwd(f'{basename}_job{i}.chk',where='chk')
            if i-1 >= 0: 
                job.add_l0_kwd(f'{basename}_job{i-1}.chk',where='oldchk')
    def enforce_same_nprocs(self,nprocs=None):
        """
        Enforces the same nprocs in all jobs. 

        Parameters
        ----------
        nprocs : int, optional
            number of processors to use in a calculation, if none is provided it 
            defaults to the value of the nprocs of the first job. 
        """
        if nprocs is None: 
            nprocs = self.jobs[0].nprocs
        for job in self.jobs: 
            job.nprocs = nprocs
    def enforce_same_mem(self,mem=None):
        """
        Enforces the same memory in all jobs. 

        Parameters
        ----------
        mem : int, optional
            Memory to use in a calculation, if none is provided it 
            defaults to the value of the mem of the first job. 
        """
        if mem is None: 
            mem = self.jobs[0].mem
        for job in self.jobs: 
            job.mem = mem
    def enforce_same_method(self,method):
        """
        Enforces the same method in all jobs. 

        Parameters
        ----------
        method : str, optional
            method to use in a calculation, if none is provided it 
            defaults to the value of the first job's method. 
        """
        if method is None: 
            method = self.jobs[0].method
        for job in self.jobs: 
            job.method = method

# TODO: Implement a class to read and manipulate the basis functions in the tail
# class BasisTail(object), whose str function returns things as it should and
# that can have a linked input file object, so that modifying the basis of this
# object will modify the input file basis in certain cases.
