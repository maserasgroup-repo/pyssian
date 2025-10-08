"""
Contains auxiliary classes, usefull in combination with the gaussianclasses.
"""
import os
import io
from warnings import warn
from math import ceil
from itertools import cycle
from pathlib import Path

import numpy as np

from .chemistryutils import PeriodicTable

class Geometry(object):
    """
    This class provides a basic interface and different construction methods 
    to simplify the extraction of the geometry of a molecule from different 
    sources 
    
    Currently: Link202, GaussianInFile, xyz file
    """
    # classattribute for object naming
    counter = 1
    def __init__(self):
        """ Initialize an empty geometry with generic title """
        self.coordinates = []
        self.atoms = []
        cls = type(self)
        self.title = str(cls.counter)
        cls.counter += 1
        self._coords_as_numpy = False

    def __repr__(self):
        return f'< Geometry {self.title} object >'

    def __str__(self):
        linef = " {:<}\t{: 0.09f}\t{: 0.09f}\t{: 0.09f}"
        if self._coords_as_numpy: 
            atoms = self.atoms.flatten().tolist()
            coordinates = self.coordinates.tolist()
        else:
            atoms = self.atoms
            coordinates = self.coordinates
        text = "\n".join([linef.format(str(atom),*xyz)
                                        for atom,xyz in zip(atoms,coordinates)])
        return text

    def write(self,fd:io.TextIOWrapper):
        """
        Write Geometry to a File

        Parameters
        ----------
        fd : str | bytes | os.PathLike
            file where the xyz formatted geometry will be written
        """
        fd.write(str(self))

    def update_atoms(self,mapping):
        """ Updates the atoms names from a mapping """
        if self._coords_as_numpy: 
            for i,atom in self.atoms.flatten().tolist(): 
                self.atoms[i] = mapping[atom]
        else:
            self.atoms = [mapping[atom] for atom in self.atoms]

    def to_numpy(self): 
        """
        Changes the coordinate representation and atoms into a numpy arrays.
        """
        self.atoms = np.array(self.atoms,dtype=object).reshape(-1,1)
        self.coordinates = np.array(self.coordinates)
        self._coords_as_numpy = True
    def to_python(self): 
        """
        Changes the coordinate representation and atoms into python lists.
        """
        self.atoms = self.atoms.flatten().tolist()
        self.coordinates = self.coordinates.tolist()
        self._coords_as_numpy = False

    @classmethod
    def from_L202(cls,L202):
        """
        Generate a Geometry instance from a Gaussian file. If it can't find a
        mapping from atomic Number to symbol it will store the string
        representation of the number
        """
        geometry = cls()
        orientation = L202.orientation
        for atom in orientation:
            geometry.coordinates.append(atom[3:])
            geometry.atoms.append(PeriodicTable[atom[1]])
        return geometry

    @classmethod
    def from_Input(cls,inputfile):
        """Generate a Geometry instance from a Gaussian Input File instance."""
        geometry = cls()
        for line in str(inputfile.geometry).split('\n'):
            atid, *xyz = line.strip().split()
            geometry.atoms.append(atid)
            geometry.coordinates.append(tuple(map(float,xyz)))
        return geometry

    @classmethod
    def from_xyz(cls,xyzfile):
        """Generate a Geometry instance from a xyz file."""
        geometry = cls()
        with open(xyzfile,'r') as F:
            _iter = F.__iter__()
            n = int(next(_iter).strip())
            _ = next(_iter)
            for _ in range(n):
                line = next(_iter)
                atid, *xyz = line.strip().split()
                geometry.atoms.append(atid)
                geometry.coordinates.append(tuple(map(float,xyz)))
        return geometry

    def to_xyz(self,title:str|None=None) -> str:
        if title is None: 
            title = f'{self.title}'
        return f'{len(self.atoms)}\n{title}\n{self}\n'

class Cube(object):
    """
    Representation of a Gaussian .cube File generated with the cubegen tool.
    This class facilitates some operations of the cubeman tool allowing
    operations of more than 2 .cube files.

    Files with more than 1 cube are accepted but correct behaviour is not 
    guaranteed.

    Attributes
    ----------
    natoms : int
        total number of atoms
    type : str
        Type of cube. Currently only {'MO','other'} are considered
    is_unrestricted : bool
        True if the line that contains the number of atoms has 5 elements
        instead of 4.
    origin : tuple
        origin of the cube.
    shape : tuple
        Number of points in each of the three directions.
    basis : list
        Vectors of the three main directions.
    atoms : list
        List of the AtNumber, Charge, (x,y,z) per each atom.
    nMO : int
        Number of molecular orbitals
        (Current implementation is designed for only one MO per Cube)
    MO_Ids : list
        Ids of the MO orbitals.
        (Current implementation is designed for only one MO per Cube)
    matrix : np.array
        Matrix of values of the Cube.
    """

    def __init__(self):
        self.natoms = 0
        self.type = 'empty'
        self.is_unrestricted = False
        self.origin = (0,0,0)
        self.shape = (0,0,0)
        self.basis = [(0,0,0),
                      (0,0,0),
                      (0,0,0)]
        self.atoms = []
        self.nMO = 0
        self.MO_Ids = []
        self.matrix = []
    def __repr__(self):
        return f'<{type(self).__name__}_{id(self)}>'

    def __add__(self,other):
        new_item = self.__class__()
        for attr in 'natoms type origin shape basis atoms nMO'.split():
            if getattr(self,attr) != getattr(self,attr):
                raise ValueError(f'{self}.{attr} != {other}.{attr}')
            else:
                setattr(new_item,attr,getattr(self,attr))
            new_item.matrix = self.matrix + other.matrix
        return new_item
    def __radd__(self,other):
        return self + other
    def __mul__(self,other):
        new_item = self.__class__()
        for attr in 'natoms type origin shape basis atoms nMO'.split():
            setattr(new_item,attr,getattr(self,attr))
        new_item.matrix = self.matrix*other
        return new_item
    def __rmul__(self,other):
        return self*other
    def __neg__(self):
        return self*(-1)
    def __sub__(self,other):
        return self + (-other)
    def __pow__(self,other):
        new_item = self.__class__()
        for attr in 'natoms type origin shape basis atoms nMO'.split():
            setattr(new_item,attr,getattr(self,attr))
        new_item.matrix = self.matrix**other
        return new_item

    # Writing functions
    def write(self,ofile):
        """
        Writes a cube to a File. 

        Parameters
        ----------
        ofile : str
            A valid path to a file. 
        """
        # Formats
        atom_format  = '{0: >5}   {1: .6f}   {2: .6f}   {3: .6f}   {4: .6f}\n'.format
        basis_format = '{0: >5}   {1: .6f}   {2: .6f}   {3: .6f}\n'.format

        # preparation
        natoms = str(self.natoms)
        sign = ' '
        if self.type == 'MO':
            sign = '-'
        x0,y0,z0 = self.origin
        (ix,jx,kx), (iy,jy,ky), (iz,jz,kz) =  self.basis
        Nx,Ny,Nz = self.shape
        MO_Ids = self.MO_Ids
        if self.nMO > 0 and not self.MO_Ids:
            MO_Ids = [1,] 
        n_MO_line = ['{: >5}'.format(str(self.nMO)),]
        n_MO_line.extend(list(map('{: >5}'.format, MO_Ids)))
        # Writing
        with open(ofile,'w') as F:
            F.write('Default Title Line\n')
            F.write('Default Other Title Line\n')
            line = f'{sign+natoms: >5}   {x0: .6f}   {y0: .6f}   {z0: .6f}\n'
            if self.is_unrestricted:
                line += '   1'
            F.write(line)
            F.write(basis_format(Nx,ix,jx,kx))
            F.write(basis_format(Ny,iy,jy,ky))
            F.write(basis_format(Nz,iz,jz,kz))
            for atom in self.atoms:
                F.write(atom_format(*atom))
            F.write(''.join(n_MO_line)+'\n')
            self.write_values(F)
    def write_values(self,fd:io.TextIOWrapper):
        """
        Handles the writing in of the values of self.matrix with the fortran77
        format

        Parameters
        ----------
        fd : io.TextIOWrapper
            the output of a open(Filename,'w')
        """
        # Formats
        val_format = ' {: .5E}'.format
        # preparation
        Nx,Ny,Nz = self.shape
        MO_Ids = self.MO_Ids
        counter = 0
        for i in range(Nx):
            for j in range(Ny):
                if len(MO_Ids) > 1:
                    for Id,_ in enumerate(MO_Ids):
                        for k in range(Nz):
                            counter += 1
                            fd.write(val_format(self.matrix[i][j][k][Id]))
                            if counter == 6:
                                fd.write('\n')
                                counter = 0
                    if counter != 0:
                        counter = 0
                        fd.write('\n')
                else:
                    for k in range(Nz):
                        counter += 1
                        fd.write(val_format(self.matrix[i][j][k]))
                        if counter == 6:
                            fd.write('\n')
                            counter = 0
                    if counter != 0:
                        counter = 0
                        fd.write('\n')

    # Constructor methos
    @classmethod
    def from_file(cls,file):
        """
        Creates a Cube instance from a cube file.

        Parameters
        ----------
        file : str
            A valid filename to an existing .cube file.

        Returns
        -------
        Cube
            Instantiated cube from the file.
        """
        newcube = cls()
        with open(file,'r') as F:
            lines = F.readlines()
        _iter = lines.__iter__()
        newcube.handle_title_lines(_iter)
        newcube.read_origin_line(_iter)
        newcube.read_basis_lines(_iter)
        newcube.read_atoms_lines(_iter)
        newcube.read_MO_lines(_iter)
        newcube.read_matrix_lines(_iter)

        return newcube

    # Parsing functions
    def handle_title_lines(self,iterable):
        # Discard first two lines
        _ = next(iterable)    # {Title1} MO={1} || {Title2} MO={2}
        _ = next(iterable)    # MO coeffic + MO coeffic
    def read_origin_line(self,iterable):
        # Read origin and number of atoms
        line = next(iterable) # -{Natoms} {xorigin} {yorigin} {zorigin}
        Aux = line.strip().split()
        natoms,x0,y0,z0 = Aux[:4]
        if len(Aux) > 4:
            self.is_unrestricted = '1' == Aux[5]
        x0,y0,z0 = map(float,(x0,y0,z0))
        if int(natoms) < 0:
            self.type = 'MO'
        else:
            self.type = 'other'
        natoms = abs(int(natoms))
        self.natoms = natoms
        self.origin = (x0,y0,z0)
    def read_basis_lines(self,iterable):
        # Read number of points and directions of the basis
        line = next(iterable) # {Nx} {ix} {jx} {kx}
        Nx, ix,jx,kx = line.strip().split()
        ix,jx,kx = map(float,(ix,jx,kx))
        line = next(iterable) # {Ny} {iy} {jy} {ky}
        Ny, iy,jy,ky = line.strip().split()
        iy,jy,ky = map(float,(iy,jy,ky))
        line = next(iterable) # {Nz} {iz} {jz} {kz}
        Nz, iz,jz,kz = line.strip().split()
        iz,jz,kz = map(float,(iz,jz,kz))
        Nx,Ny,Nz = map(int,(Nx,Ny,Nz))
        self.shape = (Nx,Ny,Nz)
        self.basis = [(ix,jx,kx),(iy,jy,ky),(iz,jz,kz)]
    def read_atoms_lines(self,iterable):
        # Read the atoms, charges, and coordinates
        atoms = []
        for _ in range(self.natoms):
            line = next(iterable)
            AtN,Chrg,x,y,z = line.strip().split()
            AtN = int(AtN)
            Chrg,x,y,z = map(float,(Chrg,x,y,z))
            atoms.append((AtN,Chrg,x,y,z))
        self.atoms = atoms
    def read_MO_lines(self,iterable):
        ## For now I'm gonna treat it as a single line
        # Now I guess this is the "nMO" and "Ids..." line/s
        line = next(iterable)
        nMO,*MO_Ids = line.strip().split()
        self.MO_Ids = MO_Ids
        self.nMO = int(nMO)
    def read_matrix_lines(self,iterable):
        Nx,Ny,Nz = self.shape
        nMO = self.nMO
        nEntries = Nx*Ny
        nVals = Nz*nMO
        nlines = ceil(nVals/6)
        if nMO > 1:
            matrix = np.zeros(shape=(Nx,Ny,Nz,nMO))
        else:
            matrix = np.zeros(shape=(Nx,Ny,Nz))

        data = []
        values = []
        for index,line in zip(cycle(range(nlines)),iterable):
            values.extend(list(map(np.float64,line.strip().split())))
            if index == nlines-1:
                data.append(values)
                values = []
        else:
            if values:
                raise RuntimeError('Wrong developer assumptions on the .cube format')

        for i in range(Nx):
            for j in range(Ny):
                index1 = i*Ny + j
                items = data[index1]
                if nMO > 1:
                    for l in range(nMO):
                        for k in range(Nz):
                            index2 = l*Nz + k
                            matrix[i,j,k,l] = items[index2]
                else:
                    for k in range(Nz):
                        matrix[i,j,k] = items[k]
        self.matrix = matrix

