================
pyssian examples
================

.. contents:: Table of Contents
   :backlinks: none
   :local:


GaussianOutFile
...............

*pyssian.GaussianOutFile* can be opened as a file object. When reading an
already finished calculation file it is recommended to use the with statement as
in this first example. If the file is being written the recommended way to open
the file is as example 2 or 3.

.. code:: python

   from pyssian import GaussianOutFile
   # Example 1
   with GaussianOutFile(FilePath) as GOF1:
      GOF1.read()
   # Example 2
   GOF2 = GaussianOutFile(FilePath)
   GOF2.read()
   GOF2.close()
   # Example 3
   File = open(FilePath,'r')
   GOF3 = GaussianOutFile(File)
   GOF3.read()
   File.close()


It also accepts another parameter that indicates which Links are of interest.
The default behaviour (without specifying any value) is to consider all links.

A second option, only used to see the structure of the file without actually
storing any information of each link is to pass -1 as the first value:

.. code:: python

   File = open(FilePath,'r')
   with GaussianOutFile(File,parselist=[-1,]) as GOF:
       GOF.update(clean=False) # read is an alias for update(clean=True)


If you specify the numbers of the links to the GOF only those will be parsed and
will not be cleaned afterwards.

.. code:: python

   File = open(FilePath,'r')
   with GaussianOutFile(File,[1,103,202,502,9999]) as GOF:
       GOF.read()


GaussianInFile
..............

*pyssian.GaussianInFile* can be instantiated either from an existing input file
or to create a new file.

.. note::

   Currently it is in an early stage as proper support for method-basis
   management as well as oniom and zmatrix support require special attention.

The following Code snippet shows how to create a new input from an existing
one changing the geometry and method but retaining the rest of the options

.. code:: python

   from pyssian import GaussianInFile

   initial_theory_file = 'InitialTheory.in'
   initial_geometry_file = 'InitialGeometry.in'
   output_file = 'Output.in'
   with GaussianInFile(initial_theory_file) as theory_file:
       theory_file.read()
   with GaussianInFile(initial_geometry_file) as geometry_file:
       geometry_file.read()
   old_geometry = theory_file.geometry # In case we want to use it somewhere else
   theory_file.geometry = geometry_file.geometry
   theory_file.method = 'b3lyp'
   with open(output_file,'w') as F:
       theory_file.write(F)


It combines fairly well with pyssian.classutils.Geometry to create inputs from
outputs. The following code snippet is an example of how to create an input to
continue an optimization that failed due to exceeding the number of optimization
steps.

.. code:: python

   from pyssian import GaussianInFile, GaussianOutFile
   from pyssian.classutils import Geometry

   with GaussianOutFile('Old_Calc.out') as GOF:
      GOF.read()

   # Get the last geometry of the calculation
   geom = Geometry.from_l202(GOF.get_links(202)[-1])

   # Get the Link1 of the GaussianOutFile
   Link1 = GOF.get_links(1)[0]

   # Extract the calculation type and commands
   commandline = Link1.commandline
   nprocs = Link1.nprocs
   mem = Link1.mem
   Link0 = Link1.link0

   # Get Charge and spin from Link101
   Link101 = GOF.get_links(101)[0]
   charge = Link101.charge
   spin = Link101.spin

   # Now write the new input file
   with GaussianInFile('New_Calc.in') as GIF:
       GIF.parse_commandline([commandline,])
       # We can instead set a dict for the variable GIF.commandline
       # "GIF.commandline = {'opt':'','freq':'NoRamman','b3lyp':''}"
       # but using parse_commandline is easier in this case.
       GIF.preprocessing = {key:'' for key in Link0}
       GIF.preprocessing['nprocshared'] = nprocs
       GIF.preprocessing['mem'] = mem
       GIF.title = 'New Title'
       GIF.spin = spin
       GIF.charge = charge
       GIF.geometry = geom
       GIF.write()


Cube
....

*pyssian.classutils.Cube* class was introduce to simplify the sometimes a bit 
bothersome usage of cubeman from gaussian to add, substract, multiply... cube 
files. You can initialize an empty cube and populate it yourself but the class 
was thought to be used as follows: 

.. code:: python

   from pyssian.classutils import Cube 
   MO_1 = Cube.from_file('MO_01.cube')
   MO_2 = Cube.from_file('MO_02.cube')
   MO_3 = Cube.from_file('MO_03.cube')
   FinalCube = MO_1*2 + MO_2 - MO_3**2
   FinalCube.write('Final.cube') 


LinkJob
.......

The *pyssian.LinkJob* class is the base class for all the LinkJob subclasses
And implements the two basic attributes of all Links, *.number* and *.text*.
Currently the specific parsers implemented are:

- Link1
- Link101
- Link103
- Link123
- Link202
- Link502
- Link508
- Link601
- Link716
- Link804
- Link913
- Link914

.. code:: python

   with GaussianOutFile(File) as GOF:
       GOF.read()
   Link = GOF[0][RandomPosition]
   # General Attributes of all LinkJob classes
   print(Link.number)
   print(Link.text)


Link1
+++++

.. code:: python

   # From the file Get the first Link1
   Link1 = GOF.get_links(1)[0]
   # Attributes of Link1
   Link1.commandline
   Link1.nprocs
   Link1.mem
   Link1.link0
   Link1.IOps
   Link1.info # Will be deprecated in the future


Link101 & Link103
+++++++++++++++++

.. code:: python

   Link101 = GOF.get_links(101)[0]
   Link101.spin
   Link101.charge

   Link103 = GOF.get_links(103)[0]
   Link103.mode
   Link103.state
   Link103.conversion
   Link103.parameters
   Link103.stepnumber
   Link103.scanpoint
   if Link103.mode == 'Iteration':
       Link103.print_convergence()

Link123
+++++++

.. code:: python

   Link123 = GOF.get_links(123)[0]
   Link123.orientation
   Link123.step
   Link123.direction
   Link123.reactioncoord


Link202
+++++++

.. code:: python

   Link202 = GOF[-1].get_links(202)[0]
   Link202.orientation
   Link202.DistanceMatrix
   Link202.print_orientation()
   Link202.get_atom_mapping()

Link502 & Link508
+++++++++++++++++

.. code:: python

   list_of_links = GOF.get_links(502,508)
   energies = [link.energy for link in list_of_links if link.energy is not None]

Link601
+++++++

.. code:: python

   Link601 = GOF[-1].get_links(601)[-1]
   Link601.mulliken
   Link601.mulliken_heavy

Link716
+++++++

.. code:: python

   Link716 = GOF[-1].get_links(716)[-1]
   Link716.mode
   Link716.dipole
   Link716.units
   Link716.zeropoint
   Link716.thermal_energy
   Link716.enthalpy
   Link716.gibbs
   Link716.EContribs
   Link716.IRSpectrum

Link804 & Link913
+++++++++++++++++

.. code:: python

   Link804 = GOF.get_links(804)[-1]
   Link804.MP2
   Link804.SpinComponents
   scs_corr = Link804.get_SCScorr()
   HF_energy = GOF.get_links(502)[-1].energy 
   scs_energy = HF_energy + scs_corr

   Link913 = GOF.get_links(913)[-1]
   Link913.MP4
   Link913.CCSDT


Link914
+++++++

.. code:: python

   Link914 = GOF.get_links(914)[-1]
   for es in Link914.excitedstates: 
       number, energy, wavelen, OStrength, s2, transitions = es
       for transition in transitions: 
           donor = transition.donor
           acceptor = transition.acceptor 
           contribution = transition.contribution
           print(f'{donor} -> {acceptor}     {contribution}')
   # which can be done for the excited states 2,5,6: 
   Link914.print_excitedstates(2,5,6,show_transitions=True)


Usage Examples
..............

Code snippet to extract the last potential energy and geometry

.. code:: python

   from pyssian import GaussianOutFile

   MyFile = 'path-to-file'
   with GaussianOutFile(MyFile) as GOF:
      GOF.read()

   final_geometry = GOF.get_links(202)[-1].orientation
   last_potential_energy = GOF.get_links(502)[-1]
   print(last_potential_energy)
   print(str(final_geometry))


Code snippet to display 'Filename HF MP2 MP2(SCS)'

.. code:: python

   from pyssian import GaussianOutFile

   MyFile = 'path-to-file'
   with GaussianOutFile(MyFile,[1,502,804]) as GOF:
      GOF.read()

   HF = GOF.get_links(502)[-1].energy
   Link804 = GOF.get_links(804)[-1]
   MP2 = Link804.MP2
   MP2scs = HF + Link804.get_SCScorr()
   print(f'{MyFile}\t{HF}\t{MP2}\t{MP2scs}')


Code Snippet to follow a file being written by gaussian

.. code:: python

   from time import sleep

   from pyssian import GaussianOutFile

   GOF = GaussianOutFile(MyFile,[-1,])
   GOF.update(clean=False)
   print(GOF[-1][-1])
   sleep(10)
   GOF.update(clean=False)
   print(GOF[-1][-1])
   GOF.close()
