# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.

# This file generates the input file (dftb_in.hsd) which is compatible with DFTB+


import numpy as np
import os


class DFTBInputGenerator:

    def __init__(self):  # Constructor
        self.xyz = 'xyz_files/50_graphene.xyz'  # Default input xyz file
        self.scc = 'yes'  # Default scc
        self.diagonalization = 'DivideAndConquer'  # Default diagonalization routine
        self.charge = 0.0  # Default charge

    def create_input(self):
        filename = self.xyz
        mod_filename = 'xyz_files/mod.xyz'
        with open(filename) as fp:
            lines = fp.readlines()
        lines[0] = lines[0].replace(' ', '')
        values = lines[0].split(' ')
        if len(values) == 1:
            del lines[0]
            del lines[0]
            for i in range(0, len(lines)):
                if lines[i] == ' \n':
                    del lines[i]
            fp = open(mod_filename, 'w')
            for i in range(0, len(lines)):
                fp.write(lines[i])
            fp.close()

        atoms = []
        geom = np.genfromtxt(mod_filename, dtype=None,
                             usecols=(1, 2, 3))  # Convert input to xyz coordinates
        with open(mod_filename) as fp:
            for ln in fp:
                values = ln.split(' ')
                atoms.append(values[0])
        natoms = len(atoms)
        elements = np.unique(atoms)  # Figure out the elements in the system
        numberedatoms = np.empty(natoms)  # List to keep the number of atoms
        mio = os.path.dirname(os.path.realpath(__file__)) + "/mio/"  # Path to the folder containing the .skf files

        # Generate the atom type list by labeling each element according to its appearance in the atoms array
        for idx, i in enumerate(elements):
            for jdx, j in enumerate(atoms):
                if i == j:
                    numberedatoms[jdx] = idx + 1

        # Now we check the position of the elements in the list and give the numbers according to how
        # They appear in the elements list. This is a requirement of the list
        geom = np.c_[numberedatoms, geom]

        # Add a column with increasing numbers to enumerate the atoms using arange+1
        geom = np.c_[np.arange(natoms) + 1, geom]

        setupg = "Geometry = GenFormat { \n" \
                 + str(natoms) + " C \n"  # This indicates that the system is a non periodic cluster
        setupg = setupg + str(" ".join(elements)) + " \n"  # Element types need to be put here
        setup = '''}
        
         Driver = ConjugateGradient {
           MovedAtoms = 1:-1               # Move all atoms in the system
           MaxForceComponent = 1.0e-4      # Stop if maximal force below 1.0e-4
           MaxSteps = 0                    # Stop after maximal 100 steps
           OutputPrefix = "geom.out"       # Final geometry in geom.out.{xyz,gen}
           ConvergentForcesOnly=no         # Stop the program from exiting badly when the scc doesnt converge
         }
         options ={
         WriteDetailedOut = No
         }   
         analysis ={
         WriteBandOut = No         
         }
         
         Hamiltonian = DFTB {
         # SCC = ''' + self.scc + '''
         # SCCTolerance = 1.0E-5
         # MaxSCCIterations = 5
         # Charge = ''' + str(self.charge) + '''
          MaxAngularMomentum = {\n'''

        for i in elements:
            if i == 'C':
                setup = setup + '''C = "p"\n'''
            elif i == 'N':
                setup = setup + '''N = "p"\n'''
            elif i == 'H':
                setup = setup + '''H = "s"\n'''
            elif i == 'O':
                setup = setup + '''O = "p"\n'''

        setup = setup + ''' 
          }
          solver = ''' + self.diagonalization + '''{}            # Specify the diagonalization routine that you want to use
          SlaterKosterFiles = Type2FileNames {    # File names from two atom type names
              Prefix = "mio/"  # Path as prefix
              Separator = "-"         # Dash between type names
              Suffix = ".skf"         # Suffix after second type name
              LowerCaseTypeName = No
          }
        }
        '''

        with open('dftb_in.hsd', 'w') as myfile:  # Write the dftb_in.hsd file
            myfile.write(setupg)
            np.savetxt(myfile, geom, fmt='%4i %2i %15.10f %15.10f %15.10f', delimiter=' ')
            myfile.write(setup)
