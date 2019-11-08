# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.

# This file contains all the DFTB+ APIs

from typing import Optional, Collection
import numpy as np
from ase import Atoms
from ctypes import CDLL, c_int, c_char_p, c_double, pointer
import ctypes
from numpy.ctypeslib import ndpointer
from ase.calculators.calculator import Calculator, all_changes

EV_PER_HARTREE = 27.2114
ANG_PER_BOHR = 0.529177


class DFTBCalculator(Calculator):

    implemented_properties = ['energy', 'forces']

    def __init__(self, atoms: Optional[Atoms] = None, method='DFTB0', **kwargs):
        super().__init__(**kwargs)
        self.atoms = atoms
        self.method = method
        try:
            self.DFTBpluslib = CDLL('./libdftb+.so')  # Load DFTB+ library
        except Expetions as error:
            print('Unable to load DFTB+ library', error)
        output_location = "tempfile.out"  # Output file
        output_mutable_string = ctypes.create_string_buffer(str.encode(output_location))
        input_location = "dftb_in.hsd"  # Input file
        input_mutable_string = ctypes.create_string_buffer(str.encode(input_location))

        try:
            self.DFTB_state = ctypes.c_double(0.0)  # Handler to output file
        except:
            print('Unable to set the DFTB+ handler for the output file')

        try:
            self.DFTB_handler = ctypes.c_double(0.0)  # Handler to input file
        except:
            print('Unable to set the DFTB+ handler for the input file')

        try:
            self.__wrap_function(self.DFTBpluslib, "dftbp_init", None, [ctypes.POINTER(ctypes.c_double),
                                                                        ctypes.c_char_p])
            self.DFTBpluslib.dftbp_init(ctypes.pointer(self.DFTB_state), output_mutable_string)  # Initialize
        except:
            print('Unable to initialize DFTB+')

        try:
            self.__wrap_function(self.DFTBpluslib, "dftbp_get_input_from_file", None,
                                 [ctypes.POINTER(ctypes.c_double), ctypes.c_char_p, ctypes.POINTER(ctypes.c_double)])
            self.DFTBpluslib.dftbp_get_input_from_file(ctypes.pointer(self.DFTB_state), input_mutable_string,
                                                       ctypes.pointer(self.DFTB_handler))  # Read input from file
        except:
            print('Unable to read the input file')

        try:
            self.__wrap_function(self.DFTBpluslib, "dftbp_process_input", None,
                                 [ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double)])
            self.DFTBpluslib.dftbp_process_input(ctypes.pointer(self.DFTB_state),
                                                 ctypes.pointer(self.DFTB_handler))  # Process input
        except:
            print('Unable to process the input file')

    def __wrap_function(self, lib, funcname, restype, argtypes):  # Wrapper function
        """Simplify wrapping ctypes functions"""
        func = lib.__getattr__(funcname)
        func.restype = restype
        func.argtypes = argtypes

    def calculate(self, atoms: Optional[Atoms] = None,
                  properties=('energy', 'forces'),
                  system_changes=all_changes):

        if atoms is None:
            atoms = self.atoms
            raise ValueError('No ASE atoms supplied to calculator, and no ASE atoms supplied with initialisation.')
        self._calculate_dftb(atoms, properties)

    def _calculate_dftb(self, atoms: Atoms, properties: Collection[str]):
        # self.__wrap_function(self.DFTBpluslib, "dftbp_get_nr_atoms", c_int, [ctypes.POINTER(ctypes.c_double)])
        # natoms = self.DFTBpluslib.dftbp_get_nr_atoms(self.DFTB_state)  # Get number of atoms
        natoms = len(atoms)
        # Generate general numpy matrices that are c compatible
        _doublepp = ndpointer(dtype=np.float64, ndim=1, flags='C')
        x = np.arange(natoms * 3)
        y = np.arange(1)
        positions_dummy = (x.__array_interface__['data'][0] + np.arange(x.shape[0]) * x.strides[0]).astype(np.float64)
        # print(positions_dummy)
        positions = atoms.positions
        positions = positions.reshape(1, natoms*3)
        positions = positions[0]
        # print(positions)
        positions = positions /ANG_PER_BOHR
        gradients_hartree_bohr = (x.__array_interface__['data'][0] + np.arange(x.shape[0])
                                  * x.strides[0]).astype(np.float64)
        energy_hartree = (y.__array_interface__['data'][0] + np.arange(y.shape[0]) * y.strides[0]).astype(np.float64)


        try:
            self.__wrap_function(self.DFTBpluslib, "dftbp_set_coords", None, [ctypes.POINTER(ctypes.c_double),
                                                                              _doublepp])
            self.DFTBpluslib.dftbp_set_coords(ctypes.pointer(self.DFTB_state), positions)  # Set coordinates
        except:
            print('Unable to set the coordinates')

        kwargs = {property_name: True for property_name in properties}

        if 'energy' in properties:
            try:
                self.__wrap_function(self.DFTBpluslib, "dftbp_get_energy", None, [ctypes.POINTER(ctypes.c_double),
                                                                                  _doublepp])
                self.DFTBpluslib.dftbp_get_energy(ctypes.pointer(self.DFTB_state), energy_hartree)  # Get energy
            except:
                print('Unable to extract the energy')
            self.results['energy'] = energy_hartree[0] * EV_PER_HARTREE

        if 'forces' in properties:
            try:
                self.__wrap_function(self.DFTBpluslib, "dftbp_get_gradients", None, [ctypes.POINTER(ctypes.c_double),
                                                                                     _doublepp])
                self.DFTBpluslib.dftbp_get_gradients(ctypes.pointer(self.DFTB_state),
                                                     gradients_hartree_bohr)   # Get forces
                gradients_hartree_bohr = gradients_hartree_bohr.reshape(natoms, 3)
            except:
                print('Unable to extract the forces')
            self.results['forces'] = - gradients_hartree_bohr * (EV_PER_HARTREE / ANG_PER_BOHR)
        return

