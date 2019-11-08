# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.

# This file launches and executes DFTB+

import numpy as np
from Input_Generator import *
from DFTB_calculator import *
import time
from ase import io


# Initialize DFTB+ input generator
DFTBin = DFTBInputGenerator()

# Set the DFTB+ parameters
DFTBin.scc = 'no'  # Choose among 'yes' or 'no'. Default value set is 'yes'
DFTBin.xyz = 'xyz_files/50_graphene.xyz'  # Default is 'Random.xyz'
DFTBin.diagonalization = 'DivideAndConquer'  # Choose among 'RelativelyRobust' or 'DivideAndConquer' or 'QR'. Default
#                                              is 'DivideAndConquer'
DFTBin.charge = 0.0  # Default is 0.0

# Generate DFTB+ input (dftb_in.hsd)
DFTBin.create_input()

# Initialize DFTB+ calculation instance
DFTB = DFTBCalculator()
atoms = io.formats.read(DFTBin.xyz)

# # Extract positions and forces from DFTB+
start_time = time.time()
for i in range(1, 1000):
    print("Step = " + str(i) + '\n')
    DFTB.calculate(atoms)  # Get positions and coordinates
    # print('\nPositions')
    # print(DFTB.xyz)  # Print coordinates
    # print('\nForces')
    # print(DFTB.forces)  # Print forces
    # print('\nEnergy')
    # print(DFTB.energy)  # Print energy
end_time = time.time()
print('Time = ' + str(end_time - start_time))

