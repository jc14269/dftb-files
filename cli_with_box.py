# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.
"""
Demonstrates interactive molecular dynamics running with DFTB+ as the engine and ASE as the integrator.

Run with:

.. code bash
    python cli_with_box.py structure.xyz

"""
import argparse
import textwrap
#added in line below for periodic conditions
import numpy as np
#added in line below for periodic coniditons
from ase import Atoms
from ase import units, io
from ase.md import MDLogger
from ase.md.nvtberendsen import NVTBerendsen
# from ase.md.langevin import Langevin
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from narupa.ase.imd_server import ASEImdServer

from narupa.DFTB.DFTB_calculator import DFTBCalculator
from narupa.DFTB.Input_Generator import DFTBInputGenerator


def handle_user_arguments() -> argparse.Namespace:
    """
    Parse the arguments from the command line.

    :return: The namespace of arguments read from the command line.
    """
    description = textwrap.dedent("""\
    Run an ASE IMD simulation from an ASE supported file using Sparrow.
    """)
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument(
        '-structure_path',
        default='50_graphene.xyz',
        help='The structure to run in a format supported by ASE.',
    )
    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Display state information.',
    )

    # parser argument added by Joe Crossley-Lewis, starting on 29/10/2019 to allow periodicity
    


    parser.add_argument('-t', '--trajectory_port', default=None)
    parser.add_argument('-i', '--imd_port', default=None)
    parser.add_argument('-a', '--address', default=None)
    parser.add_argument('-f', '--frame_interval', default=1, type=int)
    parser.add_argument('-s', '--time_step', default=5.0)
    parser.add_argument('-m', '--method', default='DFTB0')
    arguments = parser.parse_args()
    return arguments


def start_imd(input_file, address=None, traj_port=None, imd_port=None, method='DFTB0', frame_interval=5, time_step=2.0,
              verbose=False):

    #DFTBin = DFTBInputGenerator()
    #DFTBin.xyz = input_file
    #DFTBin.create_input()

    atoms = io.formats.read(input_file)
    
    #block below is to make periodic conditions
    #import numpy as np
    #from ase import Atoms
    #a = atoms
    #a.cell[:]
    #Cell([1.0650000000E+01.,  -6.1487803650E+00.,   0.0000000000E+00.],
        # [1.0650000000E+01.,   6.1487803650E+00.,   0.0000000000E+00.],
        # [0.0000000000E+00.,   0.0000000000E+00.,   1.0000000000E+02.]]



    #Jonathans comments on how to potentially input a box
    #cell = ase.cell.Cell(...)
       # atoms.set_cell(cell)
    # print(atoms)
    # print(atoms.get_positions())

    atoms.set_calculator(DFTBCalculator(method=method))
    # print(atoms)
    print(f'Running dynamics')

    print(f'ASE energy: {atoms.get_potential_energy()}')

    MaxwellBoltzmannDistribution(atoms, 300 * units.kB)

    dyn = NVTBerendsen(atoms, 1 * units.fs, 300, (time_step * units.fs))


    if verbose:
        dyn.attach(MDLogger(dyn, atoms, '-', header=True, stress=False,
                            peratom=False), interval=100)

    imd = ASEImdServer(dyn,
                       address=address,
                       frame_interval=frame_interval,
                       )
    return imd


def main():
    """
    Entry point for the command line.
    """
    arguments = handle_user_arguments()
    #added argument for periodicity
    imd = start_imd(arguments.structure_path, arguments.address, arguments.trajectory_port, arguments.imd_port,
                    arguments.method, arguments.frame_interval, arguments.time_step, arguments.verbose)
    count = 1
    while True:
        imd.run(10)


if __name__ == '__main__':
    main()

