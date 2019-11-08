#!/usr/bin/env python3

from ase import io
#change the jobname to the name of the file you want to change
atoms = io.read('MFI-st1-MePro-OC-opt.cif')
atoms.write('MFI-st1-MePro-OC-opt.gen', format='gen')
