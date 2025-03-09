'''Conversion of .cif file to either FHIaims or DFTB+ geometry input files.'''

#AUTHOR: Ariadni Boziki

import os
from ase.io import read, write
from ase.io.aims import write_aims


class GeometryConversion():

    def __init__(self, code: str):

        self.code = code
        self.path = os.getcwd()
        self.atoms = None

    def _read_cif_structure(self, input_cif_name):
        """
        Retrieve the structure from a file of '.cif' format.
        """

        self.atoms = read(input_cif_name, format='cif')

    def aims_geom_input_to_extxyz(self, input_aims_name, extxyz_file_name):
        """
        Retrieve the structure from a file of geometry.in aims file.
        """

        self.atoms = read(input_aims_name, format='aims')
        write(extxyz_file_name, self.atoms, format='extxyz')

    def _write_structure(self):
        """
        Generate the geometry input file in either FHIaims or DFTB+ format.
        """
        
        if self.code == 'aims':
            write('geometry.in', self.atoms, format='aims')
        elif self.code == 'dftb+':
            write('geo.gen', self.atoms, format='dftb')
        elif self.code == 'so3lr':
            write('so3lr.xyz', self.atoms, format='extxyz')

    def check_cif_file(self):
        """
        Verifies the presence of a file with the .cif extension within the directory.
        """
    
        cif_file = [file for file in os.listdir(self.path) if file.endswith('.cif')]

        if cif_file:
            cif_file = cif_file[0]
            print(f'ONE .cif FILE EXISTS IN THE DIRECTORY: {cif_file}')
            self._read_cif_structure(cif_file)
            self._write_structure()
