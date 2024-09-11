'''Prepares perturbed structures to facilitate the application of the finite displacements method, 
enabling the subsequent calculation of molecular vibrational frequencies.'''

#AUTHOR: Ariadni Boziki

import numpy as np
from THeSeuSS import InputsPreparation as inputs
from itertools import chain
import shutil, os
from collections import OrderedDict


class GenerateDisplacements():

    GEOMETRY_IN = 'vibrations/geometry.in'
    GEO_GEN = 'vibrations/geo.gen'

    def __init__(self, code: str):

        self.code = code
        self.path = os.getcwd()
        self.geometry_processor = None
        self.disp = None
        self.disp_inp = None
        self._set_geometry_processor()

    def _set_geometry_processor(self):
        """
        Setup the GeometryProcessor class.
        """

        if self.code == 'aims':
            geom_input = GenerateDisplacements.GEOMETRY_IN
        elif self.code == 'dftb+':
            geom_input = GenerateDisplacements.GEO_GEN

        self.geometry_processor = inputs.GeometryProcessor(geom_input, self.code)

    def _read_coordinates(self)-> list:
        """
        Returns the coordinates of the structure. 
        """
        
        return self.geometry_processor.read_coordinates()

    def _get_no_of_atoms(self)-> int:
        """
        Returns the number of atoms of the structure.
        """

        return self.geometry_processor.number_of_atoms()

    def _get_atom_type(self)-> tuple:
        """
        Returns both the atom type for each atom 
        and a corresponding numerical identifier for each atom type.
        """

        return self.geometry_processor.read_the_atom_type()

    def _define_the_displacement(self):
        """
        Defines the displacement value based on the chosen computational code ('aims' or 'dftb+').
        """

        if self.code == 'aims':
            self.disp = 0.01
        elif self.code == 'dftb+':
            self.disp = 0.005

    def _name_displaced_input(self, index: int):
        """
        Generate the input of the displaced structure.
        """

        if self.code == 'aims':
            self.disp_inp = f'geometry.in-{index:03d}'
        elif self.code == 'dftb+':
            self.disp_inp = f'geo.genS-{index:03d}'

    def _generate_displaced_input(self, coord: list):
        """
        Generates the geometry input file for both FHIaims and DFTB+,
        containing the displaced atomic coordinates.
        """
        
        no_of_atoms = self._get_no_of_atoms()
        numbers, types = self._get_atom_type()
        line = ''

        path_of_disp_str = os.path.join(self.path, 'vibrations', self.disp_inp) 
        with open(path_of_disp_str, 'w') as fh:
            if self.code == 'aims':
                for row, type_ in zip(coord, types):
                    coordinates_str = ' '.join([f'{coord}' for coord in row])
                    line += f'atom {coordinates_str} {type_}\n'
                fh.write(line)
            elif self.code == 'dftb+':
                fh.write(f'{no_of_atoms} C\n')
                list_types = list(OrderedDict.fromkeys(types))
                types_str = ' '.join(list_types)
                fh.write(f'{types_str}\n')
                for i, row, no_type in zip(range(1, no_of_atoms+1), coord, numbers):
                    coordinates_str = ' '.join([f'{coord}' for coord in row])
                    line += f' {i} {no_type} {coordinates_str}\n' 
                fh.write(line)

    def _create_supercell_file(self):
        """
        Creates a supercell input file based on the initial molecular geometry input file. 
        Since the concept of a supercell typically applies to periodic systems, 
        in the case of molecular structures, this file is created to adhere to PHONOPY's convention. 
        """

        if self.code == 'aims':
            geom_inp = 'geometry.in'
            geom_inp_supercell = 'geometry.in.supercell'
        elif self.code == 'dftb+':
            geom_inp = 'geo.gen'
            geom_inp_supercell = 'geo.genS'

        path_geom_inp = os.path.join(self.path, 'vibrations', geom_inp)
        path_geom_inp_supercell = os.path.join(self.path, 'vibrations', geom_inp_supercell)
        shutil.copy(path_geom_inp, path_geom_inp_supercell)

    def create_the_displaced_structures(self):
        """
        Generate displaced input structures by adding and subtracting the displacement value.
        """

        self._create_supercell_file()
        coordinates = self._read_coordinates()
        no_of_atoms = self._get_no_of_atoms()
        self._define_the_displacement()
        k = 1

        for i, row in enumerate(coordinates):
            for j, value in enumerate(row):
                coordinates[i][j] += self.disp
                self._name_displaced_input(k)
                self._generate_displaced_input(coordinates)

                coordinates = self._read_coordinates()
                coordinates[i][j] -= self.disp
                self._name_displaced_input(k+1)
                self._generate_displaced_input(coordinates)

                k += 2
                coordinates = self._read_coordinates()
