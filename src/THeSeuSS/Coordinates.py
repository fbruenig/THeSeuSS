'''It converts the coordinates from cartesian to fractional.'''

#AUTHOR: Ariadni Boziki

import shutil
import numpy as np
import InputsPreparation as inputs
import os


class ConversionUnitsCoordinates():

    PATH = os.getcwd()

    def __init__(self, code: str):
        """
        Attributes
            code: str -> code name
        """

        self.code = code

    def _set_geometry_processor(self, geom_input: str):
        """
        Setup the GeometryProcessor class.
        """

        self.geometry_processor = inputs.GeometryProcessor(geom_input, self.code)

    def _get_number_of_atoms(self)-> int:
        """
        Returns the number of atoms of the structure.
        """

        return self.geometry_processor.number_of_atoms()

    def _get_lattice(self)-> list:
        """
        Returns the lattice parameters, (9 elements matrix).
        """

        return self.geometry_processor.read_lattice()

    def _get_coordinates(self)-> list:
        """
        Returns the coordinates of the structure.
        """

        return self.geometry_processor.read_coordinates()

    def _get_atom_type(self)-> tuple:
        """
        Returns both the atom type for each atom 
        and a corresponding numerical identifier for each atom type.
        """

        return self.geometry_processor.read_the_atom_type()

    def _extract_structural_characteristics(self):
        """
        Extracts structural characteristics,
        including the cell, coordinates, atom type, and atom number.
        Additionally, initializes the spglib library.
        """

        self.no_atoms = self._get_number_of_atoms()
        self.lattice = self._get_lattice()
        self.positions = self._get_coordinates()
        self.numbers, self.types = self._get_atom_type()
    
    def _vector_inside_the_box(self, super_invmat: np.ndarray, m: np.ndarray, super_mat: np.ndarray)-> np.ndarray:
        """
        Moves the atoms into the simulation box.
        """

        real_vectors = np.dot(super_invmat,m)
        vectors_in_the_box = (real_vectors + 1e-5) % 1.0 - 1e-5
        dot_vectors_in_the_box = np.dot(super_mat,vectors_in_the_box)

        return dot_vectors_in_the_box

    def _lattice_inversion(self):
        """
        Lattice inversion.
        """

        self.lattice = np.array(self.lattice)
        self.super_mat = self.lattice.transpose()
        self.super_invmat = np.linalg.inv(self.super_mat)

    def _copy_files_in_the_same_directory(self, filename: str, new_name: str):
        """
        Copies a file within the same directory.
        """

        shutil.copy(filename, new_name)

    def cartesian_to_fractional_DFTBplus(self, path_of_current_directory, supercell: bool = False):
        """
        Converts cartesian coordinates to fractional - DFTB+.
        """

        if supercell:
            geo_input = 'geo.genS'
        else:
            geo_input = 'geo.gen'
        path_geo_input = os.path.join(path_of_current_directory, geo_input)
        with open(path_geo_input, 'r') as file:
            file_content = file.read()
            if 'F' in file_content:
                path_geo_gen = os.path.join(path_of_current_directory, 'geo.gen')
                path_geo_backup_gen = os.path.join(path_of_current_directory, 'geo_backup.gen')
                if os.path.isfile(path_geo_backup_gen):
                    print('BACKUP GEOMETRY FILE ALREADY EXISTS.')
                else:
                    self._copy_files_in_the_same_directory(path_geo_gen, path_geo_backup_gen)
            else:
                self._set_geometry_processor(path_geo_input)
                self._extract_structural_characteristics()

                geometry_input = open(path_geo_input, 'r')
                atom_types_dftb = geometry_input.readlines()[1]
                atom_types_dftb = atom_types_dftb.split()

                path_geo_gen = os.path.join(path_of_current_directory, 'geo.gen')
                path_geo_backup_gen = os.path.join(path_of_current_directory, 'geo_backup.gen')
                self._copy_files_in_the_same_directory(path_geo_gen, path_geo_backup_gen)

                self._lattice_inversion()

                output = os.path.join(path_of_current_directory, "geo.gen.F")
                with open(output, 'w') as ff:
                    ff.writelines('%s F\n' %(self.no_atoms))
                    for i in atom_types_dftb:
                        ff.writelines('%s ' '' %(i))
                    ff.writelines('\n')
                    for i, j, k in zip(self.positions, self.numbers, range(self.no_atoms)):

                        kk = k + 1

                        ii = np.array(i)
                        trans_conv = self._vector_inside_the_box(self.super_invmat, ii, self.super_mat)
                        trans_conv = np.dot(self.super_invmat,trans_conv)

                        ff.writelines(' ' ' %s ' ' %s ' ' %s ' ' %s ' ' %s \n' %(kk, j, trans_conv[0], trans_conv[1], trans_conv[2]))
                    ff.writelines(' ' ' 0.0 ' ' 0.0 ' ' 0.0\n')
                    ff.writelines(' ' ' %s ' ' %s ' ' %s\n' %(self.lattice[0][0], self.lattice[0][1], self.lattice[0][2]))
                    ff.writelines(' ' ' %s ' ' %s ' ' %s\n' %(self.lattice[1][0], self.lattice[1][1], self.lattice[1][2]))
                    ff.writelines(' ' ' %s ' ' %s ' ' %s\n' %(self.lattice[2][0], self.lattice[2][1], self.lattice[2][2]))

                os.rename(output, path_geo_gen)

        return path_geo_input

    def cartesian_to_fractional_FHIaims(self, path_of_current_directory, supercell: bool = False):
        """
        It converts cartesian coordinates to fractional - FHI-aims.'''
        """

        if supercell:
            geo_input = 'geometry.in.supercell'
        else:
            geo_input = 'geometry.in'
        path_geo_input = os.path.join(path_of_current_directory, geo_input)
        with open(path_geo_input, 'r') as file:
            file_content = file.read()
            if 'atom_frac' in file_content:
                path_geometry_in = os.path.join(path_of_current_directory, 'geometry.in')
                path_geometry_backup_in = os.path.join(path_of_current_directory, 'geometry_backup.in')
                if os.path.isfile(path_geometry_backup_in):
                    print('BACKUP GEOMETRY FILE ALREADY EXISTS.')
                else:
                    self._copy_files_in_the_same_directory(path_geometry_in, path_geometry_backup_in)
            else:
                self._set_geometry_processor(path_geo_input)
                self._extract_structural_characteristics()

                path_geometry_in = os.path.join(path_of_current_directory, 'geometry.in')
                path_geometry_backup_in = os.path.join(path_of_current_directory, 'geometry_backup.in')
                self._copy_files_in_the_same_directory(path_geometry_in, path_geometry_backup_in)

                self._lattice_inversion()

                output = os.path.join(path_of_current_directory, 'geometry.in.F')
                with open(output, 'w') as ff:
                    ff.writelines('\n')
                    ff.writelines(' lattice_vector   %s  %s  %s\n' %(self.lattice[0][0], self.lattice[0][1], self.lattice[0][2]))
                    ff.writelines(' lattice_vector   %s  %s  %s\n' %(self.lattice[1][0], self.lattice[1][1], self.lattice[1][2]))
                    ff.writelines(' lattice_vector   %s  %s  %s\n' %(self.lattice[2][0], self.lattice[2][1], self.lattice[2][2]))
                    ff.writelines('\n')
                    for i, j in zip(self.positions, self.types):

                        ii = np.array(i)
                        trans_conv = self._vector_inside_the_box(self.super_invmat, ii, self.super_mat)
                        trans_conv = np.dot(self.super_invmat,trans_conv)

                        ff.writelines(' atom_frac  %s  %s  %s  %s\n' %(trans_conv[0], trans_conv[1], trans_conv[2], j))

                os.rename(output, path_geometry_in)

        return path_geo_input
