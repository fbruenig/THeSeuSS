'''It takes care of equivalent atoms. It maps the polarizability and 
cartesian polarization of irreducible atoms back to reducible atoms.''' 

#AUTHOR: Ariadni Boziki

import os
import spglib
import numpy as np
from THeSeuSS import InputsPreparation as inputs
from THeSeuSS import Coordinates as coords


class spglibProcessor():

    GEOMETRY_IN =  "geometry.in"
    GEO_GEN = "geo.gen"

    def __init__(self, code: str):
        """
        Attributes
            code: str -> code name
        """

        self.code = code
        self.path = os.getcwd()
        self.geometry_processor = None
        self.dataset = None
        self.pol_extended = np.array([])
        self.cartesian_pol_extended = np.array([])
        self.pol = np.array([])
        self.cartesian_pol = np.array([])

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

        self.natoms = self._get_number_of_atoms()
        lattice = self._get_lattice()
        positions = self._get_coordinates()
        numbers, types = self._get_atom_type()
        cell = (lattice, positions, numbers)
        self.dataset = spglib.get_symmetry_dataset(cell)

    def get_international_space_group_number(self, source_path)-> int:
        """
        Returns the international space group number.
        """

        coordinates_conversion = coords.ConversionUnitsCoordinates(self.code)
        if (self.code == 'aims'):
            path_geom_inp = coordinates_conversion.cartesian_to_fractional_FHIaims(source_path, supercell=False)
                
        if (self.code == 'dftb+'):
            path_geom_inp = coordinates_conversion.cartesian_to_fractional_DFTBplus(source_path, supercell=False)

        self._set_geometry_processor(path_geom_inp)
        self._extract_structural_characteristics()

        space_group_number = self.dataset['number']

        return space_group_number

    def _set_path_of_vibrations(self)-> str:
        """
        Sets the path of the input file, while vibrations directory has been generated. 
        """

        if self.code == 'aims':
            path_input = os.path.join(self.path, 'vibrations', spglibProcessor.GEOMETRY_IN)
        elif self.code == 'dftb+':
            path_input = os.path.join(self.path, 'vibrations', spglibProcessor.GEO_GEN)
        return path_input

    def _equivalent_symmetry_atoms(self)-> np.ndarray:
        """
        Utilizes the spglib library to compute and return an array 
        mapping equivalent atoms within the crystal structure.
        """

        path_input_vibrations = self._set_path_of_vibrations()
        self._set_geometry_processor(path_input_vibrations)
        self._extract_structural_characteristics()

        equivalent_atoms = self.dataset['equivalent_atoms']

        return equivalent_atoms

    def _merge_arrays(self, pol: np.ndarray, cartesian_pol: np.ndarray, element_axis_coord: np.ndarray):
        """
        Extends the polarizability and Cartesian polarization arrays by appending two additional columns: 
        one containing the atom number and the other containing the corresponding coordinates.
        """

        element_axis_coord_float = element_axis_coord.astype(np.float64)

        self.pol_extended = np.concatenate((pol, element_axis_coord_float), axis=1)
        self.cartesian_pol_extended = np.concatenate((cartesian_pol, element_axis_coord_float), axis=1)

    def _check_similarities_of_coord(coord: list, no_of_line: int)-> bool:
        """
        Determines if two sets of coordinates are identical.
        It returns True if they are exactly the same.
        """
        
        line = coord[no_of_line]

        return all(x == line[0] for x in line[1:]) 

    def _map_properties_of_equivalent_atoms(self, pol: np.ndarray, cartesian_pol: np.ndarray, element_axis_coord: np.ndarray):
        """
        Returns the polarizability and Cartesian polarization arrays, 
        accounting for the mapping of reducible atoms to irreducible atoms (equivalent atoms). 
        This function also addresses scenarios where only one coordinate of an atom has been displaced.
        """

        self.pol = np.empty([0,8])
        self.cartesian_pol = np.empty([0,5])

        equivalent_atoms = self._equivalent_symmetry_atoms()

        self._merge_arrays(pol, cartesian_pol, element_axis_coord)

        for atom, equiv_atom in enumerate(equivalent_atoms):
            for values, values_cartesian in zip(self.pol_extended, self.cartesian_pol_extended):
                if values[6] == equiv_atom and values[7] == 1.0:
                    self.pol = np.append(self.pol,[values], axis = 0)
                if values_cartesian[3] == equiv_atom and values_cartesian[4] == 1.0:
                    self.cartesian_pol = np.append(self.cartesian_pol,[values_cartesian], axis = 0)

            for values, values_cartesian in zip(self.pol_extended, self.cartesian_pol_extended):
                if values[6] == equiv_atom and values[7] == 2.0:
                    self.pol = np.append(self.pol,[values], axis = 0)
                if values_cartesian[3] == equiv_atom and values_cartesian[4] == 2.0:
                    self.cartesian_pol = np.append(self.cartesian_pol,[values_cartesian], axis = 0)

            for values, values_cartesian in zip(self.pol_extended, self.cartesian_pol_extended):
                if values[6] == equiv_atom and values[7] == 3.0:
                    self.pol = np.append(self.pol,[values], axis = 0)
                if values_cartesian[3] == equiv_atom and values_cartesian[4] == 3.0:
                    self.cartesian_pol = np.append(self.cartesian_pol,[values_cartesian], axis = 0)

        length_pol = len(self.pol)

        path_input_vibrations = self._set_path_of_vibrations()
        self._set_geometry_processor(path_input_vibrations)
        self._extract_structural_characteristics()
        number_of_atoms = self.natoms

        if number_of_atoms != length_pol:
            for ii in range(number_of_atoms):
                count = (self.pol[:,6] == ii).sum()
                if count != 0 and count % 3 != 0:
                    check = self._check_similarities_of_coord(self.geometry_processor.read_coordinates(),ii)
                    if check and count == 1:
                        kk = 0
                        for jj in self.pol:
                            if jj[6] == ii:
                                line = kk	
                            kk += 1
                        self.pol = np.insert(self.pol, line, self.pol[line,:], 0)
                        self.pol[line+1, 7] = 2	
                        self.pol = np.insert(self.pol, line+1, self.pol[line+1,:], 0)
                        self.pol[line+2, 7] = 3	

            for ii in range(number_of_atoms):
                count = (self.cartesian_pol[:,3] == ii).sum()
                if count != 0 and count % 3 != 0:
                    check = self._check_similarities_of_coord(self.geometry_processor.read_coordinates(),ii)
                    if check and count == 1:
                        kk = 0
                        for jj in self.cartesian_pol:
                            if jj[3] == ii:
                                line = kk	
                            kk += 1
                        self.cartesian_pol = np.insert(self.cartesian_pol, line, self.cartesian_pol[line,:], 0)
                        self.cartesian_pol[line+1, 4] = 2	
                        self.cartesian_pol = np.insert(self.cartesian_pol, line+1, self.cartesian_pol[line+1,:], 0)
                        self.cartesian_pol[line+2, 4] = 3	

    def prop_equivalent_atoms(self, pol: np.ndarray, cartesian_pol: np.ndarray, element_axis_coord: np.ndarray)-> [np.ndarray, np.ndarray]:
        """
        Returns the polarizability and Cartesian polarization arrays after removing the last two columns, 
        which represent the atom number and coordinate (axis) of each atom.
        """

        self._map_properties_of_equivalent_atoms(pol, cartesian_pol, element_axis_coord)

        self.pol = np.delete(self.pol,np.s_[6:8],axis=1)
        self.cartesian_pol = np.delete(self.cartesian_pol,np.s_[3:5],axis=1)

        return self.pol, self.cartesian_pol
