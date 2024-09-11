'''Eigenvectors and eigenvalues from PHONOPY code'''

#AUTHOR: Ariadni Boziki

import yaml
import math
import numpy as np
import scipy.linalg as spl
import Constants as c
import SubmitCalculations as submit
import InputsPreparation as inputs
from mendeleev import element


class VibrationalFrequencies():

    GEOMETRY_IN = 'geometry.in'
    GEO_GEN = 'geo.gen'

    def __init__(self, code: str, output_file_of_SCFSs: str, dispersion: bool):

        self.code = code
        self.output_file_of_SCFSs = output_file_of_SCFSs
        self.dispersion = dispersion
        self.cell_dims = None
        self.commands = None
        self.geometry_processor = None
        self.no_of_atoms = None
        self.phonopy_calculator = None
        self.no_of_atoms = None
        self.conts = None
        self.displacement = None
        self.dataset = None
        self.delta_forces = None
        self.mass = None
        self._set_geometry_processor()
        self._set_phonopy_calculator()
        self._displacement()

    def _set_geometry_processor(self):
        """
        Setup the GeometryProcessor class.
        """

        if self.code == 'aims':
            geom_input = VibrationalFrequencies.GEOMETRY_IN
        elif self.code == 'dftb+':
            geom_input = VibrationalFrequencies.GEO_GEN

        self.geometry_processor = inputs.GeometryProcessor(geom_input, self.code)

    def _set_phonopy_calculator(self):
        """
        Setup the PhonopyCalculator class.
        """

        self.phonopy_calculator = submit.PhonopyCalculator(self.code, self.cell_dims, self.output_file_of_SCFSs)

    def _get_number_of_atoms(self)-> int:
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

    def _forces_output_path(self, drct: str)-> str:
        """
        Returns the path to the output file containing forces for a specific calculation.
        """

        return self.phonopy_calculator.force_identifier(drct)

    def _get_forces(self, output_path: str)-> list :
        """
        Returns the forces from the FHIaims output and results.tag DFTB+ file 
        """

        return self.phonopy_calculator.read_forces_from_output(output_path)

    def _get_sorted_directories(self)-> list:
        """
        Returns the sorted directories in the current directory that contain the substring 'Coord', 
        based on the numeric value in their names.
        """

        return self.phonopy_calculator.sort_directories()

    def _displacement(self):
        """
        Set displacement value based on code.
        For 'aims', the displacement is set to 0.01, while for 'dftb+', it is set to 0.005.
        """

        if self.code == 'aims':
            self.displacement = 0.01
        elif self.code == 'dftb+':
            self.displacement = 0.005

    def _set_forces_dictionary(self):
        """
        Set up a dictionary containing forces data for each atom.
        """

        self.no_of_atoms = self._get_number_of_atoms()
        self.conts = self._get_sorted_directories()

        self.dataset = {'natom': self.no_of_atoms,
                'atoms': []}

        for i in self.conts:
            output_path = self._forces_output_path(i)
            self.forces = self._get_forces(output_path)

            atom_number = i.split('-')[2]  
            atom_axis = i.split('-')[3]
            if i[-1] == '+':
                disp = self.displacement
            elif i[-1] == '-':
                disp = -self.displacement

            entry = {
                    'atom_number': atom_number,
                    'atom_axis': atom_axis,
                    'disp': disp,
                    'forces': self.forces
            }
            self.dataset['atoms'].append(entry)

    def _group_dispplacement_forces_dict(self):
        """
        Calculates the difference between the forces applied to atoms along the same axis and with opposite displacements.
        """

        self._set_forces_dictionary()

        self.delta_forces = {'natom': self.no_of_atoms,
                'atoms': []}

        atoms = self.dataset['atoms']
        for i, atom_i in enumerate(atoms[:-1]):
            for j, atom_j in enumerate(atoms[i+1:], start=i+1):
                if atom_i['atom_number'] == atom_j['atom_number'] and atom_i['atom_axis'] == atom_j['atom_axis']:
                    if self.code == 'aims':
                        forces_pos = np.array(next(atom['forces'] for atom in (atom_i, atom_j) if atom['disp'] == 0.01))
                        forces_neg = np.array(next(atom['forces'] for atom in (atom_i, atom_j) if atom['disp'] == -0.01))
                    elif self.code == 'dftb+':
                        forces_pos = np.array(next(atom['forces'] for atom in (atom_i, atom_j) if atom['disp'] == 0.005))
                        forces_neg = np.array(next(atom['forces'] for atom in (atom_i, atom_j) if atom['disp'] == -0.005))
                    delta_forces_entry = {
                        'number': atom_i['atom_number'],
                        'axis': atom_i['atom_axis'],
                        'difference': forces_pos - forces_neg
                    }
                    self.delta_forces['atoms'].append(delta_forces_entry)

    def _mass_matrix(self, d_forces: dict):
        """
        Computes the mass matrix product based on the atomic weights of atoms involved in the force calculation.
        Calculates the reciprocal square root of the mass product for each pair of atoms.
        """

        numbers, types = self._get_atom_type()
        masses = []
        for el in types:
            current_element = element(el)
            current_element.atomic_weight
            masses.append(current_element.atomic_weight)

        num_total_coordinates = self.no_of_atoms * 3
        self.mass = np.zeros((num_total_coordinates, num_total_coordinates))

        for atom1 in d_forces['atoms']:
            for atom2 in d_forces['atoms']:
                atom1_index = int(atom1['number'])
                atom2_index = int(atom2['number'])
                mass_product = masses[atom1_index]*masses[atom2_index]
                i = int(atom1['number']) * 3
                j = int(atom2['number']) * 3
                self.mass[i:i+3, j:j+3] = mass_product
        self.mass = 1.0/np.sqrt(self.mass)

    def _generation_of_fc_matrix(self, dforces: dict)-> np.ndarray:
        """
        Generates the force constants matrix based on the differences in forces between pairs of atoms with opposite displacements.
        """

        num_total_coordinates = self.no_of_atoms * 3
        force_constants = np.zeros((num_total_coordinates, num_total_coordinates))

        for i, atom in enumerate(dforces['atoms']):
            difference = atom['difference']
            difference_row = difference.flatten()
            force_constants[i, :] = -difference_row

        return force_constants

    def _division_of_fc_by_displacement(self, force_constants: np.ndarray)-> np.ndarray:
        """
        Divides the force constants matrix by twice the displacement value. 
        """

        force_constants /= (2*self.displacement)

        return force_constants

    def _symmetrization_of_fc_matrix(self, force_constants: np.ndarray)-> np.ndarray:
        """
        Symmetrizes the force constants matrix by averaging the non-diagonal 
        elements with their corresponding symmetric elements. 
        """

        num_total_coordinates = self.no_of_atoms * 3
        symmetric_elements = []
        for i in range(num_total_coordinates):
            for j in range(i+1, num_total_coordinates):

                average = (force_constants[i, j] + force_constants[j, i]) / 2
                force_constants[i, j] = average
                force_constants[j, i] = average

        return force_constants

    def _mass_weighted_hessian(self)-> np.ndarray:
        """
        Calculates the mass-weighted Hessian matrix / Dynamical matrix.
        """
        
        self._group_dispplacement_forces_dict()
        self._mass_matrix(self.delta_forces)
        fc = self._generation_of_fc_matrix(self.delta_forces)
        fc = self._division_of_fc_by_displacement(fc)
        fc = self._symmetrization_of_fc_matrix(fc)

        dynamical_matrix = np.multiply(self.mass, fc)
        return dynamical_matrix

    def read_eig_vec_phonopy(self, hessian: np.ndarray)-> [np.ndarray, np.ndarray, np.ndarray]:
        """
        Returns eigenvalues, eigenvectors and in turn frequencies.
        """

        eigvals, eigvecs = np.linalg.eigh(hessian)

        #   Different method for the diagonalization. The results of both methods 
        #   are almost identical
        #   eigvals_test, eigvecs_test = spl.eigh(hessian_mass_weighted)

        frequencies = np.sqrt(np.abs(eigvals)) * np.sign(eigvals)

        if self.code == 'aims':
            conversion_factor_to_THz = 15.633302
            conversion_factor_to_cm_minus_1 = 15.633302*33.356

        if self.code == 'dftb+':
            conversion_factor_to_THz = 154.10794
            conversion_factor_to_cm_minus_1 = 154.10794*33.356

        frequencies_THz = frequencies * conversion_factor_to_THz
        frequencies_in_cm_minus_1 = frequencies * conversion_factor_to_cm_minus_1

        return eigvecs, eigvals, frequencies_in_cm_minus_1
