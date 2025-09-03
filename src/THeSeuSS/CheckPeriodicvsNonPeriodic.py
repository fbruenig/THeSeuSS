'''Verifies whether the simulation of spectra is intended for a periodic or non-periodic system.'''

#AUTHOR: Ariadni Boziki

import os
import re
import numpy as np
from THeSeuSS import GeometryInputConversion as geominput
from THeSeuSS import MapAtoms
from THeSeuSS import GenerateDisplacementsMolecules as dispsmolecules
from THeSeuSS import SubmitCalculations as submit
from THeSeuSS import CheckSuccessOutput as check
from THeSeuSS import EigenvectorsFrequenciesPHONOPY as eigenvec


class PeriodicvsNonPeriodic():

    def __init__(self, code: str, cell_dims: str, output_file: str, dispersion: bool, restart: bool, commands: str, functional: str = None, subsystem_size: str = None):

        self.code = code
        self.cell_dims = cell_dims
        self.output_file = output_file
        self.dispersion = dispersion
        self.restart = restart
        self.commands = commands
        self.functional = functional
        self.subsystem_size = subsystem_size
        self.non_periodic = True
        self.path = os.getcwd()
        self.geom_input = None
        self.spglib_processor = None
        self.phonopy_calculator = None
        self.geometry_conversion = None
        self.generated_displacements_mol = None
        self.check_calculator = None
        self.calculator = None
        self.vibrational_freq = None
        self.eigvecs = None
        self.eigvals = None
        self.frequencies_in_cm_minus_1 = None
        self._initialization_of_classes()
        self._specify_geometry_inputs()
        self.check_periodic_vs_non_periodic()

    def _initialization_of_classes(self):
        """
        Initializes several classes that are needed for different aspects of the processing and submission workflow.
        """
        
        self.spglib_processor = MapAtoms.spglibProcessor(self.code)
        self.phonopy_calculator = submit.PhonopyCalculator(self.code, self.cell_dims, self.output_file, self.dispersion, self.restart, self.commands, self.functional)
        self.geometry_conversion = geominput.GeometryConversion(self.code)
        self.generated_displacements_mol = dispsmolecules.GenerateDisplacements(self.code)
        self.check_calculator = check.CheckOutputSuccess(self.code, self.output_file, self.dispersion, self.functional)
        self.calculator = submit.Calculator(self.code, self.output_file, self.dispersion, self.restart, self.functional, self.commands, self.cell_dims)
        self.vibrational_freq = eigenvec.VibrationalFrequencies(self.code, self.output_file, self.dispersion, self.cell_dims, self.restart, self.commands, self.functional, self.subsystem_size)

    def _specify_geometry_inputs(self):
        """
        Define geometry input files.
        """

        if self.code == 'aims':
            self.geom_input = 'geometry.in'
        elif self.code == 'dftb+':
            self.geom_input = 'geo.gen'
        elif 'so3lr' in self.code:
            self.geom_input = 'so3lr.xyz'

    def check_periodic_vs_non_periodic(self)-> bool:
        """
        Checks if the system is periodic or non-periodic.
        """

        cif_file = [file for file in os.listdir(self.path) if file.endswith('.cif')]
        if cif_file:
            self.non_periodic = False
        else:
            with open(self.geom_input, 'r') as geometry_input:
                for line in geometry_input:
                    if self.code == 'aims' and 'lattice_vector' in line:
                        self.non_periodic = False
                        break
                    elif self.code == 'dftb+':
                        columns = line.split()
                        if len(columns) == 3:
                            self.non_periodic = False
                            break
                    elif 'so3lr' in self.code and 'Lattice' in line:
                        self.non_periodic = False
                        break
        return self.non_periodic

    def print_periodic_vs_non_periodic(self):
        """
        Prints if the system is periodic or non-periodic.
        """

        if self.non_periodic:
            print(f'THE SYSTEM IS NON PERIODIC')
        else:
            print(f'THE SYSTEM IS PERIODIC')
    
    def code_initialization(self):
        """
        Starts the code based on the periodicity of the system. 
        If the system is periodic, and a .cif file is provided, 
        the geometry is converted and the international space group number 
        of the experimental structure is determined. 
        If the system is non-periodic, no action is taken.
        """

        if not self.non_periodic:
            self.geometry_conversion.check_cif_file()
            if self.code == 'aims' or self.code == 'dftb+':
                initial_space_group = self.spglib_processor.get_international_space_group_number(self.path)
                print(f'INTERNATIONAL SPACE GROUP NUMBER OF THE EXPERIMENTAL STRUCTURE: {initial_space_group}\n')
        else:
            pass

    def decision_submission_cell(self):
        """
        Determines the appropriate approach for generating the displaced structures.
        If the system is periodic, submits calculations using the PhonopyCalculator class. 
        If the system is non-periodic, generates displaced structures using the GenerateDisplacements class.
        """

        if not self.non_periodic:
            self.phonopy_calculator.submit_phonopy_displacements()
        else:
            self.generated_displacements_mol.create_the_displaced_structures()

    def periodic_space_group_calc(self):
        """
        Calculate the space group of the periodic structure after cell and geometry optimization.
        """

        flag_exit = self.check_calculator.single_successful_output()
        if flag_exit:
            exit()
        else:
            if self.code == 'aims' or self.code == 'dftb+':
                if not self.non_periodic:
                    self.calculator.frozen_phonon_approximation_drct()
                    source_path = os.path.join(self.path, 'vibrations')
                    space_group_of_optimized_str = self.spglib_processor.get_international_space_group_number(source_path)
                    print(f'INTERNATIONAL SPACE GROUP NUMBER OF THE OPTIMIZED STRUCTURE: {space_group_of_optimized_str}\n')
                else:
                    self.calculator.frozen_phonon_approximation_drct()

    def dyn_mat_for_eigvec_eig_val_freq(self):
        """
        Returns eigenvalues, eigenvectors and in turn frequencies.
        """
        if self.code == 'so3lr-ana':
            hessian = np.loadtxt("energies_mw-hessian.txt", skiprows=2)
            hessian = (hessian+hessian.T)/2  # symmetrize the Hessian matrix
        else:
            if self.non_periodic:
                hessian = self.vibrational_freq._mass_weighted_hessian()

            else:
                hessian = self.phonopy_calculator.disp_forces_dataset_dyn_matrix()
                hessian = np.real(hessian)
        
        self.eigvecs, self.eigvals, self.frequencies_in_cm_minus_1 = self.vibrational_freq.read_eig_vec_phonopy(hessian)
        
    def eigenvec_eigenval_freq(self)-> [np.ndarray, np.ndarray, np.ndarray]:
        """
        If successful, it returns the eigenvectors, eigenvalues and frequencies. 
        If an error occurs, it prints an error message.
        """

        try:
            self.dyn_mat_for_eigvec_eig_val_freq()
            if not self.non_periodic:
                print(f'THE FORCE CONSTANTS HAVE BEEN CALCULATED USING PHONOPY\n')
                print(f'THE DYNAMICAL MATRIX HAS BEEN CALCULATED USING PHONOPY\n')
                print('*' * 150)
            else: 
                print(f'THE FORCE CONSTANTS HAVE BEEN CALCULATED\n')
                print(f'THE DYNAMICAL MATRIX HAS BEEN CALCULATED\n')
                print('*' * 150)
             
            print(f'FREQUENCIES (1/cm)')
            print(f'{self.frequencies_in_cm_minus_1}')
            
            num_negative_freqs = np.sum(self.frequencies_in_cm_minus_1 < 0)
            print(f'NUMBER OF NEGATIVE FREQUENCIES: {num_negative_freqs}')
            num_excluded_freqs = np.sum(self.frequencies_in_cm_minus_1 < 30)
            mask = self.frequencies_in_cm_minus_1 >= 30
            self.frequencies_in_cm_minus_1 = self.frequencies_in_cm_minus_1[mask]
            np.savetxt("Frequency.txt", self.frequencies_in_cm_minus_1)

            return self.eigvecs, self.eigvals, self.frequencies_in_cm_minus_1, num_excluded_freqs

        except Exception as e:
            print(f'FORCE CONSTANTS/DYNAMICAL MATRIX / AN ERROR WAS OCCURED: {e}')
            print('*' * 150)
