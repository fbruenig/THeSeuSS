'''It checks if the single point calculations have been succesfully finished.'''

#AUTHOR: Ariadni Boziki

import os
import re


class CheckOutputSuccess():

    def __init__(self, code: str, output: str, dispersion: bool, restart: bool, functional: str = None):

        self.code = code
        self.output = output
        self.dispersion = dispersion
        self.restart = restart
        self.functional = functional
        self.path = os.getcwd()
        self.contents = None

    def single_successful_output(self)-> bool:
        """
        Checks if the cell and geometry optimization calculation has been successful.
        """

        if not os.path.exists(self.output):
            raise FileNotFoundError(f'The output file does not exist')
        else:
            with open(self.output) as f:
                if self.code == 'aims':
                    pattern = 'Have a nice day.'
                    if not pattern in f.read():
                        print(f'THE OPTIMIZATION CALCULATION WAS NOT SUCCESSFUL \ PROBLEM DURING THE FHIAIMS CALCULATION\n')
                        flag_exit = True
                    else:
                        print(f'THE OPTIMIZATION CALCULATION WAS SUCCESSFUL\n')
                        print("*" * 121)
                        flag_exit = False
                elif self.code == 'dftb+':
                    pattern = 'DFTB+ running times'
                    if not pattern in f.read():
                        print(f'THE OPTIMIZATION CALCULATION WAS NOT SUCCESSFUL \ PROBLEM DURING THE DFTB+ CALCULATION\n')
                        flag_exit = True
                    else:
                        print(f'THE OPTIMIZATION CALCULATION WAS SUCCESSFUL\n')
                        print("*" * 121)
                        flag_exit = False

        return flag_exit

    def _get_directories(self):
        """
        Returns the folders that exist under a path.
        """

        new_path = os.path.join(self.path, 'vibrations')
        self.contents = [item for item in os.listdir(new_path) if os.path.isdir(os.path.join(new_path, item))]

    def _raise_output_not_exist(self, path: str):
        """
        Raises a FileNotFoundError if the specified file does not exist in the directory.
        """

        if not os.path.exists(path):
            raise FileNotFoundError(f"OUTPUT FILE DOES NOT EXIST IN DIRECTORY.")

    def _raise_output_not_successful(self, path: str, pattern: str):
        """
        Raises a ValueError if the pattern is not found in the file.
        """

        with open(path, 'r') as f:
            if not pattern in f.read():
                raise ValueError(f"OUTPUT FILE DID NOT FINISH SUCCESSFULLY.")

    def _check_success(self, directory: str, pattern: str):
        """
        Checks if the file in the given directory exists and contains a pattern.
        FileNotFoundError: If the file does not exist in the directory.
        ValueError: If the file does not contain the pattern.
        """

        path_drct = os.path.join(self.path, 'vibrations', directory, self.output)
        if self.restart:
            if not os.path.exists(path_drct):
                return path_drct
            elif os.path.exists(path_drct): 
                with open(path_drct, 'r') as f:
                    if not pattern in f.read():
                        return path_drct
        else:
            self._raise_output_not_exist(path_drct)
            self._raise_output_not_successful(path_drct, pattern)

    def _check_success_pol(self, directory: str, pattern: str):
        """
        Checks if the file for polarizability calculation in the given directory exists and contains a specific pattern.
        FileNotFoundError: If the file does not exist in the directory.
        ValueError: If the file does not contain the specific pattern.
        """

        path_pol = os.path.join(self.path, 'vibrations', directory, 'polarizability', self.output)
        if self.restart:
            if not os.path.exists(path_pol):
                return path_pol
            elif os.path.exists(path_pol):
                with open(path_pol, 'r') as f:
                    if not pattern in f.read():
                        return path_pol
        else:
            self._raise_output_not_exist(path_pol)
            self._raise_output_not_successful(path_pol, pattern)

    def _successful_output(self, pattern: str): 
        """
        Checks if the single point calculations related to finite displacement method have been successful.
        """
        
        if self.restart:
            not_completed_calculations = []
        
        self._get_directories()
        for directory in self.contents:
            if directory.startswith('Coord'):
                if self.restart:
                    path_output = self._check_success(directory, pattern)
                    if path_output:
                        path_out = os.path.dirname(path_output)
                        not_completed_calculations.append(path_out)
                else:
                    self._check_success(directory, pattern)

        if self.restart:
            return not_completed_calculations

    def _successful_output_dispersion(self, pattern: str):
        """
        Checks if the single point calculations related to frozen phonon approximation as well as 
        polarizability calculations have been successful.
        """

        if self.restart:
            not_completed_calculations = []

        self._get_directories()
        for directory in self.contents:
            if directory.startswith('Coord'):
                if self.restart:
                    path_output = self._check_success(directory, pattern)
                    if path_output:
                        path_out = os.path.dirname(path_output)
                        not_completed_calculations.append(path_out)
                    path_output_pol = self._check_success_pol(directory, pattern)
                    if path_output_pol:
                        path_out_pol = os.path.dirname(path_output_pol)
                        not_completed_calculations.append(path_out_pol)
                else:    
                    self._check_success(directory, pattern)
                    self._check_success_pol(directory, pattern)

        if self.restart:
            return not_completed_calculations

    def check_for_success_calc_before_spectra(self):
        """
        Checks for successful completion of calculations before proceeding to spectra generation.
        FileNotFoundError: If the file does not exist.
        ValueError: If the file does not contain the expected patterns.
        """

        if self.code == 'aims' and self.functional in ['pbe', 'lda']:
            if self.restart:
                not_completed_calcs = self._successful_output('Have a nice day.')
            else:
                self._successful_output('Have a nice day.')
        if self.code == 'aims' and self.functional not in ['pbe', 'lda']:
            if self.restart:
                not_completed_calcs = self._successful_output_dispersion('Have a nice day.')
            else:
                self._successful_output_dispersion('Have a nice day.')
        if self.code == 'dftb+' and not self.dispersion:
            if self.restart:
                not_completed_calcs = self._successful_output('DFTB+ running times')
            else:
                self._successful_output('DFTB+ running times')
        if self.code == 'dftb+' and self.dispersion:
            if self.restart:
                not_completed_calcs = self._successful_output_dispersion('DFTB+ running times')
            else:
                self._successful_output_dispersion('DFTB+ running times')

        if self.restart:
            print(not_completed_calcs)
            return not_completed_calcs
