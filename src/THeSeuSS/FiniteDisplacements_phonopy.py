"""Creation of the subdirectories of the finite displaced geometries."""

#AUTHOR: Ariadni Boziki

import os
import shutil
import glob
from pathlib import Path
import re
import difflib
import itertools as it
import math
from THeSeuSS import InputOutputFiles as inputfiles
from THeSeuSS import CheckPeriodicvsNonPeriodic as pervsnonper



class FDSubdirectoriesGeneration():

    def __init__(
        self, 
        code: str, 
        kpoints: str, 
        functional: str, 
        eev: str, 
        rho: str, 
        etot: str, 
        forces: str, 
        sc_iter_limit: str, 
        species: str, 
        pol_grid: str, 
        SCC_tolerance: str, 
        max_SCC_iterations: str, 
        output_file: str,
        dispersion: bool,
        dispersion_type: str,
        restart: bool
    ):

        self.code = code
        self.kpoints = kpoints
        self.functional = functional
        self.eev = eev
        self.rho = rho
        self.etot = etot
        self.forces = forces
        self.sc_iter_limit = sc_iter_limit
        self.species = species
        self.pol_grid = pol_grid
        self.SCC_tolerance = SCC_tolerance
        self.max_SCC_iterations = max_SCC_iterations
        self.output_file = output_file
        self.dispersion = dispersion
        self.dispersion_type = dispersion_type
        self.restart = restart
        self.path = os.getcwd()
        self.cell_dims = None
        self.commands = None
        self.axis = None
        self.sign = None
        self.no_of_element = None
        self.file1 = None
        self.file2 = None
        self.chunk = None
        self.generator = None
        self.pattern = None
        self.name = None
        self.input_file_name = None

    def _split_specific_lines(self, line: str)-> list:
        """
        Extracts individual records from a line containing coordinate data and returns them as a list.
        """
        
        if self.code == 'aims':
            line = line.split()[1:4]
        if self.code == 'dftb+':
            line = line.split()[2:5]
        if self.code == 'so3lr':
            line = line.split()[1:4]
        list_line_tmp = [float(i) for i in line]
        list_line = [float("{0:0.3f}".format(i)) for i in list_line_tmp]

        return list_line

    def _diff_lines_threshold(self, final_line_1: str, final_line_2: str):
        """
        Compares two lines, each containing three records. 
        If the difference between the first two records is zero 
        and the difference between the third record is 0.02 for FHIaims - 0.01 for DFTB+, 
        the function determines and returns the axis of the displaced coordinate ('x', 'y', or 'z').
        """

        list1 = self._split_specific_lines(final_line_1)
        list2 = self._split_specific_lines(final_line_2)

        diff_1_tmp = list1[0] - list2[0]
        diff_2_tmp = list1[1] - list2[1]
        diff_3_tmp = list1[2] - list2[2]

        diff_1 = abs(round(diff_1_tmp, 2))
        diff_2 = abs(round(diff_2_tmp, 2))
        diff_3 = abs(round(diff_3_tmp, 2))

        if self.code == 'aims':
            value_disp = 0.02
        if self.code == 'dftb+':
            value_disp = 0.01
        if self.code == 'so3lr':
            value_disp = 0.02

        if diff_1 == value_disp:
            self.axis = 'x'
        elif diff_2 == value_disp:
            self.axis = 'y'
        elif diff_3 == value_disp:
            self.axis = 'z'

    def _diff_lines(self)-> int:
        """
        It compares two files line by line and returns the line
        as a number if the two lines are different. 
        """

        no_of_line = 0
        final_line_1 = ''
        final_line_2 = ''
        final_no_of_line = 0

        for line_file1, line_file2 in zip(self.file1, self.file2):
            if line_file1 != line_file2:
                final_line_1 = line_file1
                final_line_2 = line_file2
                self._diff_lines_threshold(final_line_1, final_line_2)
                final_no_of_line = no_of_line
            no_of_line += 1
        return  final_no_of_line

    def _displaced_element(self, no_of_line: int):
        """
        In the FHIaims input file 'geometry.in', cell parameters are defined at the beginning, 
        with the atoms' coordinates starting from the 6th line onward. In the DFTB+ input file 'geo.gen', 
        initial lines contain comments and information about the number of atoms, 
        with the atoms' coordinates starting from the 3rd line.
        """

        check_periodic_non_periodic = pervsnonper.PeriodicvsNonPeriodic(self.code, self.cell_dims, self.output_file, self.dispersion, self.restart, self.commands, self.functional)
        non_periodic = check_periodic_non_periodic.check_periodic_vs_non_periodic()

        if not non_periodic:
            if self.code == 'aims':
                number = 5
            elif self.code == 'dftb+':
                number = 2
        else:
            if self.code == 'aims':
                number = 0
            elif self.code == 'dftb+':
                number = 2
        if self.code == 'so3lr':
            number = 2
        self.no_of_element = no_of_line - number

    def _create_dir(self, i: str, index: int):
        """
        Creation of folders.
        """
      
        if self.code == 'aims' or self.code == 'dftb+':
            os.chdir('vibrations')
        elif self.code == 'so3lr':
            pass
        os.mkdir(f'Coord-{index}-{self.no_of_element}-{self.axis}-{i}_{self.sign}')
        if self.code == 'aims' or self.code == 'dftb+':
            os.chdir('../')
        elif self.code == 'so3lr':
            pass

    def _rename(self, filename_path: str, name_path: str):
        """
        Rename files.
        """

        os.rename(filename_path, name_path)

    def _move_file(self, src_path: str, dst_path: str, filename: str):
        """
        Move files into directories.
        """

        src_path = os.path.join(src_path, filename)
        dst_path = os.path.join(dst_path, filename)
        shutil.move(src_path, dst_path)

    def _even_odd_numbers(self, num: int):
        """
        Returns the sign of the displacement.
        """

        if (num % 2) == 0:
            self.sign = '-'
        else:
            self.sign = '+'

    def _sorting_inputs(self):
        """
        Identifies files containing the specified filename
        within the current working directory, sorts them alphabetically,
        and groups them into pairs for further processing.
        """

        dir_path = os.path.join(self.path, 'vibrations')
        if self.code == 'aims' or self.code == 'dftb+':
            contents = os.listdir(dir_path)
        elif self.code == 'so3lr':
            contents = os.listdir(self.path)
        contents_list = []
        for i in contents:
            if self.input_file_name in i:
                contents_list.append(i)
        if self.code == 'aims' or self.code == 'dftb+':
            contents_list = sorted(contents_list, key=lambda x: int(re.search(r'\d+', x).group()))
        elif self.code == 'so3lr':
            contents_list = sorted(contents_list, key=lambda x: int(re.findall(r'\d+', x)[1]))
        self.chunk = [contents_list[x:x+2] for x in range(0, len(contents_list), 2)]

    def _drct_at_final_dest(self, index: int, namefile: str):
        """
        Creates directories with a specified name and moves files to their final destination directory 
        based on a specified pattern.
        """

        self._create_dir(namefile, index)
        created_dir = f'Coord-{index}-{self.no_of_element}-{self.axis}-{namefile}_{self.sign}'

        if self.code == 'aims' or self.code == 'dftb+':
            source_path = os.path.join(self.path, 'vibrations')
        if self.code == 'so3lr':
            source_path = self.path
        dst_path = os.path.join(source_path, created_dir)
        self._move_file(source_path, dst_path, namefile)

        namefile_path = os.path.join(dst_path, namefile)
        name_path = os.path.join(dst_path, self.name)
        self._rename(namefile_path, name_path)	

    def _define_namefiles_generate_drct(self):
        """
        Defines filenames based on the code.
        """
       
        if self.code == 'aims':
            self.pattern = 'geometry.in-*'
            self.name = 'geometry.in'
            self.input_file_name = 'geometry.in-'

        elif self.code == 'dftb+':
            self.pattern = 'geo.genS-*'
            self.name = 'geo.gen'
            self.input_file_name = 'geo.genS-'

        elif self.code == 'so3lr':
            self.pattern = 'so3lr-*'
            self.name = 'so3lr.xyz'
            self.input_file_name = 'so3lr-'

    def _displacements(self):
        """
        Identifies files with displacements on the same atom and coordinate, extracting the corresponding characteristics: 
        atom number, coordinate, and axis of displacement. 
        """

        index = 0
        self._define_namefiles_generate_drct()
        self._sorting_inputs()


        for item in self.chunk:
            if self.code == 'aims' or self.code == 'dftb+':
                path_file1 = os.path.join(self.path, 'vibrations', item[0])
                path_file2 = os.path.join(self.path, 'vibrations', item[1])
            elif self.code == 'so3lr':
                path_file1 = os.path.join(self.path, item[0])
                path_file2 = os.path.join(self.path, item[1])
            self.file1 = open(path_file1, 'r')
            self.file2 = open(path_file2, 'r')

            no_of_line =  self._diff_lines()
            if self.axis == None:
                continue
            self._displaced_element(no_of_line)
            no1 = re.sub('[^0-9]', '', item[0])
            no2 = re.sub('[^0-9]', '', item[1])
            no1 = int(no1)
            no2 = int(no2)

            if self.code == 'aims' or self.code == 'dftb+': 
                dir_path = os.path.join(self.path, 'vibrations')
            elif self.code == 'so3lr':
                dir_path = self.path
            files = Path(dir_path).glob(self.pattern)
            for file in files:
                f = str(file)
                namefile = f.split("/")[-1]
                number = re.sub('[^0-9]', '', namefile)
                number = int(number)
                if no1 == number:
                    self._even_odd_numbers(no1)
                    self._drct_at_final_dest(index, namefile)
                if no2 == number:
                    self._even_odd_numbers(no2)
                    self._drct_at_final_dest(index, namefile)
            index += 1

    def _copy_files_in_dir(self, filename: str):	
        """
        Copy input files in the generated directories.
        """
        
        if self.code == 'aims' or self.code == 'dftb+':
            dir_path = os.path.join(self.path, 'vibrations')
        elif self.code == 'so3lr':
            dir_path = self.path
        contents = os.listdir(dir_path)
        for item in contents:
            path_item = os.path.join(dir_path, item)
            if os.path.isdir(path_item):
                file_path = os.path.join(dir_path, filename)
                shutil.copy(file_path, path_item)	

    def _generate_input_prms(self):
        """
        Generates input files for single point calculations based on the specified code.
        If the code is 'aims', generates the FHI-aims control file. If the code is 'dftb+',
        generates the DFTB+ parameter input file.

        Raises:
            Exception: If an error occurs during file generation.
        """

        if self.code == 'aims':
            try:
                path_sub = os.path.join(self.path, 'vibrations', 'control.in')
                self.generator.FHIaims_control_file(path_sub)
                print(f'THE CONTROL FILE (input of FHI-aims) FOR SINGLE POINT CALCULATION HAS BEEN GENERATED\n')
            except Exception as e:
                print(f'THE CONTROL FILE (input of FHI-aims) FOR SINGLE POINT CALCULATION HAS NOT BEEN GENERATED / AN ERROR WAS OCCURED: {e}')
        elif self.code == 'dftb+':
            try:
                path_sub = os.path.join(self.path, 'vibrations', 'dftb_in.hsd')
                self.generator.DFTB_parameters_input_file(path_sub)
                print(f'THE dftb_in.hsd (input of DFTB+) FOR SINGLE POINT CALCULATION HAS BEEN GENERATED\n')
            except Exception as e:
                print(f'THE dftb_in.hsd FILE (input of DFTB+) FOR SINGLE POINT CALCULATION HAS NOT BEEN GENERATED / AN ERROR WAS OCCURED: {e}')
        elif self.code == 'so3lr':
            pass

    def _pol_drct_DFTB(self):
        """
        Generates directories for polarizability calculations based on DFTB+.
        Note: This function assumes a dispersion method is used.
        """

        self.generator = inputfiles.InputsGenerator(self.code, self.kpoints, self.functional, self.eev, self.rho, self.etot, self.forces, 
            self.sc_iter_limit, self.species, True, None, None, None, self.pol_grid, None, None, 
            self.SCC_tolerance, self.max_SCC_iterations, self.output_file, None, None, self.restart)

        dir_path = os.path.join(self.path, 'vibrations')

        contents = [item for item in os.listdir(dir_path) if os.path.isdir(os.path.join(dir_path, item))]
        for i in contents:
            if 'Coord' in i:
                drct = i
            path_pol = os.path.join(dir_path, drct, 'polarizability')
            os.mkdir(path_pol)
            path_drct_geom_inp = os.path.join(dir_path, drct, 'geo.gen')
            shutil.copy(path_drct_geom_inp, path_pol)
            path_pol_file = os.path.join(path_pol, 'dftb_in.hsd')
            self.generator.DFTB_parameters_input_file(path_pol_file)
        print(f'POLARIZABILITY DIRECTORIES HAVE BEEN GENERATED WITHIN THE DIRECTORIES THAT INCLUDE THE DISPLACED STRUCTURES SINCE A DISPERSION METHOD IS USED\n')

    def _pol_drct_FHIaims(self):
        """
        Generates directories for polarizability calculations based on FHIaims.
        Note: This function assumes a hybrid functional is used. THE DFPT calculation is performed with PBE functionali, instead.
        """

        self.generator = inputfiles.InputsGenerator(self.code, self.kpoints, 'pbe', self.eev, self.rho, self.etot, self.forces, 
            self.sc_iter_limit, self.species, False, None, None, None, self.pol_grid, None, None, 
            self.SCC_tolerance, self.max_SCC_iterations, self.output_file, self.dispersion, self.dispersion_type, self.restart)

        dir_path = os.path.join(self.path, 'vibrations')

        contents = [item for item in os.listdir(dir_path) if os.path.isdir(os.path.join(dir_path, item))]
        for i in contents:
            if 'Coord' in i:
                drct = i
            path_pol = os.path.join(dir_path, drct, 'polarizability')
            os.mkdir(path_pol)
            path_drct_geom_inp = os.path.join(dir_path, drct, 'geometry.in')
            shutil.copy(path_drct_geom_inp, path_pol)
            path_pol_file = os.path.join(path_pol, 'control.in')
            self.generator.FHIaims_control_file(path_pol_file)
        print(f'POLARIZABILITY DIRECTORIES HAVE BEEN GENERATED WITHIN THE DIRECTORIES THAT INCLUDE THE DISPLACED STRUCTURES SINCE A HYBRID FUNCTIONAL IS USED. FOR THE DFPT CALCULATION PBE IS USED INSTEAD.\n')

    def iterate_over_files(self):
        """
        Processes displaced structures and generates required input files accordingly.
        Copies necessary files into directories containing displaced structures.
        Note: Exception handling is implemented to catch any errors during the iteration.
        """

        if self.code == 'so3lr':
            check_periodic_non_periodic = pervsnonper.PeriodicvsNonPeriodic(self.code, self.cell_dims, self.output_file, self.dispersion, self.restart, self.commands, self.functional)
            non_periodic = check_periodic_non_periodic.check_periodic_vs_non_periodic()

        if self.code == 'aims' or self.code == 'dftb+': 
            self.generator = inputfiles.InputsGenerator(self.code, self.kpoints, self.functional, self.eev, self.rho, self.etot, self.forces, 
                self.sc_iter_limit, self.species, True, None, None, None, self.pol_grid, None, None, 
                self.SCC_tolerance, self.max_SCC_iterations, self.output_file, self.dispersion, self.dispersion_type, self.restart)

        if self.code == 'aims':
            try:
                self._displacements()
                self._generate_input_prms()
                self._copy_files_in_dir('control.in')

                if self.functional not in ['pbe', 'lda']:
                    self._pol_drct_FHIaims()

                print(f'THE DIRECTORIES THAT CONTAIN THE DISPLACED STRUCTURES HAVE BEEN GENERATED\n')
                print('*' * 150) 
            except Exception as e:
                print(f'THE DIRECTORIES THAT CONTAIN THE DISPLACED STRUCTURES HAVE NOT BEEN GENERATED / AN ERROR WAS OCCURED: {e}')
        elif self.code == 'dftb+':
            try:
                self._displacements()
                self._generate_input_prms()
                self._copy_files_in_dir('dftb_in.hsd')

                if self.dispersion:
                    self._pol_drct_DFTB()
    
                print(f'THE DIRECTORIES THAT CONTAIN THE DISPLACED STRUCTURES HAVE BEEN GENERATED\n')
                print('*' * 150)
            except Exception as e:
                print(f'THE DIRECTORIES THAT CONTAIN THE DISPLACED STRUCTURES HAVE NOT BEEN GENERATED / AN ERROR WAS OCCURED: {e}')
        #elif self.code == 'so3lr' and non_periodic:
        #    self._displacements()
        #elif self.code == 'so3lr' and not non_periodic:
        #    pass
        elif self.code == 'so3lr':
            self._displacements()
