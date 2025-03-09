'''Preparation of inputs for cell and geometry optimization / single point calculations'''

#AUTHOR: Ariadni Boziki

import os
import sys
from THeSeuSS import InputsPreparation as inputs
from THeSeuSS import Constants as c
from THeSeuSS import CheckPeriodicvsNonPeriodic as pervsnonper
import hsd
from collections import OrderedDict


class InputsGenerator:

    GEOMETRY_IN = 'geometry.in'
    GEO_GEN = 'geo.gen'
    SO3LR_GEN = 'so3lr.xyz'

    '''Generation of FHIaims and DFTB+ inputs.

    Attributes
    ----------

    code : str
        Code
    kpoints: str 
        K-points (Three values: they are converted to integer)
    functional: str
        Functional (FHIaims)
    eev: str
        Convergence criterion for the self-consistency cycle, based on the sum of eigenvalues (FHIaims)
    rho: str
        Convergence criterion for the self-consistency cycle, based on the charge density (FHIaims)
    etot: str
        Convergence criterion for the self-consistency cycle, based on the total energy (FHIaims)
    forces: str
        Convergence criterion for the self-consistency cycle, based on energy derivatives/forces (FHIaims)
    sc_iter_limit: str
        Maximum number of s.c.f. cycles before a calculation is considered and abandoned (FHIaims)
    species: str
        Species defaults settings; basis set, integration grids, accuracy of the Hartree potential (FHIaims)
    frequencies: bool
    geometry: str
        Maximum residual force component per atom, below which the geometry relaxation is considered converged (geometry relaxation) (FHIaims)
    energy: str
        Energy amount by which a relaxation step can move upwards and is still accepted (FHIaims)
    steps: str
        Maximum number of steps after which a structure optimization will be aborted (FHIaims)
    pol_grid: str
        Polarization grid
    max_force_component: str / converted to float
        Optimization is stopped, if the force component with the maximal absolute value goes below this value (DFTB+)
    max_steps: str
        Maximum number of steps after which the optimization should stop (DFTB+)
    SCC_tolerance: str
        Stopping criteria for the scc. Tolerance for the maximum difference in any charge between two scc cycles (DFTB+)
    max_SCC_iterations: str
        Maximal number of scc cycles to reach convergence (DFTB+)
    output_file: str
        Output file
    dispersion: bool
    path: str
        Path of the directory from which the code has been executed
    package_path: str
        Path of the directory where the code is located
    dispersion_type: str
        Dispersion type (MBD_NL or TS)
    '''

    def __init__(
        self, 
        code: str,
        kpoints: str = None,
        functional: str = None,
        eev: str = None,
        rho: str = None,
        etot: str = None,
        forces: str = None,
        sc_iter_limit: str = None,
        species: str = None,
        frequencies: str = None,
        geometry: str = None,
        energy: str = None,
        steps: str = None,
        pol_grid: str = None,
        max_force_component: str = None,
        max_steps: str = None,
        SCC_tolerance: str = None,
        max_SCC_iterations: str = None,
        output_file: str = None,
        dispersion: bool = False,
        dispersion_type: str = None,
        restart: bool = False
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
        self.frequencies = frequencies
        self.geometry = geometry
        self.energy = energy
        self.steps = steps
        self.pol_grid = pol_grid
        self.max_force_component = max_force_component
        self.max_steps = max_steps
        self.SCC_tolerance = SCC_tolerance
        self.max_SCC_iterations = max_SCC_iterations
        self.output_file = output_file
        self.dispersion = dispersion
        self.restart = restart
        self.dispersion_type = dispersion_type
        self.path = os.getcwd()
        self.species_path = os.path.abspath(__file__)
        self.package_path = os.path.dirname(self.species_path)
        self.input_parameters = {}
        self.cell_dims = None
        self.commands = None
        self.non_periodic = None
        self.drct_existance = None
        self._set_GeometryProcessor()
        self._check_periodicity()
        self._check_coord_directories_existence()

    def _check_periodicity(self):
        """
        Setup PeriodicvsNonPeriodic class.
        """

        check_periodic_non_periodic = pervsnonper.PeriodicvsNonPeriodic(self.code, self.cell_dims, self.output_file, self.dispersion, self.restart, self.commands, self.functional)
        self.non_periodic = check_periodic_non_periodic.check_periodic_vs_non_periodic()

    def _set_GeometryProcessor(self):
        """
        Setup the GeometryProcessor class.
        """

        if self.code == 'aims':
            geom_input = InputsGenerator.GEOMETRY_IN
        elif self.code == 'dftb+':
            geom_input = InputsGenerator.GEO_GEN
        elif self.code == 'so3lr':
            geom_input = InputsGenerator.SO3LR_GEN
        
        self.geometry_processor = inputs.GeometryProcessor(geom_input, self.code)

    def _get_atom_type(self)-> tuple:
        """
        Returns both the atom type for each atom 
        and a corresponding numerical identifier for each atom type.
        """

        return self.geometry_processor.read_the_atom_type()

    def _print_output(self):
        """
        Outputs the results of the code execution.
        """

        for key, value in self.input_parameters.items():
            print(f' {key} {value} ')
        print(f'\n')

    def _print_header(self, pattern: str):
        """
        Prints a section header for the input parameters.
        """

        width = len(pattern) + 4
        print('*' * width)
        print(f'* {pattern} *')
        print('*' * width)

    def _check_coord_directories_existence(self)-> bool:
        """
        Checks for the existence of directories starting with "Coord" in the current directory.
        It is used in order to decide in the control.in input file of FHIaims if the commands
        for the calculation of polarizability and dipole moment will be added. It is only for non_periodic systems.
        """

        if self.code == 'aims':
            path_dst = os.path.join(self.path, 'vibrations', 'geometry.in.supercell')
        elif self.code == 'dftb+':
            path_dst = os.path.join(self.path, 'vibrations', 'geo.genS')

        if os.path.exists(path_dst):
            self.drct_existance = True
        else:
            self.drct_existance = False

    def FHIaims_control_file(self, control_file):
        """
        Generates the control input file of FHIaims.
    
        number_atype: int, list
            The number of each atom within the crystal structure
        atype: str, list
            The atom type of each atom within the crystal structure
        unique_atype: str, list
            The atom types present in the crystal structure (appearing only once)
        input_parameters: dict
            It stores the parameters found in the control file and is utilized to display the input information in the output
        species_path_updated: str
            Directory path where the species defaults settings are located
        path_of_file_atype: str
            Path of the species defaults of an atom type
        source_data: str
            Data of the species defaults of an atom type
        """

        number_atype, atype = self._get_atom_type()
        unique_atype = list(OrderedDict.fromkeys(atype))

        with open(control_file, 'w') as fh:

            lines = [
                f'xc {self.functional}\n',
                f'\n',
                f'sc_accuracy_eev {self.eev}\n',
                f'sc_accuracy_rho {self.rho}\n',
                f'sc_accuracy_etot {self.etot}\n',
                f'sc_accuracy_forces {self.forces}\n',
                f'sc_iter_limit {self.sc_iter_limit}\n',
                f'\n'
                ]

            self.input_parameters = {
                'FUNCTIONAL:': self.functional,
                'CONVERGENCE CRITERION FOR THE SELF-CONSISTENCY CYCLE, BASED ON THE SUM OF EIGENVALUES:': self.eev,
                'CONVERGENCE CRITERION FOR THE SELF-CONSISTENCY CYCLE, BASED ON THE CHARGE DENSITY:': self.rho,
                'CONVERGENCE CRITERION FOR THE SELF-CONSISTENCY CYCLE, BASED ON THE TOTAL ENERGY:': self.etot,
                'CONVERGENCE CRITERION FOR THE SELF-CONSISTENCY CYCLE, BASED ON ENERGY DERIVATIVES ("FORCES"):': self.forces,
                'MAXIMUM NUMBER OF S.C.F. CYCLES BEFORE A CALCULATION IS CONSIDERED AND ABANDONED:': self.sc_iter_limit,
                'SPECIES DEFAULTS SETTINGS (BASIS SET, INTEGRATION GRIDS, ACCURACY OF THE HARTREE POTENTIAL):': self.species
                }

            if self.non_periodic:
                pass
            else:
                if self.functional in ['pbe', 'lda']:
                    lines.extend([
                        f'k_grid {self.kpoints}\n',
                        f'k_offset 0.5 0.5 0.5\n',
                        f'relativistic atomic_zora scalar\n',
                        f'\n'
                    ])
                else:
                    lines.extend([
                        f'k_grid {self.kpoints}\n',
                        f'relativistic atomic_zora scalar\n',
                        f'\n'
                    ])
                self.input_parameters['K-POINTS:'] = self.kpoints

            if self.dispersion:
                if self.dispersion_type == 'MBDNL':
                    lines.extend([
                        f'many_body_dispersion_nl\n'
                        ])
                elif self.dispersion_type == 'TS':
                    lines.extend([
                        f'vdw_correction_hirshfeld\n'
                        ])

            if self.geometry is not None:
                if self.non_periodic:
                    lines.extend([
                        f'relax_geometry bfgs {self.geometry}\n',
                        f'energy_tolerance {self.energy}\n',
                        f'max_relaxation_steps {self.steps}\n',
                        f'\n'
                    ])
                else:
                    lines.extend([
                        f'relax_geometry bfgs {self.geometry}\n',
                        f'energy_tolerance {self.energy}\n',
                        f'max_relaxation_steps {self.steps}\n',
                        f'relax_unit_cell full\n'
                        f'\n'
                    ])

                self.input_parameters['MAXIMUM RESIDUAL FORCE COMPONENT PER ATOM (in eV/A) BELOW WHICH THE GEOMETRY RELAXATION IS CONSIDERED CONVERGED (GEOMETRY RELAXATION):'] = self.geometry
                self.input_parameters['ENERGY AMOUNT BY WHICH A RELAXATION STEP CAN MOVE UPWARDS AND IS STILL ACCEPTED:'] = self.energy
                self.input_parameters['MAXIMUM NUMBER OF STEPS AFTER WHICH A STRUCTURE OPTIMIZATION WILL BE ABORTED:'] = self.steps

                self._print_header('Cell and geometry optimization input parameters')
                self._print_output()

            if self.frequencies:
                lines.extend([
                    f'compute_forces .true.\n',
                    f'final_forces_cleaned .true.\n'
                    f'KS_method serial\n',
                ])

                if self.non_periodic:
                    if self.drct_existance:    
                        lines.extend([
                            f'output dipole\n',
                            '\n'
                            ])

                    self._print_header('Single point calculation input parameters - Forces - Dipole moment - Polarizability')
                else:
                    if self.pol_grid is not None:
                        grid0 = self.pol_grid.split()[0]
                        grid1 = self.pol_grid.split()[1]
                        grid2 = self.pol_grid.split()[2]

                        lines.extend([
                            f'output polarization 1 {grid0} 1 1\n',
                            f'output polarization 2 1 {grid1} 1\n',
                            f'output polarization 3 1 1 {grid2}\n',
                            '\n'
                            ])

                    self.input_parameters['POLARIZATION GRID:'] = self.pol_grid
                    self._print_header('Single point calculation input parameters - Forces - Cartesian Polarization - Polarizability')

                self._print_output()

            if self.non_periodic:
                if self.drct_existance and self.functional in ['pbe', 'lda']:    
                    lines.extend([
                        f'DFPT polarizability\n',
                        '\n'
                        ])
                else:
                    pass
            else:
                if self.pol_grid is not None and self.functional in ['pbe', 'lda']:
                    lines.extend([
                        f'DFPT dielectric\n',
                        '\n'
                        ])
                else:
                    pass

            fh.writelines(lines)

        species_path_updated = os.path.join(self.package_path, 'species_defaults', self.species)

        for i in unique_atype:
            for filename in os.listdir(species_path_updated):
                species_type = filename.split('_')[1]
                if (i==species_type):
                    path_of_file_atype = os.path.join(species_path_updated, filename)
                    with open(path_of_file_atype, 'r') as source_file:
                        source_data = source_file.read()
                    with open(control_file, 'a') as target_file:
                        target_file.write(source_data)

    def _dict_to_hsd(self, data: dict)-> str:
        """
        Returns the necessary information in the correct format for DFTB+ to read the data specified in GenFormat/Geometry.
        """

        hsd = ""
        for key, value in data.items():
            if isinstance(value, dict):
                if key == "Geometry":            
                    hsd += f"{key} {{"
                    hsd += self._dict_to_hsd(value)
            else:
                if key == "GenFormat":
                    hsd += f"{key} {{\n"
                    hsd += f" <<<'{value}'\n"
                    hsd += "}}\n"
        return hsd

    def DFTB_parameters_input_file(self, dftb_in_file):
        """
        Generates the dftb_in.hsd input file of DFTB+

        species_path_updated: str
            Directory path where the .skf files are located
        data: dict
            It stores the parameters found in the dftb_in.hsd file
        """

        species_path_updated = os.path.join(self.package_path, '3ob-3-1')

        number_atom_type, atom_type = self._get_atom_type()
        atype_unique_list = list(OrderedDict.fromkeys(atom_type))

        costants_initialization = c.ConstantsSpectra(self.code)
        hubbard_derivates = costants_initialization.hubbard_derivates_dict()
        max_angular_mom = costants_initialization.max_angular_mom_dict()

        if not self.non_periodic:
            kpoint0 = self.kpoints.split()[0]
            kpoint1 = self.kpoints.split()[1]
            kpoint2 = self.kpoints.split()[2]

            value0 = 1 if int(kpoint0) % 2 == 0 else 0
            value1 = 1 if int(kpoint1) % 2 == 0 else 0
            value2 = 1 if int(kpoint2) % 2 == 0 else 0

        if self.max_force_component is not None:

            data = {
                'Geometry': {
                    "GenFormat": "geo.gen"
                },
                'Driver': 
                    {'ConjugateGradient': 
                        {'MovedAtoms': '1:-1', 
                        'MaxForceComponent': float(self.max_force_component), 
                        'MaxSteps': self.max_steps, 
                        'OutputPrefix': f'"{self.output_file}"'}}, 
                'Hamiltonian': 
                    {'DFTB': {'SCC': True, 
                        'SCCTolerance': self.SCC_tolerance, 
                        'MaxSCCIterations': self.max_SCC_iterations, 
                        'SlaterKosterFiles': {'Type2FileNames': 
                            {'Prefix': f'"{species_path_updated}/"', 'Separator': '"-"', 'Suffix': '".skf"'}}, 
                        'MaxAngularMomentum': {i: f'"{max_angular_mom[i]}"' for i in atype_unique_list},
                        'ThirdOrderFull': True, 
                        'HCorrection': {'Damping': {'Exponent': 4.05}}, 
                        'HubbardDerivs': {i: hubbard_derivates[i] for i in atype_unique_list}}}, 
                'Options': {'WriteResultsTag': True}, 
                'Analysis': {'PrintForces': True}} 

            self.input_parameters = {
                'OPTIMIZATION IS STOPPED, IF THE FORCE COMPONENT WITH THE MAXIMAL ABSOLUTE VALUE GOES BELOW:': self.max_force_component,
                'MAXIMUM NUMBER OF STEPS AFTER WHICH THE OPTIMIZATION SHOULD STOP:': self.max_steps,
                'STOPPING CRITERIA FOR THE SCC. TOLERANCE FOR THE MAXIMUM DIFFERENCE IN ANY CHARGE BETWEEN TWO SCC CYCLES:': self.SCC_tolerance,
                'MAXIMAL NUMBER OF SCC CYCLES TO REACH CONVERGENCE:': self.max_SCC_iterations
                    }

            if self.dispersion:
                if self.dispersion_type == 'MBD':
                    data['Hamiltonian']['DFTB']['Dispersion'] = {'MBD': {'KGrid': [1, 1, 1], 'Beta': 0.83}}
                elif self.dispersion_type == 'TS':
                    data['Hamiltonian']['DFTB']['Dispersion'] = {'TS': {'RangeSeparation': 0.94}}

                self.input_parameters['DISPERSION'] = self.dispersion_type

            if not self.non_periodic:
                data['Hamiltonian']['DFTB']['KPointsAndWeights'] = {'SuperCellFolding': [[kpoint0, 0, 0], [0, kpoint1, 0], [0, 0, kpoint2], [value0, value1, value2]]}
                data['Driver']['ConjugateGradient']['LatticeOpt'] = True
                self.input_parameters['K-POINTS:'] = self.kpoints
                
            self._print_header('Optimization input parameters')
            self._print_output()

        else:
            data = {
                'Geometry': {
                    "GenFormat": "geo.gen"
                },
                'Hamiltonian': 
                    {'DFTB': {'SCC': True, 
                        'SCCTolerance': self.SCC_tolerance, 
                        'MaxSCCIterations': self.max_SCC_iterations, 
                        'SlaterKosterFiles': {'Type2FileNames': 
                            {'Prefix': f'"{species_path_updated}/"', 'Separator': '"-"', 'Suffix': '".skf"'}}, 
                        'MaxAngularMomentum': {i: f'"{max_angular_mom[i]}"' for i in atype_unique_list},
                        'ThirdOrderFull': True, 
                        'HCorrection': {'Damping': {'Exponent': 4.05}}, 
                        'HubbardDerivs': {i: hubbard_derivates[i] for i in atype_unique_list}}}, 
                'Options': {'WriteResultsTag': True}, 
                'Analysis': {'PrintForces': True, 'Polarisability': {}}}

            self.input_parameters = {
                'STOPPING CRITERIA FOR THE SCC. TOLERANCE FOR THE MAXIMUM DIFFERENCE IN ANY CHARGE BETWEEN TWO SCC CYCLES:': self.SCC_tolerance,
                'MAXIMAL NUMBER OF SCC CYCLES TO REACH CONVERGENCE:': self.max_SCC_iterations
                    }

            if self.dispersion:
                if self.dispersion_type == 'MBD':
                    data['Hamiltonian']['DFTB']['Dispersion'] = {'MBD': {'KGrid': [1, 1, 1], 'Beta': 0.83}}
                elif self.dispersion_type == 'TS':
                    data['Hamiltonian']['DFTB']['Dispersion'] = {'TS': {'RangeSeparation': 0.94}}
                data['Analysis'] = {'PrintForces': True}
                self.input_parameters['DISPERSION'] = self.dispersion_type

            if not self.non_periodic:
                data['Hamiltonian']['DFTB']['KPointsAndWeights'] = {'SuperCellFolding': [[kpoint0, 0, 0], [0, kpoint1, 0], [0, 0, kpoint2], [value0, value1, value2]]}
                self.input_parameters['K-POINTS:'] = self.kpoints

            self._print_header('Single point calculation')
            self._print_output()

        hsd_content = self._dict_to_hsd(data)
        data_without_genformat = data.copy()
        data_without_genformat["Geometry"].pop("GenFormat")
        data_without_genformat.pop("Geometry")

        hsd.dump(data_without_genformat, "output1.hsd")
        with open(dftb_in_file, "w") as hsd_file:
            hsd_file.write(hsd_content)

        with open("output1.hsd", "r") as hsd_file1:
            content_to_append = hsd_file1.read()

        with open(dftb_in_file, "a") as hsd_file:
            hsd_file.write(content_to_append)

        os.remove("output1.hsd")

    def file_exists(self):
        """
        Checks the existence of input files required for cell and geometry optimization in FHI-aims or DFTB+ calculations. 
        If the file does not exist, it generates the required input file.
        """

        if self.code == 'aims':
            if os.path.isfile('control.in'):
                print(f'THE CONTROL FILE (input of FHI-aims) FOR OPTIMIZATION ALREADY EXISTS\n')
                pass
            else:
                try:
                    path_sub = os.path.join(self.path, 'control.in')
                    self.FHIaims_control_file(path_sub)
                    print(f'THE CONTROL FILE (input of FHI-aims) FOR OPTIMIZATION HAS BEEN GENERATED\n')
                except Exception as e:
                    print(f'THE CONTROL FILE GENERATION WAS NOT SUCCESSFUL / AN ERROR WAS OCCURED: {e}')

        elif self.code == 'dftb+':
            if os.path.isfile('dftb_in.hsd'):
                print(f'THE dftb_in.hsd FILE (input of DFTB+) FOR OPTIMIZATION ALREADY EXISTS\n')
                pass
            else:
                try:
                    path_sub = os.path.join(self.path, 'dftb_in.hsd')
                    self.DFTB_parameters_input_file(path_sub)
                    print(f'THE dftb_in.hsd FILE (input of DFTB+) FOR OPTIMIZATION HAS BEEN GENERATED\n')
                except Exception as e:
                    print(f'THE dftb_in.hsd FILE GENERATION WAS NOT SUCCESSFUL / AN ERROR WAS OCCURED: {e}')
