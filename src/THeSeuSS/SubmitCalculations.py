'''Submission of FHI-aims and PHONOPY calculations'''

#AUTHOR: Ariadni Boziki

import os
import re
import subprocess
import shutil
import concurrent.futures
import numpy as np
import glob
from THeSeuSS import InputsPreparation as inputs
from THeSeuSS import FiniteDisplacements_phonopy as finitedisps
from THeSeuSS import EigenvectorsFrequenciesPHONOPY as eigenfreq
from THeSeuSS import Restart as restart
from THeSeuSS import GeometryInputConversion as inputconvert
from THeSeuSS import CheckPeriodicvsNonPeriodic as pervsnonper
import multiprocessing
from concurrent.futures import ThreadPoolExecutor
from phonopy import Phonopy
from phonopy.structure.atoms import PhonopyAtoms
from phonopy.interface.calculator import read_crystal_structure, write_crystal_structure, write_supercells_with_displacements, get_default_displacement_distance
from phonopy.file_IO import write_FORCE_SETS, write_FORCE_CONSTANTS
from phonopy.phonon.band_structure import get_band_qpoints_and_path_connections
from ase import Atoms
from ase.io import read, write
from so3lr import So3lrCalculator
from ase.optimize import FIRE, BFGS
from ase.filters import FrechetCellFilter


class PhonopyCalculator:
   
    def __init__(
        self, 
        code: str, 
        cell_dims: str, 
        output_file_of_SCFSs: str,
        dispersion: bool,
        restart: bool,
        commands: str,
        functional: str
    ):

        self.code = code
        self.cell_dims = cell_dims
        self.output_file_of_SCFSs = output_file_of_SCFSs
        self.dispersion = dispersion
        self.restart = restart
        self.commands = commands
        self.functional = functional
        self.path = os.getcwd()
        self.conts = []
        self.forces = None
        self.all_forces = None
        self.dataset = None
        self.dynamical_matrix = None
        self.eigvecs = None
        self.eigvals = None
        self.freq = None
        self.geometry_processor = None
        self.no_of_atoms = None
        self.forces_path = None
        self._setup_geometry_processor()

    def _read_crystal_structure(self):
        """
        Reads the crystal structure file based on the specified code.
        """

        if self.code == 'aims':
            return read_crystal_structure('vibrations/geometry.in', interface_mode='aims')
        elif self.code == 'dftb+':
            return read_crystal_structure('vibrations/geo.gen', interface_mode='dftbp')
        elif 'so3lr' in self.code:

            check_periodic_non_periodic = pervsnonper.PeriodicvsNonPeriodic(self.code, self.cell_dims, self.output_file_of_SCFSs, self.dispersion, self.restart, self.commands, self.functional)
            non_periodic = check_periodic_non_periodic.check_periodic_vs_non_periodic()
            
            atoms = read('so3lr.xyz')
            symbols = atoms.get_chemical_symbols()
            cell = atoms.get_cell()
            positions = atoms.get_positions()
            if non_periodic:
                unitcell = PhonopyAtoms(symbols = symbols,
                                positions = positions)
            else:
                unitcell = PhonopyAtoms(symbols = symbols,
                                        cell = cell,
                                        positions = positions)

            return unitcell

    def _setup_phonopy(self):
        """
        Sets up Phonopy for phonon calculations.
        Reads the crystal structure and extracts the unit cell and optional structure info. 
        Initializes a Phonopy object with the unit cell and supercell matrix. 
        Generates atomic displacements using default distances specific to the code.
        """

        if self.code == 'aims' or self.code == 'dftb+':
            unitcell, optional_structure_info = self._read_crystal_structure()
        elif 'so3lr' in self.code:
            unitcell = self._read_crystal_structure()
        cell_dims_int = [int(num) for num in self.cell_dims.split()]
        supercell_matrix = [[cell_dims_int[0], 0, 0], [0, cell_dims_int[1], 0], [0, 0, cell_dims_int[2]]]
        self.phonon = Phonopy(unitcell,supercell_matrix)
        if self.code == 'aims' or self.code == 'dftb+':
            default_displacement_for_code = get_default_displacement_distance(interface_mode=self.code)
        elif 'so3lr' in self.code:
            default_displacement_for_code = get_default_displacement_distance(interface_mode='aims')
        self.phonon.generate_displacements(distance=default_displacement_for_code)
        self.disps = self.phonon.displacements

    def _convert_aims_extxyz(self):

        files = glob.glob("geometry.in-*")
        aims_to_extxyz = inputconvert.GeometryConversion(self.code)

        for old_file in files:

            num_part = old_file.split('-')[-1]
            new_name = f"so3lr-{int(num_part):03d}.xyz"
            aims_to_extxyz.aims_geom_input_to_extxyz(old_file, new_name)
            os.remove(old_file)

    def supercell_disp_PHONOPY(self):
        """
        Retrieves the supercells with atomic displacements generated by Phonopy and the supercell. 
        Writes the supercells with displacements to output files based on the code.
        """

        self._setup_phonopy()
        supercells = self.phonon.supercells_with_displacements
        gen_supercell = self.phonon.supercell
        if 'so3lr' in self.code:
            write_supercells_with_displacements('aims', gen_supercell, supercells)
            self._convert_aims_extxyz()
        if self.code == 'aims':
            os.chdir('vibrations')
            write_supercells_with_displacements('aims', gen_supercell, supercells)
            os.chdir('../')
        elif self.code == 'dftb+':
            os.chdir('vibrations')
            write_supercells_with_displacements('dftbp', gen_supercell, supercells)
            os.chdir('../')

    def submit_phonopy_displacements(self):
        """
        Generate structures with displacements using PHONOPY. If successful, it prints a confirmation message. 
        If an error occurs it prints an error message.
        """

        try:
            self.supercell_disp_PHONOPY()
            print(f'PHONOPY GENERATED THE STRUCTURES WITH THE DISPLACEMENTS\n')
            print('*' * 150) 

        except Exception as e:
            print(f'GENERATION OF THE STRUCTURES WITH THE DISPLACEMENTS BY PHONOPY / AN ERROR WAS OCCURED: {e}')
            print('*' * 150) 

    def _setup_geometry_processor(self):
        """
        Initializes GeometryProcessor object.
        """

        if self.code == 'aims':
            geom_input = os.path.join(self.path, 'vibrations', 'geometry.in.supercell')
        elif self.code == 'dftb+':
            geom_input = os.path.join(self.path, 'vibrations', 'geo.genS')
        elif 'so3lr' in self.code:
            geom_input = os.path.join(self.path, 'so3lr.xyz')

        self.geometry_processor = inputs.GeometryProcessor(geom_input, self.code)

    def _number_of_atoms(self)-> int:
        """
        Returns the number of atoms in the supercell based on the code. Used for both periodic and non periodic systems.
        """

        return self.geometry_processor.number_of_atoms()

    def _get_coordinates(self)-> list:
        """
        Returns the coordinates of the supercell based on the code. Used for both periodic and non periodic systems.
        """

        return self.geometry_processor.read_coordinates()

    def read_forces_from_output(self, output_path: str):
        """
        Extracts the forces from the FHIaims output and results.tag DFTB+ file.
        """

        self.no_of_atoms = self._number_of_atoms()

        if self.code == 'aims':
            forces_identifier = 'Total atomic forces (unitary forces cleaned)'
        elif self.code == 'dftb+':
            forces_identifier = 'forces'

        f=open(output_path, 'r')
        found_unique_line = False

        forces_lines = []
        for line in f:
            if forces_identifier in line:
                found_unique_line = True
                continue

            if found_unique_line:
                forces_lines.append(line.strip())
                if len(forces_lines) >= self.no_of_atoms:
                    break

        pattern = r"[-+]?\d+\.\d+E[-+]?\d+"

        self.forces = []

        for item in forces_lines:
            numbers = re.findall(pattern, item)
            if len(numbers) >= 3:
                last_three_columns = [float(num) for num in numbers[-3:]]
                self.forces.append(last_three_columns)

        return self.forces

    def read_forces_from_files_ML(self, path_filename):
        """
        Reads the forces from the outputs of so3lr.
        """

        self.forces = []
        inside_forces = False
        
        with open(path_filename, 'r') as f:
            for line in f:
                line = line.strip()
                
                if line.startswith("Forces = ["):
                    inside_forces = True
                    line = line.replace("Forces = [", "").strip()
                    
                if inside_forces:
                    line = line.replace("[", "").replace("]", "")
                    
                    if line:
                        try:
                            self.forces.append(list(map(float, line.split())))
                        except ValueError:
                            pass
                        
                    if "]" in line:
                        inside_forces = False

        return self.forces

    def all_forces_from_files_ML(self):
        """
        Reads the forces from all the outputs of so3lr.
        """

        files = glob.glob("energies_forces-*")
        files = sorted(files, key=lambda x: int(x.split('-')[-1].split('.')[0]))
        self.all_forces = []
        
        for forces_file in files:
            path_filename = os.path.join(self.path, forces_file)
            if os.path.exists(path_filename):
                forces = self.read_forces_from_files_ML(path_filename)
                self.all_forces.append(forces)

        return self.all_forces

    def sort_directories(self)-> list:
        """
        Sorts directories in the current directory that contain the substring 'Coord', 
        based on the numeric value in their names.
        """

        if self.code == 'aims' or self.code == 'dftb+': 
            new_path = os.path.join(self.path, 'vibrations')
        elif 'so3lr' in self.code:
            new_path = self.path
        contents = [item for item in os.listdir(new_path) if os.path.isdir(os.path.join(new_path, item))]
        for ii in contents:
            if 'Coord' in ii:
                self.conts.append(ii)
        self.conts.sort(key=lambda f: int(''.join(filter(str.isdigit, f))))
    
        return self.conts

    def force_identifier(self, drct: str)-> str: 
        """
        Constructs the path to the output file containing forces for a specific calculation.
        """

        if self.code == 'aims':
            self.forces_path = os.path.join(self.path, 'vibrations', drct, self.output_file_of_SCFSs)
        elif self.code == 'dftb+':
            self.forces_path = os.path.join(self.path, 'vibrations', drct, 'results.tag')
        elif 'so3lr' in self.code:
            self.forces_path = os.path.join(self.path, drct, 'energies_forces.txt')

        return self.forces_path

    def disp_forces_dataset_dyn_matrix(self):
        """
        Generates a dataset for producing the force constants and calculating the dynamical matrix 
        based on displacements and forces obtained from FHIaims or DFTB+ output files.
        """

        self.no_of_atoms = self._number_of_atoms()
        self.dataset = {'natom': self.no_of_atoms,
                'first_atoms': []}

        self._setup_phonopy()

        if self.code == 'aims' or self.code == 'dftb+':
            self.conts = self.sort_directories()
            for j, jj in zip(self.disps,self.conts):

                self.forces_path = self.force_identifier(jj)
                self.forces = self.read_forces_from_output(self.forces_path)
                entry = {
                        'number': np.array([j[0]]),
                        'displacement': np.array(j[1:]),
                        'forces': np.array(self.forces)
                }
                self.dataset['first_atoms'].append(entry)
        elif 'so3lr' in self.code:
            self.all_forces = self.all_forces_from_files_ML()
            
            for j, jj in zip(self.disps,self.all_forces):
                entry = {
                        'number': np.array([j[0]]),
                        'displacement': np.array(j[1:]),
                        'forces': np.array(jj)
                }
                self.dataset['first_atoms'].append(entry)
        self.phonon.dataset = self.dataset
        self.phonon.produce_force_constants()

        q = [0.0,0.0,0.0]
        self.dynamical_matrix = self.phonon.get_dynamical_matrix_at_q(q)

        return self.dynamical_matrix

    def animate_eigenvectors(self):
        """
        Creates the xyz_jmol file to animate the vibrations using jmol package.
        """
        
        dyn_mat = self.disp_forces_dataset_dyn_matrix()
        q = [0.0,0.0,0.0]
        anime_type='jmol'
        band_index=0
        amplitude=5
        num_div=0
        shift=[0.0,0.0,0.0]
        filename='anime.xyz'
        self.phonon.write_animation(q, anime_type, band_index, amplitude, num_div, shift, filename)

    def plot_band_structure(self):
        """
        Plots the band structure using predefined paths in reciprocal space.
        Uses Phonopy to calculate and visualize the band structure.
        """

        path = [[[0, 0, 0], [0, 0.5, 0], [0, 0.5, 0.5]],
                [[0, 0, 0.5], [0, 0, 0], [-0.5, 0, 0.5], [-0.5, 0.5, 0.5]],
                [[0, 0.5, 0], [-0.5, 0.5, 0], [-0.5, 0, 0], [0, 0, 0]]]
        labels = ["$\\Gamma$", "Z", "D", "B", "$\\Gamma$", "A", "E", "Z", "$C_{2}$", "$Y_{2}$", "$\\Gamma$"]
        qpoints, connections = get_band_qpoints_and_path_connections(path, npoints=51)
        self.phonon.run_band_structure(qpoints, path_connections=connections, labels=labels)
        self.phonon.plot_band_structure().savefig("phonon_band_structure.pdf")
        self.phonon.write_yaml_band_structure()
        band_dict = self.phonon.get_band_structure_dict()



class Calculator:

    def __init__(self, code: str, output_file: str, dispersion: bool, restart: bool, functional: str = None, commands: str = None, cell_dims: str = None,
                 subsystem_size: str = None, subsystem_reference_id: str = None, optimize_only_subsystem: bool = False):
        
        self.code = code
        self.commands = commands
        self.output_file = output_file
        self.functional = functional
        self.dispersion = dispersion
        self.restart = restart
        self.cell_dims = cell_dims
        self.subsystem_size = subsystem_size
        if self.subsystem_size is not None:
            self.subsystem_indices = list(range(int(subsystem_size)))
        self.subsystem_reference_id = subsystem_reference_id
        self.optimize_only_subsystem = optimize_only_subsystem

        if code == 'aims' or code == 'dftb+':
            pass
        elif 'so3lr' in code:
            geo = read('so3lr.xyz')
            self.non_periodic = not geo.get_pbc().any()

            calculate_hessian = True if code == 'so3lr-ana' else False
            #output_intermediate_quantities = ['partial_charges'] if code == 'so3lr' else None
            output_intermediate_quantities = ['dipole_vec'] if code == 'so3lr' else None

            # Stresses are not yet implemented when calculating the hessian with so3lr
            # (but it is straightforward to implement)
            if not self.non_periodic:
                calc = So3lrCalculator(
                        lr_cutoff=12.0,
                        calculate_stress=True if code == 'so3lr' else False,
                        calculate_hessian=calculate_hessian,
                        dtype=np.float64,
                        output_intermediate_quantities=output_intermediate_quantities,
                        has_aux=True if output_intermediate_quantities is not None else False
                        )
            else:
                calc = So3lrCalculator(
                        lr_cutoff=100,
                        calculate_stress=False,
                        calculate_hessian=calculate_hessian,
                        dtype=np.float64,
                        output_intermediate_quantities=output_intermediate_quantities,
                        has_aux=True if output_intermediate_quantities is not None else False
                        )
            self.calc = calc
        else:
            raise ValueError(f"Unsupported code: {self.code}. Supported codes are 'so3lr', 'aims', and 'dftb+'.")

    def submit_job(self):
        """
        Submission of electronic structure calculations, (FHIaims and DFTB+).
        """

        number_of_cores = os.environ.get('number_of_cores')
        if number_of_cores is not None:
            number_of_cores = int(number_of_cores)
            num_threads = max(1, number_of_cores)  # Adjust the number of threads as needed
            print(f'NUMBER OF CORES ALLOCATED FOR GEOMETRY OPTIMIZATION: {num_threads}')

            parts = self.commands.split(';')
            first_part = ";".join(parts[:-1]).strip()
            last_part = parts[-1].strip()
        else:
            num_threads = multiprocessing.cpu_count()
            print(f'NUMBER OF CORES: {num_threads}')

        if number_of_cores is not None:
            if ";" not in self.commands: 
                command_tmp = f'srun --cpus-per-task 1 --ntasks {num_threads} {last_part}'
            else:
                command_tmp = f'{first_part}; srun --cpus-per-task 1 --ntasks {num_threads} {last_part}'
        else:
            command_tmp = self.commands
        
        subprocess.run(command_tmp, shell=True, executable='/bin/bash')

    def read_input_for_ML(self):
        """
        Reads the structural characteristics and set the calculator.
        """
        atoms = read('so3lr.xyz')
        positions = atoms.get_positions()
        symbols = atoms.get_chemical_symbols()
        if 'so3lr' in self.code:
            return atoms, self.calc
        else:
            raise ValueError(f"Unsupported code: {self.code}. Supported codes are 'so3lr', 'aims', and 'dftb+'.")

    def signle_point_periodic_read_input_for_ML(self, num_part):
        """
        Reads the structural characteristics and set the calculator.
        """

        filename = f"so3lr-{int(num_part):03d}.xyz"
        atoms = read(filename)
        positions = atoms.get_positions()
        symbols = atoms.get_chemical_symbols()
        cell = atoms.get_cell()
        pbc = atoms.set_pbc((True, True, True))
        atoms = Atoms(symbols=symbols, positions=positions, cell=cell, pbc=pbc)

        if 'so3lr' in self.code:
            return atoms, self.calc
        else:
            raise ValueError(f"Unsupported code: {self.code}. Supported codes are 'so3lr', 'aims', and 'dftb+'.")

    def submit_geometry_opt_so3lr(self):
        """
        Submission of geometry optimization with so3lr.
        """

        atoms, calc = self.read_input_for_ML()
        atoms.calc = calc

        # Apply constraints if specific atoms should be optimized
        if self.optimize_only_subsystem:
            assert self.subsystem_indices is not None, "subsytem must be provided via subsystem_size when optimize_only_subsystem is True."
            from ase.constraints import FixAtoms
            # Create a mask of atoms to fix (all atoms NOT in indices_to_optimize)
            fix_indices = [i for i in range(len(atoms)) if i not in self.subsystem_indices]
            if self.subsystem_reference_id is not None:
                fix_indices.append(int(self.subsystem_reference_id))
            constraint = FixAtoms(indices=fix_indices)
            atoms.set_constraint(constraint)

        if not self.non_periodic:
            ecf = FrechetCellFilter(atoms, hydrostatic_strain=False, constant_volume=False)
            optimizer = FIRE(ecf, logfile="optimization.log")
            optimizer.run(fmax=0.0005)
        elif self.optimize_only_subsystem:
            optimizer = BFGS(atoms, logfile="optimization.log")
            optimizer.run(fmax=0.0005, steps=10000)  # Adjusted force convergence criterion for subsystem optimization
        else:
            optimizer = FIRE(atoms, logfile="optimization.log")
            optimizer.run(fmax=0.0005)


        write("optimized_so3lr.xyz", atoms, format="extxyz")

        os.rename('so3lr.xyz', 'so3lr_initial_str.xyz')
        os.rename('optimized_so3lr.xyz', 'so3lr.xyz')

        print("Final Energy:", atoms.get_potential_energy())
        print("Final Forces:", atoms.get_forces())

    def submit_single_point_for_so3lr(self, num_part):
        """
        Submission of single point calculation with so3lr.
        """

        if self.non_periodic:
            atoms, calc = self.read_input_for_ML()
        else:
            atoms, calc = self.signle_point_periodic_read_input_for_ML(num_part)
        atoms.calc = calc
        calc.calculate(atoms)
        energy, forces, dipole_moment = calc.results['energy'], calc.results['forces'], calc.results['aux']['dipole_vec'][0]
        if self.subsystem_size is not None:
            forces = forces[:int(self.subsystem_size), :]

        # # Alternative approach to calculate dipole moment:
        # positions = atoms.get_positions()
        # energy, forces, charges = calc.results['energy'], calc.results['forces'], calc.results['aux']['partial_charges']
        # if self.subsystem_size is not None:
        #    forces = forces[:int(self.subsystem_size), :]
        #    positions = positions[:int(self.subsystem_size), :]
        #    charges = charges[:int(self.subsystem_size)]
        # dipole_moment = np.sum(charges[:,None]*positions,axis=0)

        if self.non_periodic:
            filename = f"energies_forces.txt"
        else:
            filename = f"energies_forces-{int(num_part):03d}.txt" 

        with open(filename, "w") as f:
            f.write(f"Energy = {energy:.8f}\n")
            f.write(f"Dipole = [{dipole_moment[0]:.8f} {dipole_moment[1]:.8f} {dipole_moment[2]:.8f}]\n")
            f.write("Forces = [")
            for i, row in enumerate(forces):
                if i == len(forces) - 1:
                    f.write(" [" + "  ".join(f"{x:14.8f}" for x in row) + " ]]")  # Last row, no extra \n
                else:
                    f.write(" [" + "  ".join(f"{x:14.8f}" for x in row) + " ]\n")

    def get_mass_matrix(self, atoms):
        """
        Computes the mass matrix product based on the atomic weights of atoms involved in the force calculation.
        Calculates the reciprocal square root of the mass product for each pair of atoms.
        """

        masses = atoms.get_masses()  # Get the masses of the atoms
        masses = np.repeat(masses, 3)  # Repeat masses for each Cartesian coordinate (x, y, z)

        return np.sqrt(masses[:,None] * masses[None,:])  # Outer product to create a mass matrix

    def energy_hessian_so3lr(self):
        """
        Submission of single point calculation with so3lr.
        """

        atoms, calc = self.read_input_for_ML()
        atoms.calc = calc

        energy = atoms.get_potential_energy()
        hessian = calc.get_property('hessian', atoms)

        no_of_atoms = len(atoms)
        mass_matrix = self.get_mass_matrix(atoms)
        if self.subsystem_size is not None:
            no_of_atoms = int(self.subsystem_size)
            # If subsystem indices are provided, only consider those atoms
            mass_matrix = mass_matrix[:(no_of_atoms*3), :(no_of_atoms*3)]
            hessian = hessian[:no_of_atoms,:,:no_of_atoms,:]

        hessian = hessian.reshape((no_of_atoms*3, no_of_atoms*3))/ mass_matrix
        filename = f"energies_mw-hessian.txt"
        np.savetxt(filename, hessian, header=f"Energy = {energy:.8f}\nMass-weighted hessian =", fmt='%.8f')

    def energy_forces_so3lr(self):
        """
        Submits multiple single point calculations with so3lr.
        """

        if self.non_periodic:
            path = os.getcwd()
            contents = [item for item in os.listdir(path) if os.path.isdir(os.path.join(path, item))]
            conts = []
            for ii in contents:
                if 'Coord' in ii:
                    conts.append(ii)
            conts.sort(key=lambda f: int(''.join(filter(str.isdigit, f))))
            for i in conts:
                os.chdir(i)
                num_part = None
                self.submit_single_point_for_so3lr(num_part)
                os.chdir('../')
        else:
            files = glob.glob("so3lr-*")
            for old_file in files:
                num_part = old_file.split('-')[-1]
                num_part = num_part.split('.')[0]
                self.submit_single_point_for_so3lr(num_part)

    def frozen_phonon_approximation_drct(self):
        """
        Creates a directory to store files associated with the frozen phonon approximation.
        """

        new_geometry_flag = False
        if (self.code == 'aims'):
            geometry = 'geometry.in'
            new_geometry = 'geometry.in.next_step'
        elif (self.code == 'dftb+'):
            geometry = 'geo.gen'
            new_geometry = f'{self.output_file}.gen'

        path = os.getcwd()
        os.mkdir('vibrations')
        output_path = os.path.join(path, 'vibrations')
        path_new_geom_input = os.path.join(path, new_geometry)

        if os.path.isfile(path_new_geom_input):
            geometry_path = os.path.join(path, new_geometry)
            new_geometry_flag = True
        else:
            geometry_path = os.path.join(path, geometry)
        shutil.copy(geometry_path, output_path)  

        if new_geometry_flag:
            geometry_path_rename = os.path.join(output_path, geometry)
            os.rename(geometry_path, geometry_path_rename) 

    def run_command(self, command: str)-> str:
        """
        Executes the given command in the shell environment and returns the output.
        """

        try:
            # Run the command and capture its output
            result = subprocess.run(command, shell=True, capture_output=True, text=True, check=True)
            # Process the result as needed
            return result.stdout
        except subprocess.CalledProcessError as e:
            # Handle command failure
            return f"Error running {command}: {e}"

    def submit_jobs_in_parallel(self):
        """
        Submits multiple jobs in parallel using a ThreadPoolExecutor.
        Note: Adjust the number of threads based on the available computational resources and the number of folders.
        """

        command_statement = []

        if self.restart:
            rst = restart.RestartCalculation(self.code, self.output_file, self.dispersion, self.restart, self.functional, self.commands, self.cell_dims)
            not_completed_calcs = rst.directory_non_completed_calculations()
            no_of_folders = len(not_completed_calcs)

        else:
            path = os.getcwd()
            new_path = os.path.join(path, 'vibrations')
            contents = [item for item in os.listdir(new_path) if os.path.isdir(os.path.join(new_path, item))]
            conts = []
            for ii in contents:
                if 'Coord' in ii:
                    conts.append(ii)
            conts.sort(key=lambda f: int(''.join(filter(str.isdigit, f))))
            no_of_folders = len(conts)

        number_of_cores = os.environ.get('number_of_cores')
        if number_of_cores is not None:
            number_of_cores = int(number_of_cores)
            num_threads = max(1, int(number_of_cores/no_of_folders))  # Adjust the number of threads as needed
            print(f'NUMBER OF CORES ALLOCATED PER SINGLE POINT CALCULATION ASSOCIATED TO FINITE DIFFERENCE METHOD: {num_threads}')

            parts = self.commands.split(';')
            first_part = ";".join(parts[:-1]).strip()
            last_part = parts[-1].strip()
        else:
            num_threads = multiprocessing.cpu_count()
            if num_threads <= no_of_folders:
                print(f'NUMBER OF CORES: {num_threads}')
            else:
                num_threads = max(1, int(num_threads/no_of_folders))  # Adjust the number of threads as needed
                print(f'NUMBER OF CORES ALLOCATED PER SINGLE POINT CALCULATION ASSOCIATED TO FINITE DIFFERENCE METHOD: {num_threads}')

        if self.restart:
            for i in not_completed_calcs:
                if number_of_cores is not None:
                    if ";" not in self.commands: 
                        command_tmp = f'cd {i}; srun --cpus-per-task 1 --ntasks {num_threads} {last_part}'
                    else:
                        command_tmp = f'cd {i}; {first_part}; srun --cpus-per-task 1 --ntasks {num_threads} {last_part}'
                else:
                    command_tmp = f'cd {i}; {self.commands}'
                command_statement.append(command_tmp)

        else:
            for i in conts:
                if number_of_cores is not None:
                    if ";" not in self.commands: 
                        command_tmp = f'cd vibrations; cd {i}; srun --cpus-per-task 1 --ntasks {num_threads} {last_part}'
                    else:
                        command_tmp = f'cd vibrations; cd {i}; {first_part}; srun --cpus-per-task 1 --ntasks {num_threads} {last_part}'
                else:
                    command_tmp = f'cd vibrations; cd {i}; {self.commands}'
                command_statement.append(command_tmp)

                if self.code == 'aims' and self.functional not in ['pbe', 'lda']:
                    if number_of_cores is not None:
                        command_tmp = f'cd vibrations; cd {i}; cd polarizability; {first_part}; srun --cpus-per-task 1 --ntasks {num_threads} {last_part}'
                    else:
                        command_tmp = f'cd vibrations; cd {i}; cd polarizability; {self.commands}'
                    command_statement.append(command_tmp)

                if self.dispersion and self.code == 'dftb+':
                    if number_of_cores is not None:
                        command_tmp = f'cd vibrations; cd {i}; cd polarizability; srun --cpus-per-task 1 --ntasks {num_threads} {last_part}'
                    else:
                        command_tmp = f'cd vibrations; cd {i}; cd polarizability; {self.commands}'
                    command_statement.append(command_tmp)
        if number_of_cores is not None:
            num_run = number_of_cores//num_threads
        else:
            num_run = multiprocessing.cpu_count()//num_threads
        with ThreadPoolExecutor(max_workers=num_run) as executor:
            executor.map(self.run_command, command_statement)
