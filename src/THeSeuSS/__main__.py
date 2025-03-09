'''Main module of the package'''

#AUTHOR: Ariadni Boziki


import os
import subprocess
from THeSeuSS import InputOutputFiles as inputfiles
from THeSeuSS import FiniteDisplacements_phonopy as disp
from THeSeuSS import TwoPointCentralDifference as cendiff
from THeSeuSS import EigenvectorsFrequenciesPHONOPY as eigvec
from THeSeuSS import IRRamanSpectra as spectra
from THeSeuSS import InputsPreparation as inputs
from THeSeuSS import SubmitCalculations as submit
from THeSeuSS import MapAtoms
from THeSeuSS import PlotSpectra as plots
from THeSeuSS import CheckSuccessOutput as check
from THeSeuSS import CheckPeriodicvsNonPeriodic as pervsnonper
from THeSeuSS import Restart as rst
import numpy as np
import argparse
import re



def main(
        cell_geometry_optimization: bool, 
        functional: str, 
        kpoints: str, 
        eev: str, 
        rho: str, 
        etot: str, 
        forces: str, 
        sc_iter_limit: str, 
        species: str, 
        geometry: str, 
        energy: str, 
        steps: str, 
        pol_grid: str, 
        supercell: bool, 
        code: str, 
        output_file: str, 
        dispersion: bool, 
        spectra_calculation: bool, 
        files_preparation: bool, 
        cell_dims: str, 
        submission_cell: bool, 
        plot_bands: bool, 
        commands: str, 
        max_force_component: str, 
        max_steps: str, 
        SCC_tolerance: str, 
        max_SCC_iterations: str,
        broadening: str,
        fwhm: str,
        dispersion_type: str,
        restart: bool
):

    calculator = submit.Calculator(code, output_file, dispersion, restart, functional, commands, cell_dims)
    check_calculator = check.CheckOutputSuccess(code, output_file, dispersion, restart, functional)
    phonopy_calculator = submit.PhonopyCalculator(code, cell_dims, output_file, dispersion, restart, commands, functional)


    message = 'Input'
    border = '*' * (len(message) + 4)
    print(f'{border}\n* {message} *\n{border}')
    print(f'CELL AND GEOMETRY OPTIMIZATION: {cell_geometry_optimization}')
    print(f'DISPERSION: {dispersion}')
    print(f'GENERATION OF DISPLACED STRUCTURES: {submission_cell}')
    print(f'PREPARATION OF DIRECTORIES FOR FINITE DISPLACEMENT METHOD: {files_preparation}')
    print(f'CALCULATION OF SPECTRA: {spectra_calculation}')
    print(f'RESTART: {restart}')
    print(f'CODE: {code}')
    print(f'OUTPUT FILE: {output_file}\n')

    check_periodic_non_periodic = pervsnonper.PeriodicvsNonPeriodic(code, cell_dims, output_file, dispersion, restart, commands, functional)
    non_periodic = check_periodic_non_periodic.check_periodic_vs_non_periodic()
    check_periodic_non_periodic.print_periodic_vs_non_periodic()
    check_periodic_non_periodic.code_initialization()

    if not non_periodic:
        print(f'SUPERCELL: {supercell}')
        print(f'DIMENSIONS OF SUPERCELL: {cell_dims}')

    if cell_geometry_optimization:
        if code == 'aims' or code == 'dftb+':
            generator = inputfiles.InputsGenerator(code, kpoints, functional, eev, rho, etot, forces, sc_iter_limit, species,
                False, geometry, energy, steps, None, max_force_component, max_steps, SCC_tolerance, max_SCC_iterations, 
                output_file, dispersion, dispersion_type, restart)
            generator.file_exists()
            calculator.submit_job()
            check_periodic_non_periodic.periodic_space_group_calc()
        elif code == 'so3lr':
            calculator.submit_geometry_opt_so3lr()

    if submission_cell:
        check_periodic_non_periodic.decision_submission_cell()

    if files_preparation:
        files_prep = disp.FDSubdirectoriesGeneration(code, kpoints, functional, eev, rho, etot, forces, 
                sc_iter_limit, species, pol_grid, SCC_tolerance, max_SCC_iterations, output_file, dispersion, dispersion_type, restart)
        files_prep.iterate_over_files()

        if code == 'aims' or code == 'dftb+':
            FHIaims_calculator = submit.Calculator(code, output_file, dispersion, restart, functional, commands, cell_dims)
            FHIaims_calculator.submit_jobs_in_parallel()
        elif code == 'so3lr':
            calculator.energy_forces_so3lr()

    if plot_bands:
        phonopy_calculator.plot_band_structure()

    if restart:
        rst_calc = rst.RestartCalculation(code, output_file, dispersion, restart, functional, commands, dimensions)
        rst_calc.restart_calculations()

    if spectra_calculation:
        check_calculator.check_for_success_calc_before_spectra()
        eigvecs, eigvals, freq, no_neg_freqs = check_periodic_non_periodic.eigenvec_eigenval_freq()

        if not non_periodic:
            phonopy_calculator.animate_eigenvectors()
        else:
            pass

        if code == 'aims' or code == 'dftb+':
            central_diff = cendiff.TwoPointCentralDiff(code, output_file, dispersion, supercell, restart, functional)
            pol, cartesian_pol = central_diff.pol_cart_pol_processor()

            intensities_calculator = spectra.IntensityCalculator(code, eigvecs, cartesian_pol, pol, output_file, no_neg_freqs, restart, functional)
            IRintensity, Ramanactivity = intensities_calculator.spectra_calculation()

            plot_vibrational_spectra = plots.SpectraPlotter(freq, IRintensity, Ramanactivity, broadening, fwhm)
            plot_vibrational_spectra.plot_spectra()

        print('*' * 150) 

def run():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', dest='input_file_name', type=str)
    args = parser.parse_args()
    input_file = args.input_file_name

    keywords = {} 

    with open(input_file, 'r') as file1:
        lines = file1.readlines()

        patt = r'(\S+)\s+(.+)'

        for line in lines:
            match = re.match(patt, line)
            if match:
                key = match.group(1)
                value = match.group(2)


                if key in ('cell_geometry_optimization', 'supercell', 'dispersion', 'spectra_calculation', 'files_preparation', 'submission', 'plot_bands', 'restart'):
                    value = eval(value)
                if key in ('functional', 'eev', 'rho', 'etot', 'forces', 'sc_iter_limit', 'species', 'energy', 'geometry', 'steps', 'code', 'output_file', 'commands', 'max_force_component', 'max_steps', 'SCC_tolerance', 'max_SCC_iterations', 'broadening', 'dispersion_type'):
                    value = value
                if key in 'fwhm':
                    value = float(value)
                if key in ('kpoints', 'dimensions', 'pol_grid'):
                    value = ' '.join(value.split()[0:])

            keywords[key] = value


    cell_geometry_optimization = keywords.get('cell_geometry_optimization')
    supercell = keywords.get('supercell')
    dispersion = keywords.get('dispersion')
    spectra_calculation = keywords.get('spectra_calculation')
    files_preparation = keywords.get('files_preparation')
    submission_cell = keywords.get('submission')
    plot_bands = keywords.get('plot_bands')
    functional = keywords.get('functional')
    eev = keywords.get('eev')
    rho = keywords.get('rho')
    etot = keywords.get('etot')
    forces = keywords.get('forces')
    sc_iter_limit = keywords.get('sc_iter_limit')
    species = keywords.get('species')
    energy = keywords.get('energy')
    geometry = keywords.get('geometry')
    steps = keywords.get('steps')
    code = keywords.get('code')
    output_file = keywords.get('output_file')
    commands = keywords.get('commands')
    kpoints = keywords.get('kpoints')
    cell_dims = keywords.get('dimensions')
    pol_grid = keywords.get('pol_grid')
    max_force_component = keywords.get('max_force_component')
    max_steps = keywords.get('max_steps')
    SCC_tolerance = keywords.get('SCC_tolerance')
    max_SCC_iterations = keywords.get('max_SCC_iterations')
    broadening = keywords.get('broadening')
    fwhm = keywords.get('fwhm')
    dispersion_type = keywords.get('dispersion_type')
    restart = keywords.get('restart')

    main(cell_geometry_optimization, functional, kpoints, eev, rho, etot, forces, sc_iter_limit, species, geometry, energy, steps, pol_grid, supercell, code, output_file, dispersion, spectra_calculation, files_preparation, cell_dims, submission_cell, plot_bands, commands, max_force_component, max_steps, SCC_tolerance, max_SCC_iterations, broadening, fwhm, dispersion_type, restart)

if __name__== '__main__':
    run()
