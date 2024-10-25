#AUTHOR: Ariadni Boziki


import unittest
from unittest.mock import MagicMock, patch, mock_open
import os
import re
from pathlib import Path
from THeSeuSS import FiniteDisplacements_phonopy  as finitedisps


class Test_FDSubdirectoriesGeneration(unittest.TestCase):

    @patch("builtins.open", new_callable=mock_open)
    @patch("os.listdir", return_value=["geo.genS-001", "geo.genS-002"])
    @patch("os.path.join", side_effect=lambda *args: "/".join(args))
    @patch("os.mkdir")
    @patch("os.rename")
    @patch("shutil.move")
    @patch("pathlib.Path.glob", return_value=[Path("geo.genS-001"), Path("geo.genS-002")])
    @patch.object(finitedisps.FDSubdirectoriesGeneration, "_define_namefiles_generate_drct")
    @patch.object(finitedisps.FDSubdirectoriesGeneration, "_diff_lines", return_value=6)
    @patch.object(finitedisps.FDSubdirectoriesGeneration, "_displaced_element")
    @patch.object(finitedisps.FDSubdirectoriesGeneration, "_even_odd_numbers")
    @patch.object(finitedisps.FDSubdirectoriesGeneration, "_drct_at_final_dest")
    def test_displacements(self, mock_drct_dest, mock_even_odd, mock_displaced_element, mock_diff_lines, mock_define_files, mock_glob, mock_move, mock_rename, mock_mkdir, mock_join, mock_listdir, mock_open):

        fd_gen = finitedisps.FDSubdirectoriesGeneration(code='dftb+', kpoints='3 3 3', SCC_tolerance='1e-3', max_SCC_iterations='100', output_file='output.txt', dispersion=True, dispersion_type='MBD', functional=None, eev=None, rho=None, etot=None, forces=None, sc_iter_limit=None, species=None, pol_grid=None)

        fd_gen.path = '/current/path'
        fd_gen.axis = 'x'
        fd_gen.no_of_element = 0
        fd_gen.sign = '+'
        fd_gen.pattern = 'geo.genS-*'
        fd_gen.input_file_name = 'geo.genS-'
        fd_gen.chunk = [['geo.genS-001', 'geo.genS-002']]

        fd_gen._displacements()

        mock_define_files.assert_called_once()
        mock_diff_lines.assert_called_once()
        mock_displaced_element.assert_called_once_with(6)
        mock_even_odd.assert_any_call(1)
        mock_even_odd.assert_any_call(2)

        mock_drct_dest.assert_any_call(0, "geo.genS-001")
        mock_drct_dest.assert_any_call(0, "geo.genS-002")

    @patch('os.listdir', return_value=['Coord-0-0-x-geometry.in-001_+', 'Coord-0-0-x-geometry.in-002_-'])
    @patch('os.path.isdir', return_value=True)
    @patch('os.path.join', side_effect=lambda *args: '/'.join(args))
    @patch('shutil.copy')
    @patch('os.mkdir')
    @patch('THeSeuSS.InputOutputFiles.InputsGenerator')
    def test_iterate_over_files_aims(self, mock_generator, mock_mkdir, mock_copy, mock_path_join, mock_isdir, mock_listdir):
        
        generator_instance = MagicMock()
        mock_generator.return_value = generator_instance

        fd_gen = finitedisps.FDSubdirectoriesGeneration(code='aims', kpoints='1 1 1', functional='pbe', eev='1E-5', rho='1E-7', etot='1E-6', forces='1E-4', sc_iter_limit='300', species='tight', pol_grid='10 10 10', SCC_tolerance=None, max_SCC_iterations=None, output_file='output.txt', dispersion=True, dispersion_type='TS')

        fd_gen._displacements = MagicMock()
        fd_gen._generate_input_prms = MagicMock()
        fd_gen._copy_files_in_dir = MagicMock()
        fd_gen._pol_drct_FHIaims = MagicMock()
        fd_gen._pol_drct_DFTB = MagicMock()

        fd_gen.iterate_over_files()

        fd_gen._displacements.assert_called_once()
        fd_gen._generate_input_prms.assert_called_once()
        fd_gen._copy_files_in_dir.assert_called_once_with('control.in')

        fd_gen._pol_drct_FHIaims.assert_not_called()
    
    @patch('os.listdir', return_value=['Coord-0-0-x-geometry.in-001_+', 'Coord-0-0-x-geometry.in-002_-'])
    @patch('os.path.isdir', return_value=True)
    @patch('os.path.join', side_effect=lambda *args: '/'.join(args))
    @patch('shutil.copy')
    @patch('os.mkdir')
    @patch('THeSeuSS.InputOutputFiles.InputsGenerator')
    def test_iterate_over_files_dftb_plus(self, mock_generator, mock_mkdir, mock_copy, mock_path_join, mock_isdir, mock_listdir):
        
        generator_instance = MagicMock()
        mock_generator.return_value = generator_instance

        fd_gen = finitedisps.FDSubdirectoriesGeneration(code='dftb+', kpoints='1 1 1', functional=None, eev=None, rho=None, etot=None, forces=None, sc_iter_limit=None, species=None, pol_grid=None, SCC_tolerance='1E-7', max_SCC_iterations='100', output_file='output.txt', dispersion=True, dispersion_type='TS')

        fd_gen._displacements = MagicMock()
        fd_gen._generate_input_prms = MagicMock()
        fd_gen._copy_files_in_dir = MagicMock()
        fd_gen._pol_drct_FHIaims = MagicMock()
        fd_gen._pol_drct_DFTB = MagicMock()

        fd_gen.iterate_over_files()

        fd_gen._displacements.assert_called_once()
        fd_gen._generate_input_prms.assert_called_once()
        fd_gen._copy_files_in_dir.assert_called_once_with('dftb_in.hsd')
        
        fd_gen._pol_drct_DFTB.assert_called_once()

        fd_gen._pol_drct_FHIaims.assert_not_called()
