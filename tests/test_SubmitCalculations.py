#AUTHOR: Ariadni Boziki


import unittest
from unittest.mock import patch, MagicMock, call
import os
from THeSeuSS import SubmitCalculations as submit


class Test_PhonopyCalculator(unittest.TestCase):

    @patch('os.listdir')
    @patch('os.path.isdir')
    def test_sort_directories(self, mock_isdir, mock_listdir):
        
        directory_names = [
            'Coord-0-0-x-geometry.in-001_+',
            'Coord-5-1-z-geometry.in-012_-',
            'Coord-0-0-x-geometry.in-002_-',
            'Coord-1-0-y-geometry.in-004_-',
            'Coord-3-1-x-geometry.in-007_+',
            'Coord-1-0-y-geometry.in-003_+',
            'Coord-3-1-x-geometry.in-008_-',
            'Coord-4-1-y-geometry.in-010_-',
            'Coord-4-1-y-geometry.in-009_+',
            'Coord-5-1-z-geometry.in-011_+',
        ]

        expected_sorted_directories = [
            'Coord-0-0-x-geometry.in-001_+',
            'Coord-0-0-x-geometry.in-002_-',
            'Coord-1-0-y-geometry.in-003_+',
            'Coord-1-0-y-geometry.in-004_-',
            'Coord-3-1-x-geometry.in-007_+',
            'Coord-3-1-x-geometry.in-008_-',
            'Coord-4-1-y-geometry.in-009_+',
            'Coord-4-1-y-geometry.in-010_-',
            'Coord-5-1-z-geometry.in-011_+',
            'Coord-5-1-z-geometry.in-012_-'
        ]

        mock_listdir.return_value = directory_names
        mock_isdir.return_value = True

        calc = submit.PhonopyCalculator(code='aims', cell_dims='2 2 2', output_file_of_SCFSs='output')
        sorted_directories = calc.sort_directories()
        self.assertEqual(sorted_directories, expected_sorted_directories)

class Test_Calculator(unittest.TestCase):

    @patch('subprocess.run')
    @patch('os.environ.get')
    def test_submit_job_with_number_of_cores(self, mock_get_env, mock_subprocess_run):
        
        mock_get_env.return_value = '4'
        commands = "echo 'Starting calculation'; run_calculation"

        calc = submit.Calculator(code='aims', output_file='output', dispersion=True, restart=False, commands=commands)
        calc.submit_job()

        expected_command = "echo 'Starting calculation'; srun --cpus-per-task 1 --ntasks 4 run_calculation"

        mock_subprocess_run.assert_called_once_with(expected_command, shell=True, executable='/bin/bash')

    @patch('subprocess.run')
    @patch('os.environ.get')
    def test_submit_job_without_number_of_cores(self, mock_get_env, mock_subprocess_run):
        
        mock_get_env.return_value = None
        commands = "run_calculation"

        with patch('multiprocessing.cpu_count', return_value=8):
            calc = submit.Calculator(code='aims', output_file='output', dispersion=True, restart=False, commands=commands)
            calc.submit_job()

            expected_command = "run_calculation"

            mock_subprocess_run.assert_called_once_with(expected_command, shell=True, executable='/bin/bash')

    @patch('subprocess.run')
    @patch('os.environ.get')
    def test_submit_job_with_single_command_and_cores(self, mock_get_env, mock_subprocess_run):
        
        mock_get_env.return_value = '4'
        commands = "run_calculation"

        calc = submit.Calculator(code='aims', output_file='output', dispersion=True, restart=False, commands=commands)
        calc.submit_job()

        expected_command = "srun --cpus-per-task 1 --ntasks 4 run_calculation"

        mock_subprocess_run.assert_called_once_with(expected_command, shell=True, executable='/bin/bash')

    @patch('os.listdir')
    @patch('os.environ.get')
    @patch('os.path.isdir')
    @patch('subprocess.run')
    def test_submit_jobs_in_parallel_with_cores_and_dispersion(self, mock_subprocess_run, mock_isdir, mock_get_env, mock_listdir):
        
        mock_get_env.return_value = '8'
        mock_listdir.return_value = ['Coord-0-0-x-geometry.in-001_+', 'Coord-0-0-x-geometry.in-002_-', 'Coord-1-0-y-geometry.in-003_+', 'Coord-1-0-y-geometry.in-004_-']
        mock_isdir.return_value = True
        mock_subprocess_run.return_value = MagicMock(stdout='Command executed')

        calc = submit.Calculator(code='dftb+', output_file='output', dispersion=True, restart=False, commands="echo 'Start job'; run_calculation")

        with patch.object(calc, 'run_command', wraps=calc.run_command) as mock_run_command:
            calc.submit_jobs_in_parallel()

            expected_commands = [
                "cd vibrations; cd Coord-0-0-x-geometry.in-001_+; echo 'Start job'; srun --cpus-per-task 1 --ntasks 2 run_calculation",
                "cd vibrations; cd Coord-0-0-x-geometry.in-001_+; cd polarizability; srun --cpus-per-task 1 --ntasks 2 run_calculation",
                "cd vibrations; cd Coord-0-0-x-geometry.in-002_-; echo 'Start job'; srun --cpus-per-task 1 --ntasks 2 run_calculation",
                "cd vibrations; cd Coord-0-0-x-geometry.in-002_-; cd polarizability; srun --cpus-per-task 1 --ntasks 2 run_calculation",
                "cd vibrations; cd Coord-1-0-y-geometry.in-003_+; echo 'Start job'; srun --cpus-per-task 1 --ntasks 2 run_calculation",
                "cd vibrations; cd Coord-1-0-y-geometry.in-003_+; cd polarizability; srun --cpus-per-task 1 --ntasks 2 run_calculation",
                "cd vibrations; cd Coord-1-0-y-geometry.in-004_-; echo 'Start job'; srun --cpus-per-task 1 --ntasks 2 run_calculation",
                "cd vibrations; cd Coord-1-0-y-geometry.in-004_-; cd polarizability; srun --cpus-per-task 1 --ntasks 2 run_calculation"
            ]

            mock_run_command.assert_has_calls([call(cmd) for cmd in expected_commands], any_order=True)
            self.assertEqual(mock_run_command.call_count, len(expected_commands))

    @patch('os.listdir')
    @patch('os.environ.get')
    @patch('os.path.isdir')
    @patch('subprocess.run')
    def test_submit_jobs_in_parallel_without_number_of_cores(self, mock_subprocess_run, mock_isdir, mock_get_env, mock_listdir):
        
        mock_get_env.return_value = None
        mock_listdir.return_value = ['Coord-0-0-x-geometry.in-001_+', 'Coord-0-0-x-geometry.in-002_-', 'Coord-1-0-y-geometry.in-003_+', 'Coord-1-0-y-geometry.in-004_-']
        mock_isdir.return_value = True
        mock_subprocess_run.return_value = MagicMock(stdout='Command executed')

        with patch('multiprocessing.cpu_count', return_value=4):
            calc = submit.Calculator(code='aims', output_file='output', dispersion=False, restart=False, commands="run_calculation", functional="pbe")

            with patch.object(calc, 'run_command', wraps=calc.run_command) as mock_run_command:
                calc.submit_jobs_in_parallel()


                expected_commands = [
                    "cd vibrations; cd Coord-0-0-x-geometry.in-001_+; run_calculation",
                    "cd vibrations; cd Coord-0-0-x-geometry.in-002_-; run_calculation",
                    "cd vibrations; cd Coord-1-0-y-geometry.in-003_+; run_calculation",
                    "cd vibrations; cd Coord-1-0-y-geometry.in-004_-; run_calculation"
                ]

                mock_run_command.assert_has_calls([call(cmd) for cmd in expected_commands], any_order=True)
                self.assertEqual(mock_run_command.call_count, len(expected_commands))

    @patch('subprocess.run')
    @patch('THeSeuSS.Restart.RestartCalculation')
    @patch('os.environ.get')
    def test_submit_jobs_in_parallel_with_restart(self, mock_get_env, mock_restart, mock_subprocess_run):

        mock_get_env.return_value = '8'

        mock_restart_instance = mock_restart.return_value
        mock_restart_instance.directory_non_completed_calculations.return_value = [
            'Coord-0-0-x-geometry.in-001_+/polarizability',
            'Coord-1-0-y-geometry.in-004_-',
            'Coord-0-0-x-geometry.in-002_-',
            'Coord-1-0-y-geometry.in-003_+/polarizability'
        ]

        mock_subprocess_run.return_value = MagicMock(stdout='Command executed')

        calc = submit.Calculator(code='aims', output_file='output', dispersion=False, restart=True, commands="echo 'Start job'; run_calculation")

        with patch.object(calc, 'run_command', wraps=calc.run_command) as mock_run_command:
            calc.submit_jobs_in_parallel()

            expected_commands = [
                "cd Coord-0-0-x-geometry.in-001_+/polarizability; echo 'Start job'; srun --cpus-per-task 1 --ntasks 2 run_calculation",
                "cd Coord-0-0-x-geometry.in-002_-; echo 'Start job'; srun --cpus-per-task 1 --ntasks 2 run_calculation",
                "cd Coord-1-0-y-geometry.in-003_+/polarizability; echo 'Start job'; srun --cpus-per-task 1 --ntasks 2 run_calculation",
                "cd Coord-1-0-y-geometry.in-004_-; echo 'Start job'; srun --cpus-per-task 1 --ntasks 2 run_calculation"
            ]

            mock_run_command.assert_has_calls([call(cmd) for cmd in expected_commands], any_order=True)
            self.assertEqual(mock_run_command.call_count, len(expected_commands))

            mock_restart.assert_called_once_with('aims', 'output', False, True, None, "echo 'Start job'; run_calculation")
            mock_restart_instance.directory_non_completed_calculations.assert_called_once()
