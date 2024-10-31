#AUTHOR: Ariadni Boziki


import unittest
import numpy as np
from unittest.mock import patch, MagicMock
from THeSeuSS import EigenvectorsFrequenciesPHONOPY as eigenvec


class Test_VibrationalFrequencies(unittest.TestCase):

    @patch.object(eigenvec.VibrationalFrequencies, '_get_number_of_atoms', return_value=2)
    @patch.object(eigenvec.VibrationalFrequencies, '_get_sorted_directories', return_value=['Coord-0-0-x-geometry.in-001_+', 'Coord-0-0-x-geometry.in-002_-'])
    @patch.object(eigenvec.VibrationalFrequencies, '_forces_output_path', side_effect=lambda drct: f'/path/to/force_{drct}.out')
    @patch.object(eigenvec.VibrationalFrequencies, '_get_forces', side_effect=lambda path: [0.1, -0.2, 0.3] if '001' in path else [-0.1, 0.2, -0.3])
    def test_set_forces_dictionary(self, mock_get_number_of_atoms, mock_get_sorted_directories, mock_forces_output_path, mock_get_forces):
        
        vib_freq = eigenvec.VibrationalFrequencies(code='aims', output_file_of_SCFSs='output', dispersion=True)
        vib_freq.displacement = 0.01
        vib_freq._set_forces_dictionary()

        expected_dataset = {
            'natom': 2,
            'atoms': [
                {
                    'atom_number': '0',
                    'atom_axis': 'x',
                    'disp': 0.01,
                    'forces': [0.1, -0.2, 0.3]
                },
                {
                    'atom_number': '0',
                    'atom_axis': 'x',
                    'disp': -0.01,
                    'forces': [-0.1, 0.2, -0.3]
                }
            ]
        }

        self.assertEqual(vib_freq.dataset, expected_dataset)

    @patch.object(eigenvec.VibrationalFrequencies, '_set_forces_dictionary')
    def test_group_dispplacement_forces_dict(self, mock_set_forces_dictionary):
        
        vib_freq = eigenvec.VibrationalFrequencies(code='aims', output_file_of_SCFSs='output', dispersion=True)
        vib_freq.no_of_atoms = 2

        vib_freq.dataset = {
            'natom': 2,
            'atoms': [
                {'atom_number': '1', 'atom_axis': 'x', 'disp': 0.01, 'forces': [0.1, 0.2, 0.3]},
                {'atom_number': '1', 'atom_axis': 'x', 'disp': -0.01, 'forces': [-0.1, -0.2, -0.3]},
                {'atom_number': '2', 'atom_axis': 'y', 'disp': 0.01, 'forces': [0.4, 0.5, 0.6]},
                {'atom_number': '2', 'atom_axis': 'y', 'disp': -0.01, 'forces': [-0.4, -0.5, -0.6]}
            ]
        }

        vib_freq._group_dispplacement_forces_dict()

        expected_delta_forces = {
            'natom': 2,
            'atoms': [
                {
                    'number': '1',
                    'axis': 'x',
                    'difference': np.array([0.2, 0.4, 0.6])
                },
                {
                    'number': '2',
                    'axis': 'y',
                    'difference': np.array([0.8, 1.0, 1.2])
                }
            ]
        }

        self.assertEqual(vib_freq.delta_forces['natom'], expected_delta_forces['natom'])
        for actual, expected in zip(vib_freq.delta_forces['atoms'], expected_delta_forces['atoms']):
            self.assertEqual(actual['number'], expected['number'])
            self.assertEqual(actual['axis'], expected['axis'])
            np.testing.assert_array_almost_equal(actual['difference'], expected['difference'])

        mock_set_forces_dictionary.assert_called_once()

    @patch('THeSeuSS.EigenvectorsFrequenciesPHONOPY.element')
    @patch.object(eigenvec.VibrationalFrequencies, '_get_atom_type', return_value=(['H', 'O'], ['H', 'O']))
    def test_mass_matrix(self, mock_get_atom_type, mock_element):
        
        mock_element.side_effect = lambda el: MagicMock(atomic_weight=1.008 if el == 'H' else 15.999)
        vib_freq = eigenvec.VibrationalFrequencies(code='aims', output_file_of_SCFSs='output', dispersion=True)
        vib_freq.no_of_atoms = 2

        d_forces = {
            'atoms': [
                {'number': '0'},
                {'number': '1'}
            ]
        }

        vib_freq._mass_matrix(d_forces)

        mass_h = 1.008
        mass_o = 15.999
        expected_mass_matrix = np.zeros((6, 6))

        for i in range(3):
            for j in range(3):
                expected_mass_matrix[i, j] = 1.0 / np.sqrt(mass_h * mass_h)
                expected_mass_matrix[i, j+3] = 1.0 / np.sqrt(mass_h * mass_o)
                expected_mass_matrix[i+3, j] = 1.0 / np.sqrt(mass_o * mass_h)
                expected_mass_matrix[i+3, j+3] = 1.0 / np.sqrt(mass_o * mass_o)

        np.testing.assert_array_almost_equal(vib_freq.mass, expected_mass_matrix, decimal=6)

        mock_get_atom_type.assert_called_once()
        mock_element.assert_any_call('H')
        mock_element.assert_any_call('O')
        self.assertEqual(mock_element.call_count, 2)

    def test_generation_of_fc_matrix(self):
        
        vib_freq = eigenvec.VibrationalFrequencies(code='aims', output_file_of_SCFSs='output', dispersion=True)
        vib_freq.no_of_atoms = 3

        dforces = {
            'atoms': [
                {'difference': np.array([[-1.24588359e+00, 1.18349774e-12, -6.08472842e-29],
                    [6.22941796e-01, -4.84793723e-01, -6.08472842e-29],
                    [6.22941796e-01, 4.84793723e-01, -6.08472842e-29]])},
                {'difference': np.array([[1.49802916e-13, -8.31995937e-01, -4.05648562e-29],
                    [-3.70617094e-01, 4.15997969e-01, 1.01412140e-29],
                    [3.70617094e-01, 4.15997969e-01, 3.04236421e-29]])},
                {'difference': np.array([[-1.04851416e-13, -6.72050089e-13, -1.17756608e-04],
                    [5.07499858e-14, 2.95249936e-13, 5.88783039e-05],
                    [5.41099283e-14, 3.76809911e-13, 5.88783039e-05]])},
                {'difference': np.array([[6.22948104e-01, -3.70742394e-01, -3.04236421e-29],
                    [-6.81210014e-01, 4.27767801e-01, 1.41976997e-28],
                    [5.82619105e-02, -5.70254075e-02, 1.41976997e-28]])},
                {'difference': np.array([[-4.84599324e-01, 4.16005097e-01, -1.01412140e-29],
                    [4.27601409e-01, -3.96421136e-01, 1.01412140e-29],
                    [5.69979149e-02, -1.95839614e-02, 2.02824281e-29]])},
                {'difference': np.array([[ 2.37710025e-13, 2.79094998e-12, 5.15945421e-05],
                    [1.55610334e-13, -1.30301000e-12, -5.21544877e-05],
                    [-3.93322500e-13, -1.48793130e-12, 5.59945616e-07]])},
                {'difference': np.array([[ 6.22948104e-01,  3.70742394e-01, -1.62259425e-28],
                    [5.82619105e-02, 5.70254075e-02, 1.21694568e-28],
                    [-6.81210014e-01, -4.27767801e-01, 0.00000000e+00]])},
                {'difference': np.array([[4.84599324e-01, 4.16005097e-01, -1.01412140e-29],
                    [-5.69979149e-02, -1.95839614e-02, 3.04236421e-29],
                    [-4.27601409e-01, -3.96421136e-01, 1.01412140e-29]])},
                {'difference': np.array([[ 9.36100156e-14,  2.59919907e-13,  5.15945421e-05],
                    [-8.18950993e-14, -9.35504996e-14,  5.59945609e-07],
                    [-1.17098171e-14, -1.66369956e-13, -5.21544877e-05]])}
                ]}

        expected_force_constants = np.array([
            [1.24588359e+00, -1.18349774e-12, 6.08472842e-29, -6.22941796e-01, 4.84793723e-01, 6.08472842e-29, -6.22941796e-01, -4.84793723e-01, 6.08472842e-29],
            [-1.49802916e-13, 8.31995937e-01, 4.05648562e-29, 3.70617094e-01, -4.15997969e-01, -1.01412140e-29, -3.70617094e-01, -4.15997969e-01, -3.04236421e-29],
            [1.04851416e-13, 6.72050089e-13, 1.17756608e-04, -5.07499858e-14, -2.95249936e-13, -5.88783039e-05, -5.41099283e-14, -3.76809911e-13, -5.88783039e-05], 
            [-6.22948104e-01, 3.70742394e-01, 3.04236421e-29, 6.81210014e-01, -4.27767801e-01, -1.41976997e-28, -5.82619105e-02, 5.70254075e-02, -1.41976997e-28],
            [4.84599324e-01, -4.16005097e-01, 1.01412140e-29, -4.27601409e-01, 3.96421136e-01, -1.01412140e-29, -5.69979149e-02, 1.95839614e-02, -2.02824281e-29],
            [-2.37710025e-13, -2.79094998e-12, -5.15945421e-05, -1.55610334e-13, 1.30301000e-12, 5.21544877e-05, 3.93322500e-13, 1.48793130e-12, -5.59945616e-07],
            [-6.22948104e-01, -3.70742394e-01, 1.62259425e-28, -5.82619105e-02, -5.70254075e-02, -1.21694568e-28, 6.81210014e-01, 4.27767801e-01, -0.00000000e+00],
            [-4.84599324e-01, -4.16005097e-01, 1.01412140e-29, 5.69979149e-02, 1.95839614e-02, -3.04236421e-29, 4.27601409e-01, 3.96421136e-01, -1.01412140e-29],
            [-9.36100156e-14, -2.59919907e-13, -5.15945421e-05, 8.18950993e-14, 9.35504996e-14, -5.59945609e-07, 1.17098171e-14, 1.66369956e-13, 5.21544877e-05]
            ])

        fc_matrix = vib_freq._generation_of_fc_matrix(dforces)
        
        np.testing.assert_array_almost_equal(fc_matrix, expected_force_constants, decimal=6)

    def test_symmetrization_of_fc_matrix(self):
        
        vib_freq = eigenvec.VibrationalFrequencies(code='aims', output_file_of_SCFSs='output', dispersion=True)
        vib_freq.no_of_atoms = 3

        force_constants = np.array([[6.22941796e+01, -5.91748872e-11, 3.04236421e-27, -3.11470898e+01, 2.42396861e+01, 3.04236421e-27, -3.11470898e+01, -2.42396861e+01, 3.04236421e-27],
            [-7.49014581e-12, 4.15997969e+01, 2.02824281e-27, 1.85308547e+01, -2.07998984e+01, -5.07060702e-28, -1.85308547e+01, -2.07998984e+01, -1.52118211e-27],
            [5.24257079e-12, 3.36025045e-11, 5.88783039e-03, -2.53749929e-12, -1.47624968e-11, -2.94391520e-03, -2.70549642e-12, -1.88404956e-11, -2.94391520e-03],
            [-3.11474052e+01, 1.85371197e+01, 1.52118211e-27, 3.40605007e+01, -2.13883901e+01, -7.09884983e-27, -2.91309553e+00, 2.85127038e+00, -7.09884983e-27],
            [2.42299662e+01, -2.08002548e+01, 5.07060702e-28, -2.13800704e+01, 1.98210568e+01, -5.07060702e-28, -2.84989574e+00, 9.79198070e-01, -1.01412140e-27],
            [-1.18855013e-11, -1.39547499e-10, -2.57972710e-03, -7.78051668e-12, 6.51505000e-11, 2.60772438e-03, 1.96661250e-11, 7.43965650e-11, -2.79972808e-05],
            [-3.11474052e+01, -1.85371197e+01, 8.11297123e-27, -2.91309553e+00, -2.85127038e+00, -6.08472842e-27, 3.40605007e+01, 2.13883901e+01, -0.00000000e+00],
            [-2.42299662e+01, -2.08002548e+01, 5.07060702e-28, 2.84989574e+00, 9.79198070e-01, -1.52118211e-27, 2.13800704e+01, 1.98210568e+01, -5.07060702e-28],
            [-4.68050078e-12, -1.29959953e-11, -2.57972710e-03, 4.09475496e-12, 4.67752498e-12, -2.79972804e-05, 5.85490857e-13, 8.31849780e-12, 2.60772438e-03]])

        expected_symmetric_fc = np.array([[6.22941796e+01, -3.33325165e-11, 2.62128540e-12, -3.11472475e+01, 2.42348262e+01, -5.94275063e-12, -3.11472475e+01, -2.42348262e+01, -2.34025039e-12],
            [-3.33325165e-11, 4.15997969e+01, 1.68012522e-11, 1.85339872e+01, -2.08000766e+01, -6.97737495e-11, -1.85339872e+01, -2.08000766e+01, -6.49799767e-12],
            [2.62128540e-12, 1.68012522e-11, 5.88783039e-03, -1.26874965e-12, -7.38124839e-12, -2.76182115e-03, -1.35274821e-12, -9.42024778e-12, -2.76182115e-03],
            [-3.11472475e+01, 1.85339872e+01, -1.26874965e-12, 3.40605007e+01, -2.13842303e+01, -3.89025834e-12, -2.91309553e+00, 2.85058306e+00, 2.04737748e-12],
            [2.42348262e+01, -2.08000766e+01, -7.38124839e-12, -2.13842303e+01, 1.98210568e+01, 3.25752500e-11, -2.85058306e+00, 9.79198070e-01, 2.33876249e-12],
            [-5.94275063e-12, -6.97737495e-11, -2.76182115e-03, -3.89025834e-12, 3.25752500e-11, 2.60772438e-03, 9.83306251e-12, 3.71982825e-11, -2.79972806e-05],
            [-3.11472475e+01, -1.85339872e+01, -1.35274821e-12, -2.91309553e+00, -2.85058306e+00, 9.83306251e-12, 3.40605007e+01, 2.13842303e+01, 2.92745429e-13],
            [-2.42348262e+01, -2.08000766e+01, -9.42024778e-12, 2.85058306e+00, 9.79198070e-01, 3.71982825e-11, 2.13842303e+01, 1.98210568e+01, 4.15924890e-12],
            [-2.34025039e-12, -6.49799767e-12, -2.76182115e-03, 2.04737748e-12, 2.33876249e-12, -2.79972806e-05, 2.92745429e-13, 4.15924890e-12, 2.60772438e-03]])

        symmetrized_fc_matrix = vib_freq._symmetrization_of_fc_matrix(force_constants)

        np.testing.assert_array_almost_equal(symmetrized_fc_matrix, expected_symmetric_fc, decimal=6)

    def test_read_eig_vec_phonopy(self):
        
        vib_freq = eigenvec.VibrationalFrequencies(code='aims', output_file_of_SCFSs='output', dispersion=True)

        hessian = np.array([[3.89362957e+00, -2.08341250e-12, 1.63840577e-13, -7.75609265e+00, 6.03480475e+00, -1.47982657e-12, -7.75609265e+00, -6.03480475e+00, -5.82754507e-13],
            [-2.08341250e-12, 2.60014981e+00, 1.05014390e-12, 4.61521750e+00, -5.17950492e+00, -1.73746225e-11, -4.61521750e+00, -5.17950492e+00, -1.61809072e-12],
            [1.63840577e-13, 1.05014390e-12, 3.68012400e-04, -3.15936097e-13, -1.83803228e-12, -6.87731420e-04, -3.36852894e-13, -2.34577115e-12, -6.87731419e-04],
            [-7.75609265e+00, 4.61521750e+00, -3.15936097e-13, 3.37901793e+01, -2.12145141e+01, -3.85938327e-12, -2.88997572e+00, 2.82795938e+00, 2.03112845e-12],
            [6.03480475e+00, -5.17950492e+00, -1.83803228e-12, -2.12145141e+01, 1.96637468e+01, 3.23167163e-11, -2.82795938e+00, 9.71426657e-01, 2.32020088e-12],
            [-1.47982657e-12, -1.73746225e-11, -6.87731420e-04, -3.85938327e-12, 3.23167163e-11, 2.58702816e-03, 9.75502233e-12, 3.69030580e-11, -2.77750800e-05],
            [-7.75609265e+00, -4.61521750e+00, -3.36852894e-13, -2.88997572e+00, -2.82795938e+00, 9.75502233e-12, 3.37901793e+01, 2.12145141e+01, 2.90422052e-13],
            [-6.03480475e+00, -5.17950492e+00, -2.34577115e-12, 2.82795938e+00, 9.71426657e-01, 3.69030580e-11, 2.12145141e+01, 1.96637468e+01, 4.12623899e-12],
            [-5.82754507e-13, -1.61809072e-12, -6.87731419e-04, 2.03112845e-12, 2.32020088e-12, -2.77750800e-05, 2.90422052e-13, 4.12623899e-12, 2.58702816e-03]])

        expected_eigvals = np.array([-9.36170907e-03, -1.40409904e-06, -1.51945143e-09, 1.40489657e-07, 2.61480324e-03, 2.92866958e-03, 9.32144257e+00, 5.05940357e+01, 5.34955148e+01])

        eigvecs, eigvals, frequencies_in_cm_minus_1 = vib_freq.read_eig_vec_phonopy(hessian)

        np.testing.assert_array_almost_equal(eigvals, expected_eigvals, decimal=6)
