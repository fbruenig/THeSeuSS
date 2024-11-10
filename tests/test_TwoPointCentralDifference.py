#AUTHOR: Ariadni Boziki


import unittest
from unittest.mock import patch, MagicMock, mock_open
import os
import numpy as np
from THeSeuSS import TwoPointCentralDifference as cendiff


class Test_TwoPointCentralDiff(unittest.TestCase):

    def setUp(self):

        with open('geo.gen', 'w') as f:
            f.write('3 S\nH N\n')

    def tearDown(self):

        os.remove('geo.gen')

    @patch('os.listdir')
    @patch('os.path.isdir')
    @patch('builtins.open', new_callable=mock_open, read_data="| Unit cell volume 100.0 ")
    @patch.object(cendiff.TwoPointCentralDiff, '_unit_cell_volume')
    def test_get_volume(self, mock_unit_cell_volume, mock_open, mock_isdir, mock_listdir):

        mock_listdir.return_value = ['Coord_1', 'Coord_2', 'OtherFolder']
        mock_isdir.side_effect = lambda path: True
        mock_unit_cell_volume.side_effect = lambda path: setattr(self.two_point_diff, 'volume', np.float64(100.0))
        self.two_point_diff = cendiff.TwoPointCentralDiff(code='aims', output_file='output_file.txt', dispersion=False, supercell=False, restart=False)
        
        result_volume = self.two_point_diff.get_volume()

        self.assertEqual(result_volume, np.float64(100.0))

    @patch('builtins.open', new_callable=mock_open, read_data="Polarizability 1.0 2.0 3.0 4.0 5.0 6.0\nTotal dipole moment 1.0 2.0 3.0\n")
    def test_find_pattern_aims_non_periodic(self, mock_open):

        two_point_diff = cendiff.TwoPointCentralDiff(code='aims', output_file='output_file.txt', dispersion=False, supercell=False, restart=False, functional='pbe')
        two_point_diff.non_periodic = True

        drcts = "mock_directory"
        pol, cartesian_pol = two_point_diff._find_pattern(drcts)

        np.testing.assert_array_equal(pol, np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0]))
        np.testing.assert_array_equal(cartesian_pol, np.array([1.0, 2.0, 3.0]))

    @patch('builtins.open', new_callable=mock_open, read_data="Dipole moment:   -0.00000000   -1.88930279   -0.00000000 Debye\n")
    @patch('os.path.join', side_effect=lambda *args: '/'.join(args))
    def test_find_pattern_dftb_plus(self, mock_path_join, mock_open):
        
        two_point_diff = cendiff.TwoPointCentralDiff(code='dftb+', output_file='output_file.txt', dispersion=False, supercell=False, restart=False)
        two_point_diff.non_periodic = False

        polarizability, cartesian_polarization = two_point_diff._find_pattern('Coord_2')

        expected_cartesian_polarization = np.array([-0.00000000, -1.88930279, -0.00000000])
        self.assertTrue(np.array_equal(cartesian_polarization, expected_cartesian_polarization))

    @patch('builtins.open', new_callable=mock_open, read_data="Static polarisability:\n1.0 0.0 0.0\n0.0 2.0 0.0\n0.0 0.0 3.0\n")
    @patch('os.path.join', side_effect=lambda *args: '/'.join(args))
    def test_find_pattern_dftb_plus_polarizability(self, mock_path_join, mock_open):
        
        two_point_diff = cendiff.TwoPointCentralDiff(code='dftb+', output_file='', dispersion=True, supercell=False, restart=False)
        two_point_diff.non_periodic = True

        polarizability, cartesian_polarization = two_point_diff._find_pattern('Coord_3')

        expected_polarizability = [1.0, 2.0, 3.0, 0.0, 0.0, 0.0]
        self.assertTrue(np.array_equal(polarizability, expected_polarizability))

    @patch('os.listdir', return_value=[
        'Coord-0-0-x-geometry.in-001_+',
        'Coord-1-0-y-geometry.in-004_-',
        'Coord-0-0-x-geometry.in-002_-',
        'Coord-1-0-y-geometry.in-003_+'
    ])
    @patch('os.path.isdir', return_value=True)
    def test_subgroups_finite_disp_drct(self, mock_isdir, mock_listdir):
        
        two_point_diff = cendiff.TwoPointCentralDiff(code='dftb+', output_file='output_file.txt', dispersion=True, supercell=False, restart=False)

        two_point_diff.path = '/current_working_directory'
        result = two_point_diff._subgroups_finite_disp_drct()

        expected_result = [
            ('Coord-0-0-x-geometry.in-001_+', 'Coord-0-0-x-geometry.in-002_-'),
            ('Coord-1-0-y-geometry.in-003_+', 'Coord-1-0-y-geometry.in-004_-')
        ]

        self.assertEqual(result, expected_result)
        mock_isdir.assert_called()

    @patch.object(cendiff.TwoPointCentralDiff, '_assign_axis_to_numbers', side_effect=[1, -1])
    @patch.object(cendiff.TwoPointCentralDiff, '_find_pattern', side_effect=[
        (np.array([9.73573484, 9.17625266, 8.56369404, 0.0837749, -0., 0.]), np.array([-4.59836487e-03, -3.85670744e-01, 2.72932663e-14])),
        (np.array([9.73573484, 9.17625266, 8.56369404, -0.0837749, -0., 0.]), np.array([4.59836487e-03, -3.85670744e-01, -1.85431178e-16]))
    ])
    @patch.object(cendiff.TwoPointCentralDiff, '_map_axis_coord')
    @patch.object(cendiff.TwoPointCentralDiff, '_subgroups_finite_disp_drct')
    def test_pol(self, mock_subgroups, mock_map_axis, mock_find_pattern, mock_assign_axis):
        
        two_point_diff = cendiff.TwoPointCentralDiff(code='dftb+', output_file='output_file.txt', dispersion=True, supercell=False, restart=False)

        two_point_diff.sign_atom_coord = [
            ('Coord-0-0-x-geometry.in-001_+', 'Coord-0-0-x-geometry.in-002_-')
        ]
        two_point_diff.pol = np.empty([0, 6])
        two_point_diff.cartesian_pol = np.empty([0, 3])
        two_point_diff.element_axis_coord = np.empty([0, 2])

        two_point_diff._pol(disp=0.01)

        expected_pol = np.array([[-5.00000397e-08, 5.00000397e-08, 0.00000000e+00, -8.37749020e+00, 0.00000000e+00, 0.00000000e+00]])
        expected_cartesian_pol = np.array([[4.59836487e-01, -9.98312544e-13, -1.37393487e-12]])

        np.testing.assert_array_almost_equal(two_point_diff.pol, expected_pol)
        np.testing.assert_array_almost_equal(two_point_diff.cartesian_pol, expected_cartesian_pol)

        mock_subgroups.assert_called_once()
        mock_find_pattern.assert_any_call('Coord-0-0-x-geometry.in-001_+')
        mock_find_pattern.assert_any_call('Coord-0-0-x-geometry.in-002_-')
        mock_map_axis.assert_called_once_with('x')
        self.assertEqual(mock_assign_axis.call_count, 2)
