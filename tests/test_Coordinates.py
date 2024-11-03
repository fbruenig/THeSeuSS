#AUTHOR: Ariadni Boziki


import unittest
from unittest.mock import patch, mock_open, MagicMock
import numpy as np
import os
from THeSeuSS import Coordinates as coords


class Test_ConversionUnitsCoordinates(unittest.TestCase):

    @patch('os.getcwd', return_value='/path')
    @patch('os.path.isfile', return_value=False)
    @patch('shutil.copy')
    @patch('os.rename')
    @patch('builtins.open', new_callable=mock_open, read_data="3 S\nH C O\n0.1 0.2 0.3\n0.4 0.5 0.6\n0.7 0.8 0.9\n")
    def test_cartesian_to_fractional_DFTBplus(self, mock_open, mock_rename, mock_copy, mock_isfile, mock_getcwd):

        with patch('THeSeuSS.InputsPreparation.GeometryProcessor') as MockGeometryProcessor:
            mock_geometry_processor = MockGeometryProcessor.return_value
            mock_geometry_processor.number_of_atoms.return_value = 3
            mock_geometry_processor.read_lattice.return_value = np.eye(3).tolist()
            mock_geometry_processor.read_coordinates.return_value = [[0.1, 0.2, 0.3], [0.4, 0.5, 0.6], [0.7, 0.8, 0.9]]
            mock_geometry_processor.read_the_atom_type.return_value = ([1, 2, 3], ["H", "C", "O"])

            conversion = coords.ConversionUnitsCoordinates(code="dftb+")
            path = '/path'
            geo_path = conversion.cartesian_to_fractional_DFTBplus(path, supercell=True)

            self.assertEqual(geo_path, os.path.join(path, 'geo.genS'))
            mock_copy.assert_called_with(os.path.join(path, 'geo.gen'), os.path.join(path, 'geo_backup.gen'))

            mock_geometry_processor.number_of_atoms.assert_called_once()
            mock_geometry_processor.read_lattice.assert_called_once()
            mock_geometry_processor.read_coordinates.assert_called_once()
            mock_geometry_processor.read_the_atom_type.assert_called_once()
            mock_rename.assert_called_once_with(os.path.join(path, "geo.gen.F"), os.path.join(path, "geo.gen"))

    @patch('os.getcwd', return_value='/path')
    @patch('os.path.isfile', return_value=False)
    @patch('shutil.copy')
    @patch('os.rename')
    @patch('builtins.open', new_callable=mock_open, read_data="lattice_vector 1.0 0.0 0.0\nlattice_vector 0.0 1.0 0.0\n")
    def test_cartesian_to_fractional_FHIaims(self, mock_open, mock_rename, mock_copy, mock_isfile, mock_getcwd):

        with patch('THeSeuSS.InputsPreparation.GeometryProcessor') as MockGeometryProcessor:
            mock_geometry_processor = MockGeometryProcessor.return_value
            mock_geometry_processor.number_of_atoms.return_value = 3
            mock_geometry_processor.read_lattice.return_value = np.eye(3).tolist()
            mock_geometry_processor.read_coordinates.return_value = [[0.1, 0.2, 0.3], [0.4, 0.5, 0.6], [0.7, 0.8, 0.9]]
            mock_geometry_processor.read_the_atom_type.return_value = ([1, 2, 3], ["H", "C", "O"])

            conversion = coords.ConversionUnitsCoordinates(code="aims")
            path = '/path'
            geo_path = conversion.cartesian_to_fractional_FHIaims(path, supercell=True)

            self.assertEqual(geo_path, os.path.join(path, 'geometry.in.supercell'))
            mock_copy.assert_called_with(os.path.join(path, 'geometry.in'), os.path.join(path, 'geometry_backup.in'))

            mock_geometry_processor.number_of_atoms.assert_called_once()
            mock_geometry_processor.read_lattice.assert_called_once()
            mock_geometry_processor.read_coordinates.assert_called_once()
            mock_geometry_processor.read_the_atom_type.assert_called_once()
            mock_rename.assert_called_once_with(os.path.join(path, "geometry.in.F"), os.path.join(path, "geometry.in"))
