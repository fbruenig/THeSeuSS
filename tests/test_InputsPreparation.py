#Author: Ariadni Boziki

import unittest
from unittest.mock import patch, mock_open, MagicMock
from THeSeuSS import InputsPreparation as inputs 

class Test_GeometryProcessor(unittest.TestCase):

    def setUp(self):
        
        self.aims_geom_input = """\
        lattice_vector 1.0 0.0 0.0
        lattice_vector 0.0 1.0 0.0
        lattice_vector 0.0 0.0 1.0
        atom 0.0 0.0 0.0 H
        atom 1.0 0.0 0.0 H\
        """
        
        self.dftb_geom_input = """\
        3 S
        H N
        1 1 0.0 0.0 0.0
        2 1 1.0 0.0 0.0
        3 2 0.0 1.0 0.0
        0.0 0.0 0.0
        1.0 0.0 0.0
        0.0 1.0 0.0
        0.0 0.0 1.0\
        """

    @patch("builtins.open", new_callable=MagicMock)
    def test_read_lattice_aims(self, mock_open):

        mock_file = MagicMock()
        mock_file.__iter__.return_value = self.aims_geom_input.splitlines()
        mock_open.return_value.__enter__.return_value = mock_file
        geom_processor = inputs.GeometryProcessor("dummy_path", "aims")
        expected_lattice = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        self.assertEqual(geom_processor.read_lattice(), expected_lattice)

    @patch("builtins.open", new_callable=MagicMock)
    def test_read_lattice_dftb(self, mock_open):

        mock_file = MagicMock()
        mock_file.readlines.return_value = self.dftb_geom_input.splitlines()
        mock_open.return_value.__enter__.return_value = mock_file
        geom_processor = inputs.GeometryProcessor("dummy_path", "dftb+")
        expected_lattice = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        self.assertEqual(geom_processor.read_lattice(), expected_lattice)

    @patch("builtins.open", new_callable=MagicMock)
    def test_read_coordinates_aims(self, mock_open):
       
        mock_file = MagicMock()
        mock_file.__iter__.return_value = self.aims_geom_input.splitlines()
        mock_open.return_value.__enter__.return_value = mock_file
        geom_processor = inputs.GeometryProcessor("dummy_path", "aims")
        expected_coords = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]
        self.assertEqual(geom_processor.read_coordinates(), expected_coords)

    @patch("builtins.open", new_callable=MagicMock)
    def test_read_coordinates_dftb(self, mock_open):
       
        mock_file = MagicMock()
        mock_file.__iter__.return_value = self.dftb_geom_input.splitlines()
        mock_open.return_value.__enter__.return_value = mock_file
        geom_processor = inputs.GeometryProcessor("dummy_path", "dftb+")
        expected_coords = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]
        self.assertEqual(geom_processor.read_coordinates(), expected_coords)

    @patch("builtins.open", new_callable=MagicMock)
    def test_number_of_atoms_aims(self, mock_open):
       
        mock_file = MagicMock()
        mock_file.__iter__.return_value = self.aims_geom_input.splitlines()
        mock_open.return_value.__enter__.return_value = mock_file
        geom_processor = inputs.GeometryProcessor("dummy_path", "aims")
        self.assertEqual(geom_processor.number_of_atoms(), 2)

    @patch("builtins.open", new_callable=MagicMock)
    def test_number_of_atoms_dftb(self, mock_open):
        
        mock_file = MagicMock()
        mock_file.__iter__.return_value = self.dftb_geom_input.splitlines()
        mock_open.return_value.__enter__.return_value = mock_file
        geom_processor = inputs.GeometryProcessor("dummy_path", "dftb+")
        self.assertEqual(geom_processor.number_of_atoms(), 3)

    @patch("builtins.open", new_callable=MagicMock)
    def test_read_the_atom_type_aims(self, mock_open):
       
        mock_file = MagicMock()
        mock_file.readlines.return_value = self.aims_geom_input.splitlines()
        mock_open.return_value.__enter__.return_value = mock_file
        geom_processor = inputs.GeometryProcessor("dummy_path", "aims")
        expected_atom_types = ([1, 1], ['H', 'H'])
        self.assertEqual(geom_processor.read_the_atom_type(), expected_atom_types)

    @patch("builtins.open", new_callable=MagicMock)
    def test_read_the_atom_type_dftb(self, mock_open):
       
        mock_file = MagicMock()
        mock_file.readlines.return_value = self.dftb_geom_input.splitlines()
        mock_open.return_value.__enter__.return_value = mock_file
        geom_processor = inputs.GeometryProcessor("dummy_path", "dftb+")
        expected_atom_types = (['1', '1', '2'], ['H', 'H', 'N'])
        self.assertEqual(geom_processor.read_the_atom_type(), expected_atom_types)
