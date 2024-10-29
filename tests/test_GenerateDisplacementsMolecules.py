#AUTHOR: Ariadni Boziki


import os
import unittest
from unittest.mock import patch, mock_open, MagicMock, call
import numpy as np
from THeSeuSS import GenerateDisplacementsMolecules as genmolecules


class Test_GenerateDisplacements(unittest.TestCase):

    def setUp(self):
        
        self.generator = genmolecules.GenerateDisplacements(code='aims')

    @patch('THeSeuSS.GenerateDisplacementsMolecules.GenerateDisplacements._get_no_of_atoms', return_value=2)
    @patch('THeSeuSS.GenerateDisplacementsMolecules.GenerateDisplacements._get_atom_type', return_value=([1, 1], ['H', 'H']))
    @patch('builtins.open', new_callable=mock_open)
    def test_generate_displaced_input_aims(self, mock_file, mock_get_atom_type, mock_get_no_of_atoms):
        
        self.generator.code = 'aims'
        self.generator.disp_inp = 'geometry.in-001'
        sample_coord = [[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]]
        expected_path = os.path.join(self.generator.path, 'vibrations', 'geometry.in-001')

        self.generator._generate_displaced_input(sample_coord)

        mock_file.assert_called_once_with(expected_path, 'w')
        handle = mock_file()
        written_content = ''.join(call[0][0] for call in handle.write.call_args_list)

        expected_content = (
            "atom 0.0 0.0 0.0 H\n"
            "atom 1.0 1.0 1.0 H\n"
        )
        
        self.assertEqual(written_content, expected_content)

    @patch('THeSeuSS.GenerateDisplacementsMolecules.GenerateDisplacements._read_coordinates')
    @patch('THeSeuSS.GenerateDisplacementsMolecules.GenerateDisplacements._get_no_of_atoms')
    @patch('THeSeuSS.GenerateDisplacementsMolecules.GenerateDisplacements._define_the_displacement')
    @patch('THeSeuSS.GenerateDisplacementsMolecules.GenerateDisplacements._name_displaced_input')
    @patch('THeSeuSS.GenerateDisplacementsMolecules.GenerateDisplacements._generate_displaced_input')
    @patch('THeSeuSS.GenerateDisplacementsMolecules.GenerateDisplacements._create_supercell_file')
    def test_create_the_displaced_structures(self, mock_create_supercell_file, mock_generate_displaced_input,
                                         mock_name_displaced_input, mock_define_displacement,
                                         mock_get_no_of_atoms, mock_read_coordinates):

        initial_coords = [[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]]
        mock_read_coordinates.side_effect = lambda: [list(coord) for coord in initial_coords]
        mock_get_no_of_atoms.return_value = 2
        self.generator.disp = 0.01

        self.generator.create_the_displaced_structures()
        mock_create_supercell_file.assert_called_once()

        expected_name_calls = [call(1), call(2), call(3), call(4), call(5), call(6)]
        mock_name_displaced_input.assert_has_calls(expected_name_calls, any_order=True)

        expected_generate_calls = [
            call([[0.01, 0.0, 0.0], [1.0, 1.0, 1.0]]),
            call([[-0.01, 0.0, 0.0], [1.0, 1.0, 1.0]]),
            call([[0.0, 0.01, 0.0], [1.0, 1.0, 1.0]]),
            call([[0.0, -0.01, 0.0], [1.0, 1.0, 1.0]]),
            call([[0.0, 0.0, 0.01], [1.0, 1.0, 1.0]]),
            call([[0.0, 0.0, -0.01], [1.0, 1.0, 1.0]]),
            call([[0.0, 0.0, 0.0], [1.01, 1.0, 1.0]]),
            call([[0.0, 0.0, 0.0], [0.99, 1.0, 1.0]]),
            call([[0.0, 0.0, 0.0], [1.0, 1.01, 1.0]]),
            call([[0.0, 0.0, 0.0], [1.0, 0.99, 1.0]]),
            call([[0.0, 0.0, 0.0], [1.0, 1.0, 1.01]]),
            call([[0.0, 0.0, 0.0], [1.0, 1.0, 0.99]])
        ]
        mock_generate_displaced_input.assert_has_calls(expected_generate_calls, any_order=False)

    @patch('THeSeuSS.GenerateDisplacementsMolecules.GenerateDisplacements._get_no_of_atoms', return_value=2)
    @patch('THeSeuSS.GenerateDisplacementsMolecules.GenerateDisplacements._get_atom_type', return_value=([1, 2], ['H', 'O']))
    @patch('builtins.open', new_callable=mock_open)
    def test_generate_displaced_input_dftb(self, mock_file, mock_get_atom_type, mock_get_no_of_atoms):
        
        self.generator.code = 'dftb+'
        self.generator.disp_inp = 'geo.genS-001'
        sample_coord = [[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]]
        expected_path = os.path.join(self.generator.path, 'vibrations', 'geo.genS-001')

        self.generator._generate_displaced_input(sample_coord)

        mock_file.assert_called_once_with(expected_path, 'w')
        handle = mock_file()
        written_content = ''.join(call[0][0] for call in handle.write.call_args_list)

        expected_content = (
            "2 C\n"
            "H O\n"
            " 1 1 0.0 0.0 0.0\n"
            " 2 2 1.0 1.0 1.0\n"
        )

        self.assertEqual(written_content, expected_content)
