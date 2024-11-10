#AUTHOR: Ariadni Boziki

import unittest
import os
import shutil
import numpy as np
from unittest.mock import patch, MagicMock
from THeSeuSS import CheckPeriodicvsNonPeriodic as checkpnp


class Test_PeriodicvsNonPeriodic(unittest.TestCase):

    def setUp(self):
        
        with open('geometry.in', 'w') as f:
            f.write("lattice_vector 1 1 1\n")

    def tearDown(self):
        
        os.remove('geometry.in')
        if os.path.exists('vibrations'):
            shutil.rmtree('vibrations')

    def test_class_initialization(self):

        pnp = checkpnp.PeriodicvsNonPeriodic(
            code='aims',
            cell_dims='1 1 1',
            output_file='output',
            dispersion=True,
            restart=False,
            commands='commands'
        )

        self.assertEqual(pnp.code, 'aims')
        self.assertEqual(pnp.cell_dims, '1 1 1')
        self.assertEqual(pnp.output_file, 'output')
        self.assertTrue(pnp.dispersion)
        self.assertEqual(pnp.commands, 'commands')
        self.assertFalse(pnp.non_periodic)

    def test_check_periodic_vs_non_periodic(self):

        pnp = checkpnp.PeriodicvsNonPeriodic(
            code='aims',
            cell_dims='1 1 1',
            output_file='output',
            dispersion=True,
            restart=False,
            commands='commands'
        )
        
        result = pnp.check_periodic_vs_non_periodic()
        self.assertFalse(result)

    @patch('THeSeuSS.SubmitCalculations.PhonopyCalculator.submit_phonopy_displacements')
    def test_decision_submission_periodic(self, mock_submit_phonopy):

        pnp = checkpnp.PeriodicvsNonPeriodic(
            code='aims',
            cell_dims='1 1 1',
            output_file='output',
            dispersion=True,
            restart=False,
            commands='commands'
        )

        pnp.non_periodic = False
        pnp.decision_submission_cell()

        mock_submit_phonopy.assert_called_once()

    @patch('THeSeuSS.CheckSuccessOutput.CheckOutputSuccess.single_successful_output', return_value=False)
    @patch('THeSeuSS.MapAtoms.spglibProcessor.get_international_space_group_number')
    def test_periodic_space_group_calc(self, mock_successful_output, mock_get_space_group):

        pnp = checkpnp.PeriodicvsNonPeriodic(
            code='aims',
            cell_dims='1 1 1',
            output_file='output',
            dispersion=True,
            restart=False,
            commands='commands'
        )

        pnp.non_periodic = False
        pnp.periodic_space_group_calc()

        mock_successful_output.assert_called_once()
        mock_get_space_group.assert_called_once()

    @patch('THeSeuSS.SubmitCalculations.PhonopyCalculator.disp_forces_dataset_dyn_matrix')
    def test_dyn_mat_for_eigvec_eig_val_freq(self, mock_dynamical_matrix):

        mock_dynamical_matrix.return_value = np.array([[1, 2], [3, 4]]) 

        pnp = checkpnp.PeriodicvsNonPeriodic(
            code='aims',
            cell_dims='1 1 1',
            output_file='output',
            dispersion=True,
            restart=False,
            commands='commands'
        )

        pnp.non_periodic = False

        # Mock the eigvec, eigvals, and frequencies returned by read_eig_vec_phonopy
        pnp.vibrational_freq = MagicMock()
        pnp.vibrational_freq.read_eig_vec_phonopy.return_value = ('eigvecs', 'eigvals', np.array([100, 200, -50]))

        pnp.dyn_mat_for_eigvec_eig_val_freq()

        mock_dynamical_matrix.assert_called_once()
