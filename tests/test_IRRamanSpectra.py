#Author: Ariadni Boziki


import unittest
from unittest.mock import patch, MagicMock
import numpy as np
import os
from THeSeuSS import IRRamanSpectra as IRRaman


class Test_IntensityCalculator(unittest.TestCase):

    def setUp(self):

        self.eig_vec = np.array([
            [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
            [2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0],
            [3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0],
            [4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0],
            [5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0],
            [6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0],
            [7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0],
            [8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0],
            [9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0]
            ])
        self.cartesian_pol = np.array([[1.0, 1.0, 1.0], 
            [2.0, 2.0, 2.0],
            [3.0, 3.0, 3.0],
            [4.0, 4.0, 4.0],
            [5.0, 5.0, 5.0],
            [6.0, 6.0, 6.0],
            [7.0, 7.0, 7.0],
            [8.0, 8.0, 8.0],
            [9.0, 9.0, 9.0]
            ])
        self.pol = np.array([[1.0, 1.0, 1.0, 1.0, 1.0, 1.0], 
            [2.0, 2.0, 2.0, 2.0, 2.0, 2.0],
            [3.0, 3.0, 3.0, 3.0, 3.0, 3.0],
            [4.0, 4.0, 4.0, 4.0, 4.0, 4.0],
            [5.0, 5.0, 5.0, 5.0, 5.0, 5.0],
            [6.0, 6.0, 6.0, 6.0, 6.0, 6.0],
            [7.0, 7.0, 7.0, 7.0, 7.0, 7.0],
            [8.0, 8.0, 8.0, 8.0, 8.0, 8.0],
            [9.0, 9.0, 9.0, 9.0, 9.0, 9.0]
            ])
        self.output_file = "output"
        self.no_negfreqs = 2

    @patch('THeSeuSS.CheckPeriodicvsNonPeriodic.PeriodicvsNonPeriodic')
    @patch('THeSeuSS.Constants.ConstantsSpectra')
    def test_IRintensity_non_periodic(self, mock_constants, mock_pnp):
        
        mock_pnp_instance = mock_pnp.return_value
        mock_pnp_instance.check_periodic_vs_non_periodic.return_value = True

        mock_constants_instance = mock_constants.return_value
        mock_constants_instance.constants_definition.return_value = (1.0, 1.0, 1.0, 1.0)

        intensity_calculator = IRRaman.IntensityCalculator(
            code='aims',
            eig_vec=self.eig_vec,
            cartesian_pol=self.cartesian_pol,
            pol=self.pol,
            output_file=self.output_file,
            no_negfreqs=self.no_negfreqs,
            restart=False
        )

        intensity_calculator.IRintensity()
        expected_intensity = np.sum(np.dot(np.transpose(self.cartesian_pol * 1.0), self.eig_vec)**2, axis=0) * 1.0
        np.testing.assert_array_almost_equal(intensity_calculator.IR_intensity, expected_intensity)

    @patch('THeSeuSS.TwoPointCentralDifference.TwoPointCentralDiff')
    @patch('THeSeuSS.CheckPeriodicvsNonPeriodic.PeriodicvsNonPeriodic')
    @patch('THeSeuSS.Constants.ConstantsSpectra')
    def test_IRintensity_periodic(self, mock_constants, mock_pnp, mock_cendiff):
        
        mock_pnp_instance = mock_pnp.return_value
        mock_pnp_instance.check_periodic_vs_non_periodic.return_value = False

        mock_constants_instance = mock_constants.return_value
        mock_constants_instance.constants_definition.return_value = (1.0, 1.0, 1.0, 1.0)
        
        mock_cendiff_instance = mock_cendiff.return_value
        mock_cendiff_instance.get_volume.return_value = 100

        intensity_calculator = IRRaman.IntensityCalculator(
            code='dftb+',
            eig_vec=self.eig_vec,
            cartesian_pol=self.cartesian_pol,
            pol=self.pol,
            output_file=self.output_file,
            no_negfreqs=self.no_negfreqs,
            restart=False
        )

        intensity_calculator.IRintensity()
        grad_cart_pol = self.cartesian_pol
        expected_intensity = np.sum(np.dot(np.transpose(grad_cart_pol), self.eig_vec)**2, axis=0)
        np.testing.assert_array_almost_equal(intensity_calculator.IR_intensity, expected_intensity)

    @patch('THeSeuSS.InputsPreparation.GeometryProcessor')
    @patch('THeSeuSS.Constants.ConstantsSpectra')
    @patch('THeSeuSS.CheckPeriodicvsNonPeriodic.PeriodicvsNonPeriodic')
    def test_Ramanactivity(self, mock_pnp, mock_constants, mock_geometry):
        
        mock_pnp_instance = mock_pnp.return_value
        mock_pnp_instance.check_periodic_vs_non_periodic.return_value = True 

        mock_constants_instance = mock_constants.return_value
        mock_constants_instance.constants_definition.return_value = (1.0, 1.0, 1.0, 1.0)

        mock_geometry_instance = mock_geometry.return_value
        mock_geometry_instance.number_of_atoms.return_value = 3

        intensity_calculator = IRRaman.IntensityCalculator(
            code='aims',
            eig_vec=self.eig_vec,
            cartesian_pol=self.cartesian_pol,
            pol=self.pol,
            output_file=self.output_file,
            no_negfreqs=self.no_negfreqs,
            restart=False
        )

        intensity_calculator.Ramanactivity()
        alphas = np.dot(np.transpose(self.pol), self.eig_vec)
        
        xx = np.zeros(mock_geometry_instance.number_of_atoms.return_value*3)
        yy = np.zeros(mock_geometry_instance.number_of_atoms.return_value*3)
        zz = np.zeros(mock_geometry_instance.number_of_atoms.return_value*3)
        xy = np.zeros(mock_geometry_instance.number_of_atoms.return_value*3)
        xz = np.zeros(mock_geometry_instance.number_of_atoms.return_value*3)
        yz = np.zeros(mock_geometry_instance.number_of_atoms.return_value*3)

        for step in range(mock_geometry_instance.number_of_atoms.return_value*3):	
            xx[step] = alphas[0,step]
            yy[step] = alphas[1,step]
            zz[step] = alphas[2,step]
            xy[step] = alphas[3,step]
            xz[step] = alphas[4,step]
            yz[step] = alphas[5,step]

        alpha = (xx + yy + zz) * (1./3)
        beta = (xx - yy)**2 + (xx - zz)**2 + (yy - zz)**2 + 6*(xy**2 + xz**2 + yz**2)
        expected_raman_activity = 45 * (alpha**2) + (7./2) * beta

        np.testing.assert_array_almost_equal(intensity_calculator.raman_activity, expected_raman_activity)
