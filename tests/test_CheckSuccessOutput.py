#AUTHOR: Ariadni Boziki

import unittest
from unittest.mock import patch, mock_open
from THeSeuSS import CheckSuccessOutput as check

class Test_CheckOutputSuccess(unittest.TestCase):

    @patch('os.path.exists')
    @patch('builtins.open', new_callable=mock_open, read_data="Have a nice day.")
    def test_single_successful_output(self, mock_file, mock_exists):
        
        mock_exists.return_value = True
        obj = check.CheckOutputSuccess(code='aims', output='output', dispersion=True, functional='pbe')

        result = obj.single_successful_output()
        self.assertFalse(result)

    @patch('THeSeuSS.CheckSuccessOutput.CheckOutputSuccess._successful_output')
    @patch('THeSeuSS.CheckSuccessOutput.CheckOutputSuccess._successful_output_dispersion')
    def test_check_for_success_calc_before_spectra(self, mock_success_output_dispersion, mock_success_output):

        obj = check.CheckOutputSuccess(code='aims', output='output', dispersion=True, functional='pbe')
        obj.check_for_success_calc_before_spectra()
        mock_success_output.assert_called_once_with('Have a nice day.')
        mock_success_output_dispersion.assert_not_called()

        obj.functional = 'pbe0'
        obj.check_for_success_calc_before_spectra()
        mock_success_output_dispersion.assert_called_once_with('Have a nice day.')
