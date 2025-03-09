'''Restart module of the package'''

#AUTHOR: Ariadni Boziki


from THeSeuSS import CheckSuccessOutput as checksuccess
from THeSeuSS import SubmitCalculations as submit


class RestartCalculation():

    def __init__(self, code: str, output: str, dispersion: bool, restart: bool, functional: str=None, commands: str=None, cell_dims: str=None):

        self.code = code
        self.output = output
        self.dispersion = dispersion
        self.restart = restart
        self.functional = functional
        self.commands = commands
        self.cell_dims = cell_dims
        self.not_completed_calcs = None
        self._initialization_of_checkoutputsuccess()
        self._initialization_of_submitcalculations()

    def _initialization_of_checkoutputsuccess(self):
        """
        Setup the CheckOutputSuccess class. 
        """

        self.success_calcs = checksuccess.CheckOutputSuccess(self.code, self.output, self.dispersion, self.restart, self.functional)

    def directory_non_completed_calculations(self)-> list:
        """
        Returns the directories of single-point calculations that failed to complete successfully.
        """

        self.not_completed_calcs = self.success_calcs.check_for_success_calc_before_spectra()

        return self.not_completed_calcs

    def _initialization_of_submitcalculations(self):
        """
        Setup the Calculator class.             
        """

        self.calc = submit.Calculator(self.code, self.output, self.dispersion, self.restart, self.functional, self.commands, self.cell_dims)

    def restart_calculations(self):
        """
        Restart the submission of single-point calculations that failed to complete successfully.
        """

        self.calc.submit_jobs_in_parallel()
