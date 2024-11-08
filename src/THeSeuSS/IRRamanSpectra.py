'''IR and Raman intensities. For periodic systems.'''

#AUTHOR: Ariadni Boziki

import numpy as np
from THeSeuSS import Constants as c
from THeSeuSS import InputsPreparation as inputs
from THeSeuSS import TwoPointCentralDifference as cendiff
from THeSeuSS import CheckPeriodicvsNonPeriodic as pervsnonper


class IntensityCalculator():

    def __init__(self, code: str, eig_vec: np.ndarray, cartesian_pol: np.ndarray, pol: np.ndarray, output_file: str, no_negfreqs: int, restart: bool, functional: str = None):

        self.code = code
        self.eig_vec = eig_vec
        self.cartesian_pol = cartesian_pol
        self.pol = pol
        self.output_file = output_file
        self.no_negfreqs = no_negfreqs
        self.restart = restart
        self.functional = functional
        self.dispersion = None
        self.supercell = None
        self.volume = None
        self.cell_dims = None
        self.dispersion = None
        self.commands = None
        self.non_periodic = None
        self.geom_input = None
        self.cart_pol_factor = None 
        self.dipole_factor = None 
        self.ir_factor = None 
        self.borh2angstrom = None
        self.borh2ang_for_IR_DFTB = None
        self.borh2ang_for_Raman_DFTB = None
        self.borh2ang_for_Raman_FHIaims = None
        self.IR_intensity = np.empty(0, dtype=np.float64)
        self.raman_activity = np.empty(0, dtype=np.float64)
        self._check_periodicity()

    def _check_periodicity(self):
        """
        Setup PeriodicvsNonPeriodic class.
        """

        check_periodic_non_periodic = pervsnonper.PeriodicvsNonPeriodic(self.code, self.cell_dims, self.output_file, self.dispersion, self.restart, self.commands, self.functional)
        self.non_periodic = check_periodic_non_periodic.check_periodic_vs_non_periodic()

    def _read_volume(self):
        """
        Reads and returns the volume of the system if the system is periodic.
        """

        central_diff = cendiff.TwoPointCentralDiff(self.code, self.output_file, self.dispersion, self.supercell, self.restart, self.functional)
        self.volume = central_diff.get_volume()

    def _set_geometry_inputs(self):
        """
        Sets the input geometry file name based on the simulation parameters (periodicity and code).
        """

        if self.non_periodic:
            if self.code == 'aims':
                self.geom_input = 'vibrations/geometry.in'
            if self.code == 'dftb+':
                self.geom_input = 'vibrations/geo.gen'
        else:
            if self.code == 'aims':
                self.geom_input = 'vibrations/geometry_backup.in'
            if self.code == 'dftb+':
                self.geom_input = 'vibrations/geo_backup.gen'

    def _set_constants(self):
        """
        Initializes the ConstantsSpectra class.
        """

        constant_initialization = c.ConstantsSpectra(self.code)
        self.cart_pol_factor, self.dipole_factor, self.ir_factor, self.borh2angstrom = constant_initialization.constants_definition()
        self.borh2ang_for_IR_DFTB = 1.0/(self.borh2angstrom*self.borh2angstrom) # 3.57106431203 - D^2/(A^2*amu) length units in dftb+: bohr.
        self.borh2ang_for_Raman_DFTB = self.borh2angstrom**4 # 0.07841599489 - (bohr^6/bohr^2 to ang^4) ang^4/amu
        self.borh2ang_for_Raman_FHIaims = self.borh2angstrom**6 #  0.02195865620442408 - (bohr^6/ang^2 to ang^4) ang^4/amu

    def IRintensity(self):
        """
        Calculates the IR intensity.
        """ 

        self._set_constants()

        if self.non_periodic:
            if self.code == 'aims':
                grad_dip = self.dipole_factor * self.cartesian_pol # D
                self.IR_intensity = np.sum(np.dot(np.transpose(grad_dip), self.eig_vec)**2, axis=0) * self.ir_factor # D^2/(A^2*amu)
            if self.code == 'dftb+':
                grad_dip = self.cartesian_pol # D
                self.IR_intensity = np.sum(np.dot(np.transpose(grad_dip), self.eig_vec)**2, axis=0) * self.ir_factor * self.borh2ang_for_IR_DFTB # D^2/(A^2*amu)
        else:
            self._read_volume()
            if self.code == 'aims':
                grad_cart_pol = self.cart_pol_factor * self.dipole_factor * self.volume * self.cartesian_pol # D
                self.IR_intensity = np.sum(np.dot(np.transpose(grad_cart_pol), self.eig_vec)**2, axis=0) * self.ir_factor # D^2/(A^2*amu)
            if self.code == 'dftb+':
                grad_cart_pol = self.cartesian_pol # D
                self.IR_intensity = np.sum(np.dot(np.transpose(grad_cart_pol), self.eig_vec)**2, axis=0) * self.ir_factor * self.borh2ang_for_IR_DFTB # D^2/(A^2*amu)

    def Ramanactivity(self):
        """
        Calculates the Raman Activity.
        """

        self._set_constants()

        alphas = np.dot(np.transpose(self.pol), self.eig_vec)

        self._set_geometry_inputs()
        geometry_processor = inputs.GeometryProcessor(self.geom_input, self.code) 
        natoms = geometry_processor.number_of_atoms()
   
        xx = np.zeros(natoms*3)
        yy = np.zeros(natoms*3)
        zz = np.zeros(natoms*3)
        xy = np.zeros(natoms*3)
        xz = np.zeros(natoms*3)
        yz = np.zeros(natoms*3)

        for step in range(natoms*3):	
            xx[step] = alphas[0,step]
            yy[step] = alphas[1,step]
            zz[step] = alphas[2,step]
            xy[step] = alphas[3,step]
            xz[step] = alphas[4,step]
            yz[step] = alphas[5,step]

        alpha = (xx + yy + zz) * (1./3)
        beta = (xx - yy)**2 + (xx - zz)**2 + (yy - zz)**2 + 6*(xy**2 + xz**2 + yz**2)

        self.raman_activity = 45 * (alpha**2) + (7./2) * beta
        if self.code == 'aims':
            self.raman_activity = self.raman_activity * self.borh2ang_for_Raman_FHIaims # ang^4/amu
        if self.code == 'dftb+':
            self.raman_activity = self.raman_activity * self.borh2ang_for_Raman_DFTB # ang^4/amu

    def spectra_calculation(self)-> [np.ndarray, np.ndarray]:
        """
        Saves the IR intensity and Raman activity.
        """

        try:
            self.IRintensity()
            print(f'IR INTENSITY (D\N{SUPERSCRIPT TWO}/(A\N{SUPERSCRIPT TWO}amu))')
            print(f'{self.IR_intensity}')
            self.IR_intensity = self.IR_intensity[self.no_negfreqs:]
            np.savetxt("IRintensity.txt", self.IR_intensity)
        except Exception as e:
            print(f'IR INTENSITY HAS NOT BEEN CALCULATED')

        try:
            self.Ramanactivity()
            print(f'RAMAN ACTIVITY (A\N{SUPERSCRIPT FOUR}/amu)')
            print(f'{self.raman_activity}')
            self.raman_activity = self.raman_activity[self.no_negfreqs:]
            np.savetxt("Ramanactivity.txt", self.raman_activity)
        except Exception as e:
            print(f'RAMAN ACTIVITY HAS NOT BEEN CALCULATED')

        return self.IR_intensity, self.raman_activity
