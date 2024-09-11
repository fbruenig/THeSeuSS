'''Constants'''

#AUTHOR: Ariadni Boziki

#import scipy as sc
import math
from scipy import constants


class ConstantsSpectra():

    def __init__(self, code: str):

        self.code = code
        self.cart_pol_factor = None
        self.dipole_factor = None
        self.ir_factor = None
        self.avogadro = None
        self.borh2angstrom = None
        self.hubbard_derivates = None
        self.max_angular_mom = None

    def constants_definition(self)-> [float, float, float, float]:
        """
        Constants
        """

        self.avogadro = 6.0221413e+23
        self.borh2angstrom = 0.529177249        
        self.cart_pol_factor = 0.062415091 # C/(m^2) -> e/(A^2)
        self.dipole_factor = 4.80320546720059098640 # eAng -> D

        self.pi = constants.pi
        self.avogadro = constants.Avogadro
        self.speed_of_light = constants.speed_of_light
        self.electric_constant = constants.epsilon_0

        factor = (1/(4*self.pi*self.electric_constant))*((self.avogadro*self.pi)/(3*self.speed_of_light))

        self.ir_factor = 1

        return self.cart_pol_factor, self.dipole_factor, self.ir_factor, self.borh2angstrom

    def hubbard_derivates_dict(self):
        """
        Dictionary containing all atomic Hubbard derivates (atomic units). 
        Retrieved by https://dftb.org/parameters/download/3ob/3ob-3-1-cc.
        """

        self.hubbard_derivates = {
            "Br": -0.0573,
            "C": -0.1492,
            "Ca": -0.0340,
            "Cl": -0.0697,
            "F": -0.1623,
            "H": -0.1857,
            "I": -0.0433,
            "K": -0.0339,
            "Mg": -0.02,
            "N": -0.1535,
            "Na": -0.0454,
            "O": -0.1575,
            "P": -0.14,
            "S": -0.11,
            "Zn": -0.03
        }

        return self.hubbard_derivates

    def max_angular_mom_dict(self):
        """
        Dictionary containing maximum angular momenta.
        """

        self.max_angular_mom = {
            'Br': "d",
            'C': "p",
            'Ca': "p",
            'Cl': "d",
            'F': "p",
            'H': "s",
            'I': "d",
            'K': "p",
            'Mg': "p",
            'N': "p",
            'Na': "p",
            'O': "p",
            'P': "d",
            'S': "d",
            'Zn': "d"
        }

        return self.max_angular_mom
