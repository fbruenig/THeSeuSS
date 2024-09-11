'''It applies the spectrum broadening and prepares the IR and
Raman plots.'''

#AUTHOR: Ariadni Boziki

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import cauchy, norm


class SpectraPlotter():

    def __init__(self, freq: np.ndarray, IRintensity: np.ndarray, Ramanactivity: np.ndarray, broadening: str, fwhm: float):

        self.freq = freq
        self.IRintensity = IRintensity
        self.Ramanactivity = Ramanactivity
        self.broadening = broadening
        self.fwhm = fwhm

    def _spectrum_broadening(self, x: np.ndarray, x0: np.ndarray, y0: np.ndarray, fit_points: bool =True)-> np.ndarray:
        """
        Applies spectrum broadening using either Gaussian or Lorentzian distribution.
        """
        
        if self.broadening == 'gaussian':
            distribution = norm
        elif self.broadening == 'lorentzian':
            distribution = cauchy
            
        s = np.sum([yp * distribution.pdf(x, xp, scale=self.fwhm) for xp, yp in zip(x0, y0)], axis=0)

        if fit_points:
            s_max = np.max(s)
            if s_max == 0.0:
                s_max = 1.0
            return s * np.max(y0) / s_max

    def plot_spectrum_IR(self):
        """
        Prepares the plot for the IR spectrum.
        """

        x = np.linspace(10,6000, num=500, endpoint=True)
        y_IR = self._spectrum_broadening(x, self.freq, self.IRintensity, fit_points=True)

        np.savetxt("IRspectrum.txt", np.column_stack([x, y_IR]))	

        plt.title("IR spectrum")
#       plt.xlabel(r'Frequency (THz)')
        plt.xlabel(r'Frequency (cm$^{\rm -1}$)')
        plt.ylabel(r'IR Intensity (D$^{\rm 2}$/${\rm \AA}$$^{\rm 2}$ amu)')
        plt.plot(x, y_IR, color='r')
        plt.stem(self.freq, self.IRintensity, markerfmt=',', basefmt=" ") 
        plt.savefig("IR_spectrum.png")
        plt.show()

    def plot_spectrum_Raman(self):
        """
        Prepares the plot for the Raman spectrum.
        """

        x = np.linspace(10,6000, num=500, endpoint=True)
        y_Raman = self._spectrum_broadening(x, self.freq, self.Ramanactivity, fit_points=True)	

        np.savetxt("Ramanspectrum.txt", np.column_stack([x, y_Raman]))

        plt.title("Raman spectrum") 
#       plt.xlabel(r'Frequency (THz)') 
        plt.xlabel(r'Frequency (cm$^{\rm -1}$)')
        plt.ylabel(r'Raman Activity (${\rm \AA}$$^{\rm 4}$/amu)')
        plt.plot(x, y_Raman, color='r')
        plt.stem(self.freq, self.Ramanactivity, markerfmt=',', basefmt=" ") 
        plt.savefig("Raman_spectrum.png")
        plt.show()

    def plot_spectra(self):
        """
        Plots the IR and Raman spectra.
        """

        plt.figure()
        self.plot_spectrum_IR()
        plt.figure()
        self.plot_spectrum_Raman()
