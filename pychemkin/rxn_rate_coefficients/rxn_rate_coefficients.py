
"""Classes for reaction rate coefficients."""

import numbers
import numpy
import warnings
from pychemkin.parsers import SQLParser


class ReactionCoeff:
    """Class for (forward) reaction rate coefficient."""
    def __init__(self, k_parameters, T=None):
        """Initializes forward reaction rate coefficient.

        Args:
        =====
        k_parameters : dictionary, required
            Dictionary of parameters to compute reaction rate coefficient
        T : int or float
            Temperature of reaction, in Kelvin (default: None,
            in case rxn rate coefficient is temperature-independent)

        Attributes:
        ===========
        k : float
            Computed forward reaction rate coefficient
        """
        self.k_parameters = k_parameters
        self.T = T
        self.k = self.compute_coefficient(self.k_parameters, self.T)

    def compute_coefficient(self, k_parameters, T):
        """Computes reaction rate coefficient using passed parameters.

        Args:
        =====
        k_parameters : dictionary
            Dictionary of parameters to compute rxn rate coefficient
        T : int or float
            Temperature of the reaction, in Kelvin
        
        Returns:
        ========
        k : int or float
            Reaction rate coefficient

        NOTES:
        ------
        PRE:
            - Raise ValueError if customized reaction rate coefficient depends on T
        POST:
            - Raises NotImplementedError if dictionary of k parameters is not recognized
            - Options to alter values of R (to change units) but strongly discouraged
            - Raises ValueError if valid T not inputed/set for Arrhenius and modified Arrhenius
        """
        # Constant
        if "k" in k_parameters and len(k_parameters) == 1:
            return self.const(k_parameters['k'])
        
        # Arrhenius
        elif (("A" in k_parameters and "E" in k_parameters and "b" not in k_parameters) and
              (len(k_parameters)==2 or len(k_parameters)==3)):

            if T == None:
                raise ValueError("Temperature has not been set in the reaction. Please use set_temperature() method.")

            if "R" in k_parameters:
                return self.arr(A=k_parameters['A'],
                                E=k_parameters['E'],
                                T=T,
                                R=k_parameters['R'])
            else:
                if len(k_parameters) == 2:
                    return self.arr(A=k_parameters['A'],
                                    E=k_parameters['E'],
                                    T=T)
                else:
                    raise NotImplementedError("The combination of parameters entered is not supported "
                                              "for the calculation of reaction rate coefficient.")

        
        # Modified Arrhenius
        elif (("A" in k_parameters and "E" in k_parameters  and "b" in k_parameters) and
              (len(k_parameters) == 3 or len(k_parameters) == 4)):
            if T == None:
                raise ValueError("Temperature has not been set in the reaction. Please use set_temperature() method.")

            if "R" in k_parameters:
                return self.mod_arr(A=k_parameters['A'],
                                    E=k_parameters['E'],
                                    R=k_parameters['R'],
                                    b=k_parameters['b'],
                                    T=T)
            else:
                if len(k_parameters) == 3:
                    return self.mod_arr(A=k_parameters['A'],
                                        E=k_parameters['E'],
                                        b=k_parameters['b'],
                                        T=T)
                else:
                    raise NotImplementedError("The combination of parameters entered is not supported "
                                              "for the calculation of reaction rate coefficient.")

        else:
            raise NotImplementedError("The combination of parameters entered is not supported for the calculation of Reaction Rate Coefficient.")


    def const(self, k):
        """Returns constant reaction rate coefficient.

        INPUTS:
        -------
        k : numeric type 
            Constant reaction rate coefficient

        RETURNS:
        --------
        k : numeric type
            Constant reaction rate coefficients.

        NOTES:
        ------
        POST:
            - Raises ValueError if k is non-positive!
        """
        if k <= 0:
            raise ValueError("Reaction rate must be positive.")
        return k

    def arr(self, A, E, T, R=8.3144598):
        """Returns Arrhenius reaction rate coefficient.

        INPUTS:
        -------
        A: float, strictly positive, no default value
           The Arrhenius prefactor
        E: float, no default value
           The Arrhenius parameter
        T: float, strictly positive, no default value
           Temperature T, asuuming a Kelvin scale
        R: float, default value is 8.314, cannot be changed except to convert units
           The ideal gas constant

        RETURNS:
        --------
        k: Arrhenius reaction rate coefficients k,
           floats
           unless A or T is not postive
           in which case a ValueError exception is raised

        NOTES:
        ------
        POST:
            - Raises ValueError if A, T, or R is non-positive
            - Raises Warning if user changes value of R
        """
        if (A <= 0):
            raise ValueError("Arrhenius prefactor 'A' must be positive.")

        if (T <= 0):
            raise ValueError("Temperatures 'T' must be positive.")

        if (R <= 0):
            raise ValueError("Gas constant 'R' must be positive.")

        if not numpy.isclose(R, 8.3144598):
            warnings.warn("Please do not change the value of"
                          " universal gas constant 'R' unless you are converting units.")
        k = A * numpy.exp(-E / R / T)
        return k

    def mod_arr(self, A, b, E, T, R=8.3144598):
        """Returns Arrhenius reaction rate coefficients k.

        INPUTS:
        -------
        A: float, strictly positive, no default value
           The Arrhenius prefactor
        b: real, no default value
           Modified Arrhenius parameter
        E: float, no default value
           The Arrhenius parameter
        T: float, strictly positive, no default value
           Temperature T, asuuming a Kelvin scale
        R: float, default value is 8.314, cannot be changed except to convert units
           The ideal gas constant

        RETURNS:
        --------
        k: Arrhenius reaction rate coefficients k,
           floats
           unless A or T is not postive
           in which case a ValueError exception is raised
           Or b is not a real number
           in which case a TypeError exception is raised

        NOTES:
        ------
        POST:
            - Raises ValueError if A, T, or R is non-positive
            - Raises TypeError if b is not real
            - Raises Warning if user changes value of R
        """

        if (A <= 0):
            raise ValueError("Parameter 'A' must be positive.")
        
        if (T <= 0):
            raise ValueError("Parameter 'T' must be positive.")
        
        if (isinstance(b, numbers.Real)) == False:
            raise TypeError("Parameter 'b' must be a real number.")
        
        if (R <= 0):
            raise ValueError("Gas constant 'R' must be positive.")

        if not numpy.isclose(R, 8.3144598):
            warnings.warn("Please do not change the value of"
                          " universal gas constant 'R' unless you are converting units.")

        k = A * T ** b * numpy.exp(-E / R / T)
        return k


class BackwardReactionCoeff:
    """Base class for backward reaction rate coefficients.
    
    Methods:
    ========
    compute_backward_coeffs() : Calculates and returns the backward reaction
        rate coefficient. To be implemented by subclasses.
    """
    def __init__(self):
        pass

    def compute_backward_coeffs(self):
        """Calculates and returns the backward eaction rate coefficient.
        
        Notes:
        ======
        - Not implemented in this base class but should be
            implemented by its subclasses.
        """
        raise NotImplementedError('Subclass must implement this method!')


class NASA7BackwardCoeffs(BackwardReactionCoeff):
    """Class for computing backward reaction rate
    coefficients for reversible reactions."""
    def __init__(self, nui, nasa7_coeffs, p0=1e5, R=8.3144598):
        """Initializes backward coefficients using 7th order NASA polynomials.

        Args:
        =====
        nui : numpy.ndarray
            stoichiometric coefficient difference (stoich_products - stoich_reactants)
                for a single reversible reaction
        nasa7_coeffs : numpy.ndarray
            NASA polynomial coefficients (from appropriate temperature range)
                corresponding to species in reversible reaction

        Attributes:
        ===========
        p0 : float
            pressure of reaction, in Pascals
        R : float
            gas constant, in J / mol / K
        gamma : numpy.ndarray
            sum of stoichiometric coefficient difference 

        Notes:
        ======
        PRE:
            - Assuming ordering of stochiometric coefficient difference consistent
            with ordering of species in all involved functions
        """
        super().__init__()
        self.nui = nui
        self.nasa7_coeffs = nasa7_coeffs
        self.p0 = p0
        self.R = R
        if not numpy.isclose(self.R, 8.3144598):
            warnings.warn("Please do not change the value of"
                          " universal gas constant 'R' unless you are converting units.")
        self.gamma = numpy.sum(self.nui)

    def H_over_RT(self, T):
        """Helper function that returns the enthalpy of each specie given by
        the NASA polynomials.

        Args:
        =====
        T : float, required
            temperature of reaction

        Returns:
        ========
        H_RT : numpy.ndarray
            enthalpy values for each specie

        Notes:
        ======
        PRE:
            - Raises ValueError if inputed temperature is non-positive
        """
        if T <= 0:
            raise ValueError("Temperature has to be a positive value!")

        a = self.nasa7_coeffs
        H_RT = (a[:, 0] + (0.5 * a[:, 1] * T) + (a[:, 2] * T ** 2.0) / 3.0
                + (a[:, 3] * T ** 3.0) / 4.0 + (a[:, 4] * T ** 4.0) / 5.0
                + a[:, 5] / T)
        return H_RT

    def S_over_R(self, T):
        """Helper function that returns the entropy of each specie given by
        the NASA polynomials.

        Args:
        =====
        T : float
            temperature of reaction

        Returns:
        ========
        S_R : numpy.ndarray
            entropy values for each specie

        Notes:
        ======
        PRE:
            - Raises ValueError if inputed temperature is non-positive
        """
        if T <= 0:
            raise ValueError("Temperature has to be a positive value!")

        a = self.nasa7_coeffs
        S_R = (a[:, 0] * numpy.log(T) + a[:, 1] * T + (a[:, 2] * T ** 2.0) / 2.0
               + (a[:, 3] * T ** 3.0) / 3.0 + (a[:, 4] * T ** 4.0) / 4.0 + a[:, 6])
        return S_R

    def compute_backward_coeffs(self, kf, T):
        """Returns the backward reaction rate
        coefficient for each specie.

        Args:
        =====
        kf : numpy.ndarray[float], required
            array of forward reaction rate coefficients for each specie in reaction
        T : float, required
            temperature of reaction

        Returns:
        ========
        kb : backward reaction rate coefficient for each specie

        Notes:
        ======
        PRE:
            - Raises ValueError if inputed temperature is non-positive
        """
        if T <= 0:
            raise ValueError("Temperature has to be a positive value!")

        # Change in enthalpy and entropy for each reaction
        delta_H_over_RT = numpy.dot(self.nui, self.H_over_RT(T))
        delta_S_over_R = numpy.dot(self.nui, self.S_over_R(T))

        # Negative of change in Gibbs free energy for each reaction
        delta_G_over_RT = delta_S_over_R - delta_H_over_RT

        # Prefactor for equilibrium coefficient
        fact = self.p0 / self.R / T

        # Compute equilibrium coefficient
        ke = (fact ** self.gamma) * numpy.exp(delta_G_over_RT)
        return kf / ke
