
"""Classes for reaction rate coefficients."""

import numpy
from pychemkin.parsers import SQLParser
from pychemkin.pychemkin_errors import PyChemKinError


class RxnCoefficientBase:
    """Base class for reaction rate coefficient.
    
    Attributes:
    ===========
    k : float
        Reaction rate coefficient, initialized to None

    Methods:
    ========
    get_coefficient() : Calculates and returns the reaction
        rate coefficient. To be implemented by subclasses.
    """
    def __init__(self):
        self.k = None

    def __repr__(self):
        return 'RxnCoefficientBase()'

    def get_coefficient(self):
        """Calculates and returns the reaction rate coefficient.
        
        Notes:
        ======
        - Not implemented in this base class but should be
            implemented by its subclasses.
        """
        raise NotImplementedError('Subclass must implement this method!')


class ConstantCoefficient(RxnCoefficientBase):
    """Class for constant reaction rate coefficient.

    Attributes:
    ===========
    k : float, required
        Constant reaction rate coefficient
    """
    def __init__(self, k):
        """Initializes ConstantCoefficient.

        Notes:
        ======
        - Raises a ValueError if inputed constant reaction
            rate coefficient is negative.
        """
        super().__init__()
        if k <= 0:
            raise ValueError("Reaction rate coefficients must be positive!")
        self.k = k

    def __repr__(self):
        return 'ConstantCoefficient(k={})'.format(self.k)

    def get_coefficient(self):
        """Returns the reaction rate coefficient."""
        return self.k


class ArrheniusCoefficient(RxnCoefficientBase):
    """Class for Arrhenius reaction rate coefficient.

    Attributes:
    ===========
    A: float, required
        Arrhenius prefactor
    E: float, required
        Activation energy
    T: float, required
        Temperature, in Kelvin
    R: float, optional, default value = 8.314
        Ideal gas constant
        Default value of R should not be changed except for unit conversion
    k: float
        Reaction rate coefficient, calculated by get_coefficient()
    """
    def __init__ (self, A, E, T, R=8.314):
        """Initializes ArrheniusCoefficient.

        Notes:
        ======
        - Raises a ValueError if ...
            Arrhenius prefactor is non-positive
            Temperature is non-positive
            Ideal gas constant is non-positive
        """
        super().__init__()
        if A <= 0.0:
            raise ValueError(
                  "A = {0:18.16e}:  Non-positive Arrhenius prefactor is "
                  "prohibited!".format(A))

        if T <= 0.0:
            raise ValueError(
                  "T = {0:18.16e}:  Non-positive temperatures are "
                  "prohibited!".format(T))

        if R <= 0.0:
            raise ValueError(
                  "R = {0:18.16e}:  Non-positive ideal gas constant is "
                  "prohibited!".format(R))

        self.A = A
        self.E = E
        self.T = T
        self.R = R
        self.k = None

    def __repr__ (self):
        return ('ArrheniusCoefficient(A={}, E={}, T={}, R={})'.format(
                    self.A, self.E, self.T, self.R))

    def get_coefficient(self):
        """"Returns the reaction rate coefficient.
        
        Notes:
        ======
        - Raises an OverflowError if resulting coefficient
            is too large/small (numerical instability)
        """
        try:
            self.k = self.A * numpy.exp(-self.E / (self.R * self.T))
        except Warning:
            raise OverflowError("The result is too large/small.")
        return self.k


class ModifiedArrheniusCoefficient(ArrheniusCoefficient):
    """Class of Modified Arrhenius reaction rate coefficient.

    Attributes:
    ===========
    A: float, required
        Arrhenius prefactor
    b: float, required
        Modified Arrhenius parameter
    E: float, required
        Activation energy
    T: float, required
        Temperature, in Kelvin
    R: float, optional, default value = 8.314
        Ideal gas constant
        Default value of R should not be changed except for unit conversion
    k: float
        Reaction rate coefficient
    """
    def __init__ (self, A, b, E, T, R=8.314):
        """Initializes ModifiedArrheniusCoefficient.

        Notes:
        ======
        - Raises a ValueError if ...
            Arrhenius prefactor is non-positive
            Temperature is non-positive
            Ideal gas constant is non-positive
            Modified Arrhenius parameter b is not real
        """
        super().__init__(A, E, T, R)

        if not numpy.isreal(b):
            raise ValueError('Modified Arrhenius parameter b must be real!')
        self.b = b
        self.k = None

    def __repr__ (self):
        return('ModifiedArrheniusCoefficient(A={}, b={}, E={}, T={}, R={})'.format(
                    self.A, self.b, self.E, self.T, self.R))

    def get_coefficient(self):
        """"Returns the reaction rate coefficient.
        
        Notes:
        ======
        - Raises an OverflowError if resulting coefficient
            is too large/small (numerical instability)
        """
        try:
            self.k = (self.A * numpy.power(self.T, self.b) *
                      numpy.exp(-self.E / (self.R * self.T)))
        except Warning:
            raise OverflowError("The result is too large/small.")
        return self.k


class BackwardRxnCoefficientBase:
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


class NASA7BackwardCoeffs(BackwardRxnCoefficientBase):
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
        """
        super().__init__()
        self.nui = nui
        self.nasa7_coeffs = nasa7_coeffs
        self.p0 = p0
        self.R = R
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
        kf : numpy.ndarray[float]
            array of forward reaction rate coefficients for each specie in reaction
        T : float
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

        # Prefactor in k_e (equilibrium coefficient)
        fact = self.p0 / self.R / T

        ke = (fact ** self.gamma) * numpy.exp(delta_G_over_RT)

        return kf / ke
