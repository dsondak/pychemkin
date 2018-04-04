
"""Classes for reactions."""

import numpy

from pychemkin.rxn_rate_coefficients.rxn_rate_coefficients import *


class ElementaryReactionError(Exception):
    """Class for ElementaryReaction-related errors."""
    pass

class ElementaryReaction:
    """Base class for an elementary reaction.

    Note:
    =====
        - This class is meant to serve as a framework
        for specific types of elementary reactions"""
    def __init__(self, rxn_type, is_reversible, rxn_equation, species_list,
                 rate_coeffs_components, rate_coeffs_type, reactant_stoich_coeffs, product_stoich_coeffs):
        """Initializes an elementary reaction.
    
        Args:
        =====
        rxn_type : str, required
            Type of reaction (e.g. "Elementary")
        is_reversible : bool, required
            True if reaction is reversible
        rate_coeffs_components : dict, required
            Dictionary of components (e.g. 'A', 'b', 'E')
            to compute reaction rate coefficients
        rxn_rate_coeffs_type : str, required
            Type of reaction rate coefficient (e.g. 'constant', 'arrhenius', 'modified arrhenius')
        rxn_equation : str, required
            String representation of reaction equation
        reactant_stoich_coeffs : dict[int], required
            Dictionary for reactant stoichiometric coefficients
        product_stoich_coeffs : dict[int], required
            Dictionary for product stoichiometric coefficients
        species_list : list, required
            List of chemical species from original xml file (useful for ordering)

        Attributes:
        ===========
        unique_species : list
            List of unique chemical species
        temperature : int or float
            Temperature of reaction, in Kelvin
        concentrations : list
            Concentrations of species involved in reaction
        rxn_rate_coeff : float
            Reaction rate coefficient
        """
        self.rxn_type = rxn_type
        self.is_reversible = is_reversible
        self.rate_coeffs_components = rate_coeffs_components
        self.rate_coeffs_type = rate_coeffs_type
        self.rxn_equation = rxn_equation
        self.reactant_stoich_coeffs = reactant_stoich_coeffs
        self.product_stoich_coeffs = product_stoich_coeffs
        self.unique_species = self.get_unique_species()
        self.species_list = species_list

        # Pad non-participating species with 0
        for specie in self.species_list:
            if specie not in self.reactant_stoich_coeffs:
                self.reactant_stoich_coeffs[specie] = 0
            if specie not in self.product_stoich_coeffs:
                self.product_stoich_coeffs[specie] = 0

        self.temperature = None
        self.concentrations = []
        self.rxn_rate_coeff = None

    def __str__(self):
        """Returns user-friendly string representation of reaction.

        Returns:
        ========
        info : str
            String representation of reaction (reaction equation)
        """
        info = "Reaction : {}".format(self.rxn_equation)
        return info

    def __len__(self):
        """Returns number of unique species in reaction.

        Returns:
        ========
        n_species : int
            Number of unique species involved in the reaction
        """
        return len(self.unique_species)

    def get_unique_species(self):
        """Helper function to return unique species involved
        in the reaction.
        
        Returns:
        ========
        unique_species : list
            List of unique species in reaction
        """
        reactant_species = self.reactant_stoich_coeffs.keys()
        product_species = self.product_stoich_coeffs.keys()
        unique_species = list(set(reactant_species) | set(product_species))
        return unique_species

    def set_temperature(self, T):
        """Sets temperature of the reaction.

        Args:
        =====
        T : float, required
            Temperature of reaction

        Notes:
        ======
        POST:
            - Updates self.temperature
            - Raises ValueError if inputed temperature is non-positive
        """
        if T <= 0:
            raise ValueError("Temperature has to be a positive value!")
        self.temperature = T

    def set_concentrations(self, X):
        """Sets concentrations of the reaction from a dictionary.

        Args:
        =====
        X : dict, required
            Dictionary with species and corresponding concentrations

        Notes:
        ======
        PRE:
            - Raises KeyError if input dictionary of species and concentrations
                contains name of species that was not in original xml file
            - Raises ValueError if any of concentrations are negative
        """
        try:
            ordered_concentrations = self.get_ordered_list(X)
        except KeyError:
            raise KeyError("Invalid concentration entered!")

        if (numpy.array(ordered_concentrations) < 0).any():
            raise ValueError("Invalid (e.g. negative) concentration(s).")
        self.concentrations = numpy.array(ordered_concentrations)

    def set_concentrations_from_array(self, X):
        """Sets the concentrations of the reaction from an array.

        Args:
        =====
        X : numpy.array, required
        """
        self.concentrations = X

    def get_ordered_list(self, dictionary):
        """Helper/utilities function to output a list of a dictionary's keys
        in the order of species_list. This is to ensure a consistent ordering
        scheme (useful for ordering concentrations and
        stoichiometric coefficients).

        Args:
        =====
        dictionary : dict, required
            dictionary to order 

        Returns:
        ========
        list_of_interest : list
            list of dictionary's keys in order of species_list
        """
        index_map = ({value: index for index, value
                     in enumerate(self.species_list)})
        sorted_tuple_list = sorted(dictionary.items(),
                                   key=lambda pair: index_map[pair[0]])
        list_of_interest = [element[1] for element in sorted_tuple_list]
        return list_of_interest

    def compute_reaction_rate_coeff(self, T=None):
        """Computes reaction rate coefficient of reaction.

        Returns:
        ========
        k : numeric type (or list of numeric type)
            Reaction rate coefficient

        Notes:
        ======
        POST:
            - Raises NotImplementedError (user must define this function)
        """
        raise NotImplementedError

    def compute_progress_rate(self, T=None):
        """Computes progress rate of reaction.

        Returns:
        ========
        omega_array : numpy.ndarray
            Array of progress rates of reaction

        Notes:
        ======
        POST:
            - Raises NotImplementedError (user must define this function)
        """
        raise NotImplementedError
        
    def compute_reaction_rate(self, T=None):
        """Computes reaction rate of reaction.

        Returns:
        ========
        rxn_rate_array: numpy.ndarray
            Array of reaction rates of reaction

        Notes:
        ======
        POST:
            - Raises ValueError if reactant or product stoichiometric coefficients
                are negative
            - Raises ReactionError if compute_reaction_rate_coeff nor compute_progress_rate
                has not been implemented. This Reaction class is meant to serve as a framework
                for specific types of reactions!
        """
        reactant_stoich_coeffs = numpy.array(self.get_ordered_list(self.reactant_stoich_coeffs))
        product_stoich_coeffs = numpy.array(self.get_ordered_list(self.product_stoich_coeffs))
        concen_array = self.concentrations

        if (reactant_stoich_coeffs < 0).any():
            raise ValueError("Reactant stoichiometric coefficients must be positive!")
        
        if (product_stoich_coeffs < 0).any():
            raise ValueError("Product stoichiometric coefficients must be positive!")

        try:
            progress_rate = self.compute_progress_rate(T)
            nu_i = product_stoich_coeffs - reactant_stoich_coeffs

            reaction_rate_1_eq = progress_rate * nu_i
            return reaction_rate_1_eq

        except NotImplementedError:
            raise ElementaryReactionError('''You must first implement the functions to
                                compute the reaction rate coefficients and progress rates!''')


class IrrevElemReactionError(Exception):
    """Error for misclassified IrreversibleElementaryReaction."""
    pass

class IrrevElemReaction(ElementaryReaction):
    """Class for irreversible elementary reaction"""
    def __init__(self, rxn_type, is_reversible, rxn_equation, species_list, rate_coeffs_components,
                 reactant_stoich_coeffs, rate_coeffs_type, product_stoich_coeffs):
        """Initializes reaction that is irreversible and elementary.

        NOTES:
        ------
        PRE:
            - Raises IrrevElemReactionError if reaction type is not irreversible OR if reaction
                is not elementary (must satisfy both!)
        """
        super(IrrevElemReaction, self).__init__(rxn_type, is_reversible, rxn_equation,
                                                species_list, rate_coeffs_components,
                                                rate_coeffs_type,
                                                reactant_stoich_coeffs, product_stoich_coeffs)
        if not (rxn_type == "Elementary" and is_reversible == False):
            raise IrrevElemReactionError("This reaction is not irreversible nor elementary!") 

    def compute_reaction_rate_coeff(self, T=None):
        """Computes reaction rate coefficients of reaction.

        INPUTS:
        -------
        T : float
            temperature of reaction, in K

        RETURNS:
        --------
        k : numeric type (or list of numeric type)
            Reaction rate coefficient
        """
        T = self.temperature
        rxn_rate_coeff_obj = determine_rxn_rate_coeff_type(self.rate_coeffs_type)
        k = rxn_rate_coeff_obj(self.rate_coeffs_components,
                               T=self.temperature).k
        self.rxn_rate_coeff = k
        return k

    def compute_progress_rate(self, T=None):
        """Computes progress rates of reaction.

        INPUTS:
        -------
        T : float
            temperature of reaction, in K

        RETURNS:
        --------
        omega_array : numpy.ndarray
            Array of progress rates of reaction

        NOTES:
        ------
        PRE:
            - Raises ValueError if concentrations have not been set
        """
        T = self.temperature
        reactant_stoich_coeffs = numpy.array(self.get_ordered_list(self.reactant_stoich_coeffs))
        concen_array = self.concentrations

        if len(concen_array) == 0:
            raise ValueError("You must set the concentrations first!")

        k = self.compute_reaction_rate_coeff(T)

        concen_powered_j = concen_array ** reactant_stoich_coeffs

        progress_rate = k * numpy.prod(concen_powered_j)
        return progress_rate


class RevElemReactionError(Exception):
    """Error for misclassified ReversibleElementaryReaction."""
    pass

class RevElemReaction(ElementaryReaction):
    """Class for reversible reaction"""
    def __init__(self, rxn_type, is_reversible, rxn_equation, species_list, rate_coeffs_components, rate_coeffs_type,
                 reactant_stoich_coeffs, product_stoich_coeffs, bkwd_coeff_type="NASA7"):
        """Initializes reaction that is reversible and elementary.

        NOTES:
        ------
        PRE:
            - Raises RevElemReactionError if reaction type is not reversible OR if reaction
                is not elementary (must satisfy both!)
        """
        super(RevElemReaction, self).__init__(rxn_type, is_reversible, rxn_equation,
                                              species_list, rate_coeffs_components, rate_coeffs_type,
                                              reactant_stoich_coeffs, product_stoich_coeffs)

        self.bkwd_coeff_type = bkwd_coeff_type
        self.NASA_poly_coefs_dict = None
        self.NASA_poly_coefs = None

        if not (rxn_type == "Elementary" and is_reversible == True):
            raise RevElemReactionError("This reaction is not reversible nor elementary!") 

    def compute_reaction_rate_coeff(self, T=None):
        """Computes reaction rate coefficients of reaction.

        Args:
        =====
        T : float, optional
            Temperature of reaction, in Kelvin

        Returns:
        ========
        kf : numeric type (or list of numeric type)
            forward reaction rate coefficient
        kb : numeric type (or list of numeric type)
            backward reaction rate coefficient

        Notes:
        ======
        PRE:
            - Raises ValueError if NASA polynomial coefficients have not been
                set before trying to compute reaction rate coefficients
                (See class function set_NASA_poly_coefs())
        """
        rxn_rate_coeff_obj = determine_rxn_rate_coeff_type(self.rate_coeffs_type)
        coeffs = rxn_rate_coeff_obj(self.rate_coeffs_components, T=self.temperature)
        self.forward_rxn_rate_coeff = coeffs.k
        
        reactant_stoich_coeffs = numpy.array(self.get_ordered_list(self.reactant_stoich_coeffs))
        product_stoich_coeffs = numpy.array(self.get_ordered_list(self.product_stoich_coeffs))
        nui = product_stoich_coeffs - reactant_stoich_coeffs
        
        if (self.NASA_poly_coefs is None):
            raise ValueError("Must set NASA polynomial coefficients before computing rxn rate coefficients!")
        
        if self.bkwd_coeff_type == "NASA7":
            back_coeffs = NASA7BackwardCoeff(nui, self.NASA_poly_coefs)
        else:
            raise NotImplementedError("This type of coefficients for computing backward "
                                      "reaction rate coefficients is not currently suppoerted.")

        self.backward_rxn_rate_coeff = back_coeffs.compute_bkwd_coefficient(self.forward_rxn_rate_coeff,
                                                                           self.temperature)
        return self.forward_rxn_rate_coeff, self.backward_rxn_rate_coeff

    def set_NASA_poly_coefs(self, coefs):
        """Sets NASA polynomial coefficients.
        
        INPUTS:
        -------
        coefs : dict
            dictionary of NASA polynomial coefficients
            (For format, see Parser class' get_NASA_poly_coefs() function)
        """
        # order NASA polynomial coefficients by species (order specified by species_list)
        index_map = {v: i for i, v in enumerate(self.species_list)}
        sorted_tuple_list = sorted(coefs.items(), key=lambda pair: index_map[pair[0]])
        list_of_interest = [element[1] for element in sorted_tuple_list]
        
        # lists of high and low T polynomial coefficients IN ORDER of species
        list_of_interest = numpy.array(list_of_interest)

        self.NASA_poly_coefs_dict = coefs
        self.NASA_poly_coefs = list_of_interest

    def compute_progress_rate(self, T=None):
        """Computes progress rates of reaction.

        Args:
        =====
        T : float, optional
            Temperature of reaction, in Kelvin

        Returns:
        ========
        omega_array : numpy.ndarray
            Array of reaction progress rates

        Notes:
        ======
        PRE:
            - Raises ValueError if concentrations have not been not set
        """
        reactant_stoich_coeffs = numpy.array(self.get_ordered_list(self.reactant_stoich_coeffs))
        product_stoich_coeffs = numpy.array(self.get_ordered_list(self.product_stoich_coeffs))
        concen_array = self.concentrations

        if len(concen_array) == 0:
            raise ValueError("You must set the concentrations first!")

        # Compute forward (and backward, if applicable) rxn rate coefficients
        self.compute_reaction_rate_coeff(T)
  
        concen_powered_forward = concen_array ** reactant_stoich_coeffs
        concen_powered_backward = concen_array ** product_stoich_coeffs

        # Compute progress rates
        progress_rate_forward = self.forward_rxn_rate_coeff * numpy.prod(concen_powered_forward)
        progress_rate_backward = self.backward_rxn_rate_coeff * numpy.prod(concen_powered_backward)

        return progress_rate_forward - progress_rate_backward
