
"""Base class for a reaction."""


class BaseReaction:
    """Base class for a reaction."""
    def __init__(self, rxn_type, is_reversible, rxn_equation, species_list,
                 rate_coeffs_components, reactant_stoich_coeffs, product_stoich_coeffs):
    """Initializes base reaction.

    Args:
    -----
    rxn_type : str
        type of reaction (e.g. "Elementary")
    is_reversible : bool
        True if reaction is reversible
    rate_coeffs_components : dict
        dictionary of components (e.g. 'A', 'b', and/or 'E')
        to compute reaction rate coefficients
    rxn_equation : str
        string representation of reaction equation
    reactant_stoich_coeffs : dict
        dictionary of integers for reactant stoichiometric coefficients
    product_stoich_coeffs : dict
        dictionary of integers for product stoichiometric coefficients
    unique_species : list
        list of unique chemical species from original xml file (could be a bunch
        of reactions sometimes containing common species)
    species_list : list
        list of chemical species from original xml file (useful for ordering)
    
    Attributes:
    -----------
    temperature : int or float
        temperature of reaction, in Kelvin
    concentrations : list
        concentrations of species involved in reaction
    rxn_rate_coeff : float
        reaction rate coefficient

    progress_rate :
    reaction_rate : 

    """
    # FINISH!!!!!!! 


# TODO: FINISH!
    def __str__(self):
        """Returns user-friendly representation of reaction.
        
        Returns:
        --------
        rxn_info : str
            string representation of reaction (reaction equation)
        """
        rxn_info = ''
        #"Reaction Information: \n{0}\n{1}\n".format(self.rxn_equation, )
        return rxn_info

# TODO: FINISH!
    def __repr__ (self):

        rxn_repr = 'ReactionBase(rxn_type)'
        return rxn_repr

    def __len__(self):
        """Returns number of chemical species in reaction.
        
        Returns:
        --------
        n_species : int
            Number of chemical species involved in the reaction
        """
        n_species = len(self.reactant_stoich_coeffs)
        return n_species

    def compute_progress_rate(self):
        """Computes progress rates of reaction.
        
        Returns:
        --------
        omega_array : numpy.ndarray
            Array of progress rates of reaction

        Notes:
        ------
            - Raises NotImplementedError.
        """
        raise NotImplementedError('Subclass must implement this method!')
        
    def compute_reaction_rate(self):
        """Computes reaction rate of reaction.
        
        Returns:
        --------
        rxn_rate_array: numpy.ndarray
            Array of reaction rates of reaction
        
        Notes:
        ------
            - Raises NotImplementedError.
        """
        raise NotImplementedError('Subclass must implement this method!')
