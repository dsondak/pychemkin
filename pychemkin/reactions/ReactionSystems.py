
"""Class for a system of reactions."""

import numpy

from pychemkin.reactions.Reactions import *

# from scipy.integrate import ode


class ReactionSystem:
    """Class for a system of reactions."""
    def __init__(self, reaction_list, NASA_poly_coefs, temperature, concentrations):
        """Initializes a reaction system.
        
        Args:
        =====
        reaction_list : list[Reaction] or list[Reaction-inherited]
            list of Reaction or Reaction-inherited objects
        NASA_poly_coefs :  

        temperatures : numeric type or array of numeric types
            temperatures of reaction
        concentrations : dict
            dictionary of concentrations (in molar) with key as species name

        Attributes:
        ===========
        involved_species : list[str]
            list of all species in reaction system
        """
        self.reaction_list = reaction_list

        if temperature <= 0:
            raise ValueError("Temperature has to be a positive value!")

        self.temperature = temperature
        self.NASA_matrix = self.get_nasa_matrix(NASA_poly_coefs)

        # NOTE: Not sure if this is good enough!
        self.involved_species = reaction_list[0].species_list
        self.vis_concentrations = concentrations

        # Set up each reaction
        for r in self.reaction_list:
            r.set_concentrations(concentrations)
            r.set_temperature(self.temperature)
            
            if isinstance(r, RevElemReaction):
                r.set_NASA_poly_coefs(self.NASA_matrix)


        # self.concentrations = reaction_list[0].concentrations
        # # ODE integrator possible choices of solver: dopri5, dop853(both are explicit ruggi-kutta method) vode, zvode
        # self.r = ode(self.compute_reaction_rate).set_integrator("dop853")
        # self.r.set_initial_value(self.concentrations, 0)


    def get_reaction_rate(self):
        """Returns reaction rate for each reaction.

        Returns:
        ========
        list_rxn_rates : list[float]
            List of reaction rates of reactions in the system
        """
        reaction_rate_list = [rxnObj.compute_reaction_rate() for rxnObj in self.reaction_list]
        reaction_rate_list = numpy.array(reaction_rate_list)
        rxnrates = numpy.sum(reaction_rate_list, axis=0)
        return rxnrates

    def sort_reaction_rates(self):
        """Sorts the reaction rates and returns them as
        a dictionary.

        Returns:
        ========
        rxn_rates_dict: Sorted dictionary of reaction rate
        """
        rxn_rates_dict = {}
        list_species_ordered = list(self.involved_species)
        rxnrate = self.get_reaction_rate()
        for i in range(len(rxnrate)):
            rxn_rates_dict[list_species_ordered[i]] = rxnrate[i]
        return rxn_rates_dict

    def get_nasa_matrix(self, NASA_poly_coeff):
        """Computes array of NASA polynomial coefficients.

        INPUTS:
        -------
        NASA_poly_coef : list[dict]
            list of dictionaries of NASA polynomial coefficients
                labeled by temperature range

        RETURNS:
        --------
        NASA_array : numpy.ndarray
            array of NASA polynomial coefficients for given temperature range
        """

        NASA = {}
        for specie in NASA_poly_coeff:
            specie_dict = NASA_poly_coeff[specie]
<<<<<<< HEAD
            if self.temperature <= numpy.float(specie_dict["Tmid"]): # get the low temperature
=======
            if self.temperature <= specie_dict["Tmid"]: # get the low temperature
>>>>>>> origin/develop
                NASA[specie] = specie_dict["low"]
            else:
                NASA[specie] = specie_dict["high"]
        return NASA


    # def compute_reaction_rate(self, a, concentrations):
    #     '''
    #     Ordinary differential equation for the ODE solver. It gets the current concentration,
    #     modifies the concentration in each reaction and computes the current reaction rate(first-order derivative)
    #     INPUTS:
    #     -------
    #     a : float
    #         placeholder of time to match the format of the ode solver.
    #         In our ODE the derivative doesn't depend on the current time, but only on the current state(concentration)
    #     concentration: list[float]
    #         current state (current concentration)

    #     RETURNS:
    #     --------
    #     self.get_reaction_rate() : tuple
    #         current reaction rate
    #     '''
    #     #update the concentration of the species in each reaction and also in the system
    #     for r in self.reaction_list:
    #         r.set_concentrations_from_array(concentrations)
    #     self.concentrations = concentrations
    #     return self.get_reaction_rate()

    # def step(self, dt):
    #     """Solve the ODE， get the state after dt time

    #     INPUTS:
    #     -------
    #     dt : float
    #         timestep of the next state

    #     RETURNS:
    #     --------
    #     (self.r.t, self.r.y) : tuple
    #         current time and current concentration"""
    #     self.r.integrate(self.r.t + dt)
    #     return self.r.t, self.r.y
