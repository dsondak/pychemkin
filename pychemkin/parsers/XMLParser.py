
"""Classes for preprocessing: parsing xml files, identifying reaction type."""

from enum import Enum
import numpy
import os
import xml.etree.ElementTree as ET
from pychemkin.reactions.Reactions import *


class XMLParser:
    """Parser for input XML files to retrieve and
    preprocess reaction data."""
    def __init__(self, xml_filename):
        """Initializes parser for input XML files.
        
        Args:
        =====
        xml_filename : str, required
            name of input XML file
        
        Atrributes:
        ===========
        rxns : xml.etree.ElementTree.Element
            root node of xml tree
        reaction_list : list[Reaction]
            list of Reaction (or Reaction-inherited) objects
        species : dict[str]
            dictionary containing names of species

        Notes:
        ======
        POST:
            - Raises IOError if XML file not found
        """
        if os.path.isfile(xml_filename):
            self.xml_filename = xml_filename
            tree = ET.parse(self.xml_filename)
            self.rxns = tree.getroot()
        else:
            raise IOError("Reaction (xml) file not found!")

        self.reaction_list = []
        self.species = {}
        self.get_species()
        self.get_reaction_list()

    def get_species(self):
        """Populates and returns dictionary of species from
        species data in XML file.
        
        Returns:
        ========
        species: dict[str]  
            dictionary of the form, key: name of species, value: None
        """
        species_list = []
        for specie_i in self.rxns.find('phase'):
            new_specie = specie_i.text.strip().split()
            species_list.extend(new_specie)
            
        for specie in species_list:
            self.species[specie] = None
        return self.species

    def get_rxn_type(self, reaction):
        """Helper function that returns reaction type from input.
        
        Args:
        =====
        reaction : xml.etree.ElementTree.Element, required
            <reaction> XML element containing information about a reaction
        
        Returns:
        ========
        rxn_type : str
            string describing reaction type (e.g. "elementary")
        """
        rxn_type = reaction.get('type')
        return rxn_type

    def get_is_reversible(self, reaction):
        """Helper function that returns information about whether the reaction is reversible.
        
        Args:
        =====
        reaction : xml.etree.ElementTree.Element, required
            <reaction> XML element containing information about a reaction
        
        RETURNS:
        --------
        is_reversible : bool
            True if reversible, False if irreversible
        """
        return (reaction.get('reversible') == "yes")

    def get_rxn_equation(self, reaction):
        """Helper function that returns reaction equation.
        
        Args:
        =====
        reaction : xml.etree.ElementTree.Element, required
            <reaction> XML element containing information about a reaction
        
        RETURNS:
        --------
        rxn_equation : str
            a string representation of reaction equation
        """
        rxn_equation = reaction.find('equation').text
        return rxn_equation

    def get_rate_coeffs_components(self, reaction):
        """Helper function that returns reaction rate coefficient components
        based on type of coefficient.
        
        Args:
        =====
        reaction : xml.etree.ElementTree.Element, required
            <reaction> XML element containing information about a reaction
        
        RETURNS:
        --------
        rate_coeffs_components : dict
            dictionary of the form {coefficient component name: coefficient component value}. 
        """
        for coeff in reaction.findall('rateCoeff'):

            # constant-type
            if coeff.find('Constant') is not None:
                for arr in coeff.findall('Constant'):
                    if arr.find('k') is None or arr.find('k').text is None:
                        raise ValueError("Missing component 'k' for constant type rxn rate coefficient.")
                    else:
                        k = float(arr.find('k').text)
                rate_coeffs_components = {"k": k}

            # Arrhenius-type
            elif coeff.find('Arrhenius') is not None:
                for arr in coeff.findall('Arrhenius'):
                    if arr.find('A') is None or arr.find('A').text is None:
                        raise ValueError("Missing component 'A' for modified Arrhenius type rxn rate coefficient.")
                    else:
                        A = float(arr.find('A').text)
                    if arr.find('E') is None or arr.find('E').text is None:
                        raise ValueError("Missing component 'E' for modified Arrhenius type rxn rate coefficient.")
                    else:
                        E = float(arr.find('E').text)
                rate_coeffs_components = {"A": A, "E": E}

            # modified Arrhenius-type
            elif coeff.find('modifiedArrhenius') is not None:
                for arr in coeff.findall('modifiedArrhenius'):
                    if arr.find('A') is None or arr.find('A').text is None:
                        raise ValueError("Missing component 'A' for modified modified Arrhenius type rxn rate coefficient.")
                    else:
                        A = float(arr.find('A').text)
                    if arr.find('b') is None or arr.find('b').text is None:
                        raise ValueError("Missing component 'b' for modified modified Arrhenius type rxn rate coefficient.")
                    else:
                        b = float(arr.find('b').text)
                    if arr.find('E') is None or arr.find('E').text is None:
                        raise ValueError("Missing component 'E' for modified modified Arrhenius type rxn rate coefficient.")
                    else:
                        E = float(arr.find('E').text)
                rate_coeffs_components = {"A": A, "b": b, "E": E}

            return rate_coeffs_components

    def get_reactant_stoich_coeffs(self, reaction):
        """Helper function that returns reactant stoichiometric coefficients from input reaction.
        
        Args:
        =====
        reaction : xml.etree.ElementTree.Element, required
            <reaction> XML element containing information about a reaction
        
        Returns:
        ========
        reactant_stoich_coeffs : dict[str, int]
            dictionary in the form {reactant name: stoich coefficient}
        """
        reactant_stoich_coeffs = {}
        for reactant in reaction.find('reactants').text.split():
            name = reactant.split(":")[0]
            stoich_coeff = reactant.split(":")[1]
            reactant_stoich_coeffs[name] = int(stoich_coeff)
        return reactant_stoich_coeffs
    
    def get_product_stoich_coeffs(self, reaction):
        """Helper function that returns product stoichiometric coefficients from input.
        
        Args:
        =====
        reaction : xml.etree.ElementTree.Element, required
            <reaction> XML element containing information about a reaction
        
        Returns:
        ========
        product_stoich_coeffs : dict[str, int]
            dictionary in the form {product name: stoich coefficient}
        """
        product_stoich_coeffs = {}
        for product in reaction.find('products').text.split():
            name = product.split(":")[0]
            stoich_coeff = product.split(":")[1]
            product_stoich_coeffs[name] = int(stoich_coeff)
        return product_stoich_coeffs

    def get_reaction_list(self):
        """Loops thru reactions in input XML file to populate
        a list of Reaction or Reaction-inherited objects containing
        information about corresponding reaction."""
        for reactionData in self.rxns.findall('reactionData'):
            for reaction in reactionData.findall('reaction'):
                species = self.get_species()
                is_reversible = self.get_is_reversible(reaction) 
                rxn_type = self.get_rxn_type(reaction)
                rxn_equation = self.get_rxn_equation(reaction)
                rate_coeffs_components = self.get_rate_coeffs_components(reaction)
                reactant_stoich_coeffs = self.get_reactant_stoich_coeffs(reaction)
                product_stoich_coeffs = self.get_product_stoich_coeffs(reaction)
                

### TO FIX!!!!!!!!!!!!!!!!

                # IRREVERSIBLE elementary reaction case
                if is_reversible == False and rxn_type == "Elementary":
                    rxn = IrreversibleElementaryReaction(rxn_type, is_reversible, rxn_equation,
                                           species, rate_coeffs_components,
                                           reactant_stoich_coeffs, product_stoich_coeffs)
                    self.reaction_list.append(rxn)

                # REVERSIBLE elementary reaction case
                elif is_reversible == True and rxn_type == "Elementary":
                    rxn = ReversibleElementaryReaction(rxn_type, is_reversible, rxn_equation,
                                         species, rate_coeffs_components,
                                         reactant_stoich_coeffs, product_stoich_coeffs)
                    self.reaction_list.append(rxn)

                # Unhandled reaction case
                else:
                    raise NotImplementedError("This type of reaction has not been implemented yet!")












# class XMLParser:
#     """Parser for input XML files to retrieve and
#     preprocess reaction data.
#     """
#     def __init__(self, filename):
#         """Initializes XML file parser.

#         Args:
#         =====
#         filename : str, required
#             name of input XML file
#         """
#         if filename[-4:] != '.xml':
#             filename += '.xml'
#         self.filename = filename

#     def load(self):
#         """Parses XML file contents and loads into lists of species and
#         RxnData objects representing reactions in the file.

#         Returns:
#         ========
#         species : list[str]
#             list of names of species involved in system of reactions
#         list_RxnData : list[RxnData]
#             list of RxnData objects, each representing a reaction
#         """
#         tree = ET.parse(self.filename)
#         root = tree.getroot()

#         species = []
#         for species_i in root.find('phase'):
#             new_species = species_i.text.strip().split()
#             species.extend(new_species)

#         list_RxnData = []
#         for rxn in root.find('reactionData').findall('reaction'):
#             rxn_data = self.__extract_data_from_reaction_element(rxn)
#             list_RxnData.append(rxn_data)
#         return species, list_RxnData


# ### TODO: break this up into smaller functions?
#     def __extract_data_from_reaction_element(self, rxn):
#         """Returns RxnData object containing data from <reaction> XML element.

#         Args:
#         =====
#         rxn : xml.etree.ElementTree.Element, required
#             <reaction> XML element containing information about a reaction

#         Returns:
#         ========
#         rxn_data : RxnData
#             RxnData object containing information from input, rxn

#         Notes:
#         ======
#             - Raises PyChemKinError for invalid attribute/element values
#         """
#         rxn_data = RxnData()
#         rxn_data.rxn_id = rxn.get('id')

#         # Get info about rxn reversibility
#         reversible = rxn.get('reversible').lower().strip()
#         if reversible in ['no', 'n', 'false', 'f']:
#             rxn_data.is_reversible = False
#         elif reversible in ['yes', 'y', 'true', 't']:
#             rxn_data.is_reversible = True
#         else:
#             raise PyChemKinError(
#                   'XMLParser.load()',
#                   'Invalid reversibility attribute in reaction {}'.format(
#                         result.rxn_id))

#         # Get info about rxn type
#         rnx_type = rxn.get('type').lower().strip()
#         if rnx_type == 'elementary':
#             rxn_data.type = RxnType.Elementary
#         else:
#             raise PyChemKinError(
#                   'XMLParser.load()',
#                   'Reaction {} is non-elementary.'.format(
#                         rxn_data.rxn_id))

#         # Get info about rxn rate coefficient
#         rate_coeff = rxn.find('rateCoeff')
#         if rate_coeff is None:
#             raise PyChemKinError(
#                     'XMLParser.load()',
#                     'No <rateCoeff> element found in one of the '
#                     'reactions.')

#         if rate_coeff.find('Arrhenius') is not None:
#             arrhenius = rate_coeff.find('Arrhenius')
#             A = float(arrhenius.find('A').text.strip())
#             if A < 0:
#                 raise PyChemKinError(
#                         'XMLParser.load()',
#                         'A coeff < 0 in reaction with '
#                         'id = {}'.format(rxn_data.rxn_id))
#             E = float(arrhenius.find('E').text.strip())
#             rxn_data.rate_coeff = [A, E]

#         elif rate_coeff.find('modifiedArrhenius') is not None:
#             mod_arrhenius = rate_coeff.find('modifiedArrhenius')
#             A = float(mod_arrhenius.find('A').text.strip())
#             if A < 0:
#                 raise PyChemKinError(
#                       'A coeff < 0 in reaction with id = {}'.format(
#                             rxn_data.rxn_id))
#             b = float(mod_arrhenius.find('b').text.strip())
#             E = float(mod_arrhenius.find('E').text.strip())
#             rxn_data.rate_coeff = [A, b, E]

#         elif rate_coeff.find('Constant') is not None:
#             const = rate_coeff.find('Constant')
#             rxn_data.rate_coeff = float(const.find('k').text.strip())

#         else:
#             raise PyChemKinError(
#                     'XMLParser.load()',
#                     'No recognized child of <rateCoeff> found '
#                     'from which to parse coefficients.')

#         # Map reaction species to corresponding stoichiometric coefficients
#         rxn_data.reactants = self.__map_stoichcoeff_to_species(
#                                 rxn.find('reactants'))
#         rxn_data.products = self.__map_stoichcoeff_to_species(
#                                 rxn.find('products'))
#         rxn_data.rxn_equation = rxn_data.get_equation()
#         return rxn_data

#     @staticmethod
#     def __map_stoichcoeff_to_species(tag):
#         """Creates dictionary mapping species to corresponding
#         stoichiometric coefficients.

#         Args:
#         =====
#         tag : xml.etree.ElementTree.Element, required
#             <reactants> or <products> XML element

#         Returns:
#         ========
#         species_info_dict : dict
#             dictionary mapping species to stoichiometric coefficients
#         """
#         species_info_dict = {}
#         for item in tag.text.strip().split():
#             species, stoichcoeff = item.strip().split(':')
#             stoichcoeff = int(stoichcoeff)
#             species = species.upper()
#             species_info_dict[species] = stoichcoeff
#         return species_info_dict


# # TODO: break this up using helper functions?
#     def populate_parsed_data_list(self, Ti):
#         """Returns a list of dictionaries, with each containing
#         reaction parameters at a given temperature.

#         Args:
#         =====
#         Ti : list[float], required
#             Reaction temperatures, in Kelvin

#         Returns:
#         ========
#         parsed_data_dict_list: list[dict]
#             List of dictionaries containing reaction parameters

#         Notes:
#         ======
#             - Returns a list of dictionaries, each with following attributes:
#                 parsed_data_dict['species'] : a list of reaction species
#                 parsed_data_dict['ki'] : a list of reaction rate coefficients,
#                     i-th item for i-th reaction
#                 parsed_data_dict['sys_vi_r'] : a list of stoichiometric
#                     coefficients of the reactants
#                 parsed_data_dict['sys_vi_p'] : a list of stoichiometric
#                     coefficients of the products
#                 parsed_data_dict['is_reversible'] : a boolean indicating
#                     whether the reaction is reversible
#                 parsed_data_dict['T'] : a float indicating reaction temperature
#                     in Kelvin
#         """
#         species, rxn_data_list = self.load()
#         n_species = len(species)

#         # Build dictionary mapping species name to index for ordering
#         species_idx_dict = {}
#         for i, s in enumerate(species):
#             species_idx_dict[s] = i

#         parsed_data_dict_list = []

#         # Loop thru reaction temperatures
#         for T in Ti:

#             # Initialized for set of reactions in XML
#             stoichcoeffs_reactants = []
#             stoichcoeffs_products = []
#             ki = []
#             is_reversible_list = []
#             equations = []

#             for rxn_data in rxn_data_list:

#                 if rxn_data.type != RxnType.Elementary:
#                     raise PyChemKinError(
#                             'XMLParser.populate_parsed_data_list(Ti)',
#                             'Non-elementary reactions cannot be '
#                             'parsed now.')

#                 if rxn_data.is_reversible:
#                     is_reversible_list.append(True)
#                 else:
#                     is_reversible_list.append(False)

#                 rxn_id = rxn_data.rxn_id

#                 rxn_vi_r = numpy.zeros(n_species)
#                 for s, vi in rxn_data.reactants.items():
#                     idx = species_idx_dict[s]
#                     rxn_vi_r[idx] = vi
#                 stoichcoeffs_reactants.append(list(rxn_vi_r))

#                 rxn_vi_p = numpy.zeros(n_species)
#                 for s, vi in rxn_data.products.items():
#                     idx = species_idx_dict[s]
#                     rxn_vi_p[idx] = vi
#                 stoichcoeffs_products.append(list(rxn_vi_p))

#                 coef_params = rxn_data.rate_coeff
#                 if isinstance(coef_params, list):
#                     # modified Arrhenius
#                     if len(coef_params) == 3:
#                         A = coef_params[0]
#                         b = coef_params[1]
#                         E = coef_params[2]
#                         ki.append(
#                               ModifiedArrheniusCoefficient(A, b, E,
#                                                            T).get_coef())
#                     # Arrhenius
#                     else:
#                         A = coef_params[0]
#                         E = coef_params[1]
#                         ki.append(ArrheniusCoefficient(A, E, T).get_coef())

#                 # Constant
#                 else:
#                     ki.append(ConstantCoefficient(coef_params).get_coef())

#                 if rxn_data.rxn_equation is None:
#                     rxn_data.rxn_equation = "Reaction equation not specified"
#                 equations.append(rxn_data.rxn_equation)

#             # Populate dictionary for set of reactions
#             parsed_data_dict = {}
#             parsed_data_dict['equations'] = equations
#             parsed_data_dict['species'] = species
#             parsed_data_dict['ki'] = ki
#             parsed_data_dict['sys_vi_r'] = sys_vi_r
#             parsed_data_dict['sys_vi_p'] = sys_vi_p
#             parsed_data_dict['is_reversible'] = is_reversible_list
#             parsed_data_dict['T'] = T

#             try:
#                 b_ki = BackwardCoefficient(
#                         species, T, ki,
#                         is_reversible_list,
#                         sys_vi_r, sys_vi_p).get_backward_coefs()
#                 parsed_data_dict['b_ki'] = b_ki
#             except PyChemKinError as err:
#                 parsed_data_dict['b_ki'] = 'Not Defined'

#             parsed_data_dict_list.append(parsed_data_dict)

#         return parsed_data_dict_list


# class RxnData:
#     """Class for storing data for a single reaction.

#     Attributes:
#     ===========
#     rxn_id : str
#         id attribute of <reaction> XML element
#     is_reversible : bool
#         True if reaction is reversible
#     reactants : dict[str: int]
#         Mapping of reactant species to corresponding
#         stoichiometric coefficients
#     products : dict[str: int]
#         Mapping of product species to corresponding
#         stoichiometric coefficients
#     rate_coeff : list[float] or float
#         Reaction rate coefficients, represented by list or float depending
#         on type of reaction rate coefficient
#     rxn_equation : str
#         Equation representation of reaction
#     rxn_type : RxnType
#         Enum value for reaction type
#     """
#     def __init__(self, rxn_id=None, is_reversible=None,
#                  reactants=None, products=None, rate_coeff=None,
#                  rxn_equation=None, rxn_type=None):
#         self.rxn_id = rxn_id
#         self.is_reversible = is_reversible
#         self.reactants = reactants
#         self.products = products
#         self.rate_coeff = rate_coeff
#         self.rxn_equation = rxn_equation
#         self.rxn_type = rxn_type

#     def get_equation(self):
#         """Returns equation representation of reaction.

#         Returns:
#         ========
#         full_rxn_eqn : str
#             Equation representation of reaction, as a string

#         Notes:
#         ======
#             - Species are listed in alphabetical order on both
#                 reactant and product sides.
#         """
#         reactant_side = self.__build_equation_side(self.reactants)
#         product_side = self.__build_equation_side(self.products)

#         if self.is_reversible:
#             full_rxn_eqn = '{0} [=] {1}'.format(reactant_side, product_side)
#         else:
#             full_rxn_eqn = '{0} =] {1}'.format(reactant_side, product_side)
#         return full_rxn_eqn

#     @staticmethod
#     def __build_equation_side(species_dict):
#         """Builds one side of reaction equation (reactant or product side).
#         Helper function for function get_equation().

#         Args:
#         =====
#         species_dict : dict[str, int], required
#             dictionary mapping species to stoichiometric coefficients

#         Returns:
#         ========
#         rxn_eqn : str
#             One side of reaction equation as a string
#         """
#         rxn_eqn = ''
#         species_count = 0
#         for species, stoich_coeff in sorted(species_dict.items()):
#             species_count += 1

#             # Handle first species of each reaction equation
#             # differently than rest, then proceed to next species.
#             if species_count == 1:
#                 if stoich_coeff == 1:
#                     stoich_coeff_part = ''
#                 else:
#                     stoich_coeff_part = stoich_coeff
#                 rxn_eqn += '{0}{1}'.format(stoich_coeff_part, species)
#                 continue

#             # Non-first species.
#             operator = ' + '
#             if stoich_coeff == 1:
#                 stoich_coeff_part = ''
#             else:
#                 stoich_coeff_part = abs(stoich_coeff)
#             rxn_eqn += '{0}{1}{2}'.format(operator, stoich_coeff_part, species)
#         return rxn_eqn
