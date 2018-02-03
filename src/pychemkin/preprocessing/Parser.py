
"""Classes for preprocessing: parsing xml files, identifying reaction type."""

from enum import Enum
import numpy
import xml.etree.ElementTree as ET
from pychemkin.pychemkin_errors import PyChemKinError


class RxnType(Enum):
    Elementary = 1


class XmlParser():
    """Class for parsing input XML files to retrieve and
    preprocess reaction data.
    """
    def __init__(self, path):
        """Initializes XML file parser.

        Args:
        -----
        path : str
            file path of input XML file
        """

        if path[-4:] != '.xml':
            path += '.xml'
        self.path = path

    def load(self):
        """Parses XML file contents and loads into lists of species and
        RxnData objects representing reactions in the file.

        Returns:
        --------
        species : list[str]
            list of names of species involved in system of reactions
        list_RxnData : list[RxnData]
            list of RxnData objects, each representing a reaction
        """
        tree = ET.parse(self.path)
        root = tree.getroot()

        species = []
        for species_i in root.find('phase'):
            new_species = species_i.text.strip().split()
            species.extend(new_species)

        list_RxnData = []
        for rxn in root.find('reactionData').findall('reaction'):
            rxn_data = self.__extract_data_from_reaction_element(rxn)
            list_RxnData.append(rxn_data)
        return species, list_RxnData

    def __extract_data_from_reaction_element(self, rxn):
        """Returns RxnData object containing data from <reaction> XML element.

        Args:
        -----
        rxn : xml.etree.ElementTree.Element
            <reaction> XML element containing information about a reaction

        Returns:
        --------
        rxn_data : RxnData
            RxnData object containing information from input, rxn

        Notes:
        ------
            - Raises PyChemKinError for invalid attribute/element values
        """
        rxn_data = RxnData()
        rxn_data.rxn_id = rxn.get('id')

        # Get info about rxn reversibility
        reversible = rxn.get('reversible').lower().strip()
        if reversible in ['no', 'n', 'false', 'f']:
            rxn_data.is_reversible = False
        elif reversible in ['yes', 'y', 'true', 't']:
            rxn_data.is_reversible = True
        else:
            raise PyChemKinError(
                  'XmlParser.load()',
                  'Invalid reversibility attribute in reaction {}'.format(
                        result.rxn_id))

        # Get info about rxn type
        rnx_type = rxn.get('type').lower().strip()
        if rnx_type == 'elementary':
            rxn_data.type = RxnType.Elementary
        else:
            raise PyChemKinError(
                  'XmlParser.load()',
                  'Reaction {} is non-elementary.'.format(
                        rxn_data.rxn_id))

        # Get info about rxn rate coefficient
        rate_coeff = rxn.find('rateCoeff')
        if rate_coeff is None:
            raise PyChemKinError(
                    'XmlParser.load()',
                    'No <rateCoeff> element found in one of the '
                    'reactions.')

        if rate_coeff.find('Arrhenius') is not None:
            arrhenius = rate_coeff.find('Arrhenius')
            A = float(arrhenius.find('A').text.strip())
            if A < 0:
                raise PyChemKinError(
                        'XmlParser.load()',
                        'A coeff < 0 in reaction with '
                        'id = {}'.format(rxn_data.rxn_id))
            E = float(arrhenius.find('E').text.strip())
            rxn_data.rate_coeff = [A, E]

        elif rate_coeff.find('modifiedArrhenius') is not None:
            mod_arrhenius = rate_coeff.find('modifiedArrhenius')
            A = float(mod_arrhenius.find('A').text.strip())
            if A < 0:
                raise PyChemKinError(
                      'A coeff < 0 in reaction with id = {}'.format(
                            rxn_data.rxn_id))
            b = float(mod_arrhenius.find('b').text.strip())
            E = float(mod_arrhenius.find('E').text.strip())
            rxn_data.rate_coeff = [A, b, E]

        elif rate_coeff.find('Constant') is not None:
            const = rate_coeff.find('Constant')
            rxn_data.rate_coeff = float(const.find('k').text.strip())

        else:
            raise PyChemKinError(
                    'XmlParser.load()',
                    'No recognized child of <rateCoeff> found '
                    'from which to parse coefficients.')

        # Map reaction species to corresponding stoichiometric coefficients
        rxn_data.reactants = self.__map_stoichcoeff_to_species(
                                rxn.find('reactants'))
        rxn_data.products = self.__map_stoichcoeff_to_species(
                                rxn.find('products'))
        rxn_data.rxn_equation = rxn_data.get_equation()
        return rxn_data

    @staticmethod
    def __map_stoichcoeff_to_species(tag):
        """Creates dictionary mapping species to corresponding
        stoichiometric coefficients.

        Args:
        -----
        tag : xml.etree.ElementTree.Element
            <reactants> or <products> XML element

        Returns:
        --------
        species_info_dict : dict
            dictionary mapping species to stoichiometric coefficients
        """
        species_info_dict = {}
        for item in tag.text.strip().split():
            species, stoichcoeff = item.strip().split(':')
            stoichcoeff = int(stoichcoeff)
            species = species.upper()
            species_info_dict[species] = stoichcoeff
        return species_info_dict


# TODO: break this up using helper functions?
    def populate_parsed_data_list(self, Ti):
        """Returns a list of dictionaries, with each containing
        reaction parameters at a given temperature.

        Args:
        -----
        Ti : list[float]
            Reaction temperatures, in Kelvin

        Returns:
        --------
        parsed_data_dict_list: list[dict]
            List of dictionaries containing reaction parameters

        Notes:
        ------
            - Returns a list of dictionaries, each with following attributes:
                parsed_data_dict['species'] : a list of reaction species
                parsed_data_dict['ki'] : a list of reaction rate coefficients,
                    i-th item for i-th reaction
                parsed_data_dict['sys_vi_r'] : a list of stoichiometric
                    coefficients of the reactants
                parsed_data_dict['sys_vi_p'] : a list of stoichiometric
                    coefficients of the products
                parsed_data_dict['is_reversible'] : a boolean indicating
                    whether the reaction is reversible
                parsed_data_dict['T'] : a float indicating reaction temperature
                    in Kelvin
        """
        species, rxn_data_list = self.load()
        n_species = len(species)

        # Build dictionary mapping species name to index for ordering
        species_idx_dict = {}
        for i, s in enumerate(species):
            species_idx_dict[s] = i

        parsed_data_dict_list = []

        # Loop thru reaction temperatures
        for T in Ti:

            # Initialized for set of reactions in XML
            stoichcoeffs_reactants = []
            stoichcoeffs_products = []
            ki = []
            is_reversible_list = []
            equations = []

            for rxn_data in rxn_data_list:

                if rxn_data.type != RxnType.Elementary:
                    raise PyChemKinError(
                            'XmlParser.populate_parsed_data_list(Ti)',
                            'Non-elementary reactions cannot be '
                            'parsed now.')

                if rxn_data.is_reversible:
                    is_reversible_list.append(True)
                else:
                    is_reversible_list.append(False)

                rxn_id = rxn_data.rxn_id

                rxn_vi_r = numpy.zeros(n_species)
                for s, vi in rxn_data.reactants.items():
                    idx = species_idx_dict[s]
                    rxn_vi_r[idx] = vi
                stoichcoeffs_reactants.append(list(rxn_vi_r))

                rxn_vi_p = numpy.zeros(n_species)
                for s, vi in rxn_data.products.items():
                    idx = species_idx_dict[s]
                    rxn_vi_p[idx] = vi
                stoichcoeffs_products.append(list(rxn_vi_p))

                coef_params = rxn_data.rate_coeff
                if isinstance(coef_params, list):
                    # modified Arrhenius
                    if len(coef_params) == 3:
                        A = coef_params[0]
                        b = coef_params[1]
                        E = coef_params[2]
                        ki.append(
                              ModifiedArrheniusCoefficient(A, b, E,
                                                           T).get_coef())
                    # Arrhenius
                    else:
                        A = coef_params[0]
                        E = coef_params[1]
                        ki.append(ArrheniusCoefficient(A, E, T).get_coef())

                # Constant
                else:
                    ki.append(ConstantCoefficient(coef_params).get_coef())

                if rxn_data.rxn_equation is None:
                    rxn_data.rxn_equation = "Reaction equation not specified"
                equations.append(rxn_data.rxn_equation)

            # Populate dictionary for set of reactions
            parsed_data_dict = {}
            parsed_data_dict['equations'] = equations
            parsed_data_dict['species'] = species
            parsed_data_dict['ki'] = ki
            parsed_data_dict['sys_vi_r'] = sys_vi_r
            parsed_data_dict['sys_vi_p'] = sys_vi_p
            parsed_data_dict['is_reversible'] = is_reversible_list
            parsed_data_dict['T'] = T

            try:
                b_ki = BackwardCoefficient(
                        species, T, ki,
                        is_reversible_list,
                        sys_vi_r, sys_vi_p).get_backward_coefs()
                parsed_data_dict['b_ki'] = b_ki
            except PyChemKinError as err:
                parsed_data_dict['b_ki'] = 'Not Defined'

            parsed_data_dict_list.append(parsed_data_dict)

        return parsed_data_dict_list


class RxnData():
    """Class for storing reaction data.

    Attributes:
    -----------
    rxn_id : str
        id attribute of <reaction> XML element
    is_reversible : bool
        True if reaction is reversible
    reactants : dict[str: int]
        Mapping of reactant species to corresponding
        stoichiometric coefficients
    products : dict[str: int]
        Mapping of product species to corresponding
        stoichiometric coefficients
    rate_coeff : list[float] or float
        Reaction rate coefficients, represented by list or float depending
        on type of reaction rate coefficient
    rxn_equation : str
        Equation representation of reaction
    rxn_type : RxnType
        Enum value for reaction type
    """
    def __init__(self, rxn_id=None, is_reversible=None,
                 reactants=None, products=None, rate_coeff=None,
                 rxn_equation=None, rxn_type=None):
        self.rxn_id = rxn_id
        self.is_reversible = is_reversible
        self.reactants = reactants
        self.products = products
        self.rate_coeff = rate_coeff
        self.rxn_equation = rxn_equation
        self.rxn_type = rxn_type

    def get_equation(self):
        """Returns equation representation of reaction.

        Returns:
        --------
        full_rxn_eqn : str
            Equation representation of reaction, as a string

        Note:
        -----
            - Species are listed in alphabetical order on both
                reactant and product sides.
        """
        reactant_side = self.__build_equation_side(self.reactants)
        product_side = self.__build_equation_side(self.products)

        if self.is_reversible:
            full_rxn_eqn = '{0} [=] {1}'.format(reactant_side, product_side)
        else:
            full_rxn_eqn = '{0} =] {1}'.format(reactant_side, product_side)
        return full_rxn_eqn

    @staticmethod
    def __build_equation_side(species_dict):
        """Builds one side of reaction equation (reactant or product side).

        Args:
        -----
        species_dict : dict[str, int]
            dictionary mapping species to stoichiometric coefficients

        Returns:
        --------
        rxn_eqn : str
            One side of reaction equation as a string
        """
        rxn_eqn = ''
        species_count = 0
        for species, stoich_coeff in sorted(species_dict.items()):
            species_count += 1

            # Handle first species of each reaction equation
            # differently than rest, then proceed to next species.
            if species_count == 1:
                if stoich_coeff == 1:
                    stoich_coeff_part = ''
                else:
                    stoich_coeff_part = stoich_coeff
                rxn_eqn += '{0}{1}'.format(stoich_coeff_part, species)
                continue

            # Non-first species.
            operator = ' + '
            if stoich_coeff == 1:
                stoich_coeff_part = ''
            else:
                stoich_coeff_part = abs(stoich_coeff)
            rxn_eqn += '{0}{1}{2}'.format(operator, stoich_coeff_part, species)
        return rxn_eqn
