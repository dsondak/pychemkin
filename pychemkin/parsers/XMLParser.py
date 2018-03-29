
"""Class for parsing XML files containing reaction data."""

import csv
import numpy
import os
import xml.etree.ElementTree as ET
from pychemkin.config import UNITS_DIRECTORY
from pychemkin.reactions.Reactions import *


class XMLParser:
    """Parser for XML files containing reaction data."""
    def __init__(self, xml_filename, units_file=UNITS_DIRECTORY, convert_to_SI_units=False):
        """Initializes XML file parser.
        
        Args:
        =====
        xml_filename : str, required
            Name of input XML file
        units_file : str, optional (default: UNITS_DIRECTORY)
            File path to csv file containing SI unit conversion factors
        convert_to_SI_units : bool, optional (default: False)
            If True and XML file contains units, will convert them to SI units
        
        Atrributes:
        ===========
        rxns_node : xml.etree.ElementTree.Element
            Root node of xml tree
        reaction_list : list[Reaction]
            List of Reaction (or Reaction-inherited) objects
        species : dict[str]
            Dictionary containing names of species

        Notes:
        ======
        POST:
            - Raises IOError if XML file not found
        """
        if os.path.isfile(xml_filename):
            self.xml_filename = xml_filename
            tree = ET.parse(self.xml_filename)
            self.rxns_node = tree.getroot()
        else:
            raise OSError("Reaction (xml) file not found!")
        self.units_file = units_file
        self.reaction_list = []
        self.species = {}
        self.convert_to_SI_units = convert_to_SI_units
        self.populate_reaction_list()

    def get_species(self):
        """Populates and returns dictionary of
        reaction species in XML file.
        
        Returns:
        ========
        species: dict[str, None]
            Dictionary with name of species as key
            and None as value
        """
        species_list = []
        for specie_i in self.rxns_node.find('phase'):
            new_specie = specie_i.text.strip().split()
            species_list.extend(new_specie)
        for specie in species_list:
            self.species[specie] = None
        return self.species

    def get_rxn_type(self, reaction):
        """Helper function that returns reaction type
        of input reaction XML element.
        
        Args:
        =====
        reaction : xml.etree.ElementTree.Element, required
            <reaction> XML element containing information about a reaction
        
        Returns:
        ========
        rxn_type : str
            String description of reaction type (e.g. "elementary")
        """
        return reaction.get('type')

    def get_is_reversible(self, reaction):
        """Helper function that returns True if
        reaction is reversible.
        
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
        """Helper function that returns
        reaction equation of input reaction XML element.
        
        Args:
        =====
        reaction : xml.etree.ElementTree.Element, required
            <reaction> XML element containing information about a reaction
        
        RETURNS:
        --------
        rxn_equation : str
            a string representation of reaction equation
        """
        return reaction.find('equation').text

    def access_units(self):
        """Helper function that accesses the units.csv
        file for SI unit conversion constants.

        Returns:
        ========
        unit_conversion : dict
            Dictionary containing unit as key
            and conversion constant as value


        """
        with open(self.units_file, 'r') as unit:
            next(unit)
            unit_dict = dict(csv.reader(unit))
        unit_conversion = dict((unit, float(conversion))
                                    for unit, conversion in unit_dict.items())
        return unit_conversion

    def get_arrhenius_based_components(self, rate_coeff, unit_list=None):
        """Returns reaction rate coefficients of the
        type Arrhenius-based (e.g. Arrhenius, modified Arrhenius)

        Args:
        =====
        rate_coeff : xml.etree.ElementTree.Element
            <rateCoeff> XML element containing info about
            rate coefficient of a particular reaction

        Returns:
        ========
        rxn_rate_coeffs_components : dict
            Components for Arrhenius-type reaction
            rate coefficient
        """
        try:
            A_has_units = 'units' in rate_coeff.find('A').attrib
            E_has_units = 'units' in rate_coeff.find('E').attrib
        except AttributeError:
            raise ValueError("A or E not found in the XML file.")

        # if both A and E have units in XML file
        if (A_has_units and E_has_units):

            if self.convert_to_SI_units:

                A_unit = rate_coeff.find('A').attrib['units'].split('/')
                A_conversion_constants = []
                
                for unit in A_unit:
                    try:
                        A_conversion_constants.append(unit_list[unit])
                    except:
                        raise NotImplementedError("{0} not implemented for A".format(unit))

                A_conversion = numpy.prod(numpy.array(A_conversion_constants))

                E_unit = rate_coeff.find('E').attrib['units'].split('/')
                E_conversion_constants = []
                for unit in E_unit:
                    try:
                        E_conversion_constants.append(unit_list[unit])
                    except:
                        raise NotImplementedError("{0} not implemented for E".format(unit))
                E_conversion = numpy.prod(numpy.array(E_conversion_constants))

            else:
                A_conversion = 1.
                E_conversion = 1.

        # if only one of them has units
        elif (A_has_units or E_has_units):
            raise ValueError("A or E has units but not both. Please fix accordingly in the XML file!")

        # if no units provided
        else:
            if self.convert_to_SI_units:
                raise ValueError("Cannot convert to SI units. No units provided in XML file.")

            else:
                A_conversion = 1.
                E_conversion = 1.

        try:
            A = float(rate_coeff.find('A').text) * A_conversion
            E = float(rate_coeff.find('E').text) * E_conversion
        except:
            raise ValueError("Conversion failed. "
                             "Resulting units have not been converted.")

        if rate_coeff.find('R') is not None:
            try:
                R = float(rate_coeff.find('R').text)
            except:
                raise ValueError("R failed to be added. Please check your XML file.")
            rxn_rate_coeffs_components = {'A': A, 'E': E, 'R': R}
        else:
            rxn_rate_coeffs_components = {'A': A, 'E': E}
        return rxn_rate_coeffs_components


    def get_arrhenius_components(self, rate_coeff, unit_list=None):
        """Returns reaction rate coefficients of the
        Arrhenius type.

        Args:
        =====
        rate_coeff : xml.etree.ElementTree.Element
            <rateCoeff> XML element containing info about
            rate coefficient of a particular reaction

        Returns:
        ========
        rxn_rate_coeffs_components : dict
            Components for Arrhenius-type reaction
            rate coefficient
        """
        rxn_rate_coeffs_components = self.get_arrhenius_based_components(rate_coeff=rate_coeff,
                                                                         unit_list=unit_list)

        if rate_coeff.find('b') is not None:
            raise ValueError("Cannot use 'b' in Arrhenius type.")
        return rxn_rate_coeffs_components

    def get_mod_arrhenius_components(self, rate_coeff, unit_list=None):
        """Returns reaction rate coefficients of the
        Modified Arrhenius type.

        Args:
        =====
        rate_coeff : xml.etree.ElementTree.Element
            <rateCoeff> XML element containing info about
            rate coefficient of a particular reaction

        Returns:
        ========
        rxn_rate_coeffs_components : dict
            Components for Arrhenius-type reaction
            rate coefficient
        """
        rxn_rate_coeffs_components = self.get_arrhenius_based_components(rate_coeff=rate_coeff,
                                                                         unit_list=unit_list)
        try:
            b = float(rate_coeff.find('b').text)
            rxn_rate_coeffs_components['b'] = b
        except:
            raise ValueError("Conversion failed. "
                             "Resulting units have not been converted.")

        if rate_coeff.find('R'):
            try:
                R = float(rate_coeff.find('R').text)
            except:
                raise ValueError("R failed to be added. Please check your XML file.")
            rxn_rate_coeffs_components['R'] = R
        return rxn_rate_coeffs_components


    def get_rate_coeffs_components(self, reaction):
        """Helper function that returns reaction rate coefficient components
        based on type of coefficient.

        Args:
        =====
        reaction : xml.etree.ElementTree.Element, required
            <reaction> XML element containing information about a reaction
        convert_to_SI_units : boolean, optional (default False)
            converts units if True


        RETURNS:
        --------
        rxn_rate_coeffs_components : dict
            dictionary of the form {coefficient component name: coefficient component value}. 
        """
        unit_conversion = self.access_units()

        rateCoeffs = reaction.find('rateCoeff')
        for rateCoeff in rateCoeffs:

            if rateCoeff.tag == 'Arrhenius':
                rxn_rate_coeffs_components = self.get_arrhenius_components(rateCoeff, unit_list=unit_conversion)
                rxn_rate_coeffs_type = 'arrhenius'
            elif rateCoeff.tag in ('modifiedArrhenius', 'Kooij'):
                rxn_rate_coeffs_components = self.get_mod_arrhenius_components(rateCoeff, unit_list=unit_conversion)
                rxn_rate_coeffs_type = 'modified arrhenius'
            elif rateCoeff.tag == 'Constant':
                try:
                    k = float(rateCoeff.find('k').text)
                    rxn_rate_coeffs_components = {'k': k}
                except:
                    raise ValueError("Constant rxn rate coefficient, k failed to be added. "
                                     "Please check your XML file.")
                rxn_rate_coeffs_type = 'constant'

            else:
                raise NotImplementedError("{0} is not implemented.".format(rateCoeff.tag))

        return rxn_rate_coeffs_components, rxn_rate_coeffs_type

    def get_stoich_coefficients(self, species_type, reaction):
        """Helper function that returns stoichiometric coefficients
        of a given species type ('reactants' or 'products').

        Args:
        =====
        species_type : str, required
            Type of reaction species ('reactants' or 'products')
        reaction : xml.etree.ElementTree.Element, required
            <reaction> XML element containing information about a reaction

        Returns:
        ========
        stoich_coeffs : dict[str, int]
            Dictionary in the form {species name: stoich coefficient}

        Notes:
        ======
            - Used by get_reactant_stoich_coeffs and
                get_product_stoich_coeffs
        """
        stoich_coeffs = {}
        for specie in reaction.find(species_type).text.split():
            name = specie.split(":")[0]
            stoich_coeff = specie.split(":")[1]
            stoich_coeffs[name] = int(stoich_coeff)
        return stoich_coeffs


    def get_reactant_stoich_coeffs(self, reaction):
        """Helper function that returns reactant stoichiometric coefficients from input reaction.
        
        Args:
        =====
        reaction : xml.etree.ElementTree.Element, required
            <reaction> XML element containing information about a reaction
        
        Returns:
        ========
        reactant_stoich_coeffs : dict[str, int]
            Dictionary in the form {reactant name: stoich coefficient}
        """
        return self.get_stoich_coefficients(species_type='reactants',
                                            reaction=reaction)
    
    def get_product_stoich_coeffs(self, reaction):
        """Helper function that returns product stoichiometric coefficients from input.
        
        Args:
        =====
        reaction : xml.etree.ElementTree.Element, required
            <reaction> XML element containing information about a reaction
        
        Returns:
        ========
        product_stoich_coeffs : dict[str, int]
            Dictionary in the form {product name: stoich coefficient}
        """
        return self.get_stoich_coefficients(species_type='products',
                                            reaction=reaction)

    def populate_reaction_list(self):
        """Populates/updates a list of Reaction or Reaction-inherited objects
        containing information about corresponding reactions."""
        for reactionData in self.rxns_node.findall('reactionData'):
            for reaction in reactionData.findall('reaction'):

                species = self.get_species()
                is_reversible = self.get_is_reversible(reaction) 
                rxn_type = self.get_rxn_type(reaction)
                rxn_equation = self.get_rxn_equation(reaction)
                rate_coeffs_components, rate_coeffs_type = self.get_rate_coeffs_components(reaction)
                reactant_stoich_coeffs = self.get_reactant_stoich_coeffs(reaction)
                product_stoich_coeffs = self.get_product_stoich_coeffs(reaction)

                if is_reversible == False and rxn_type == "Elementary":
                    rxn = IrrevElemReaction(rxn_type=rxn_type, is_reversible=is_reversible,
                                            rxn_equation=rxn_equation,
                                            species_list=species, rate_coeffs_components=rate_coeffs_components,
                                            rate_coeffs_type=rate_coeffs_type,
                                            reactant_stoich_coeffs=reactant_stoich_coeffs,
                                            product_stoich_coeffs=product_stoich_coeffs)

                elif is_reversible == True and rxn_type == "Elementary":
                    rxn = RevElemReaction(rxn_type=rxn_type, is_reversible=is_reversible,
                                            rxn_equation=rxn_equation,
                                            species_list=species, rate_coeffs_components=rate_coeffs_components,
                                            rate_coeffs_type=rate_coeffs_type,
                                            reactant_stoich_coeffs=reactant_stoich_coeffs,
                                            product_stoich_coeffs=product_stoich_coeffs)

                # Unhandled reaction case
                else:
                    raise NotImplementedError("This type of reaction has not been implemented yet!")

                self.reaction_list.append(rxn)
