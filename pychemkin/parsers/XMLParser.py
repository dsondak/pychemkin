
"""Classes for preprocessing: parsing xml files, identifying reaction type."""

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
            raise OSError("Reaction (xml) file not found!")

        self.reaction_list = []
        self.species = {}
        self.get_species()
        self.populate_reaction_list()

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

            else:
                raise NotImplementedError("This reaction rate coefficient type has not been handled!")

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

    def populate_reaction_list(self):
        """Populates a list of Reaction or Reaction-inherited objects
        containing information about corresponding reactions."""
        for reactionData in self.rxns.findall('reactionData'):
            for reaction in reactionData.findall('reaction'):
                
                species = self.get_species()
                is_reversible = self.get_is_reversible(reaction) 
                rxn_type = self.get_rxn_type(reaction)
                rxn_equation = self.get_rxn_equation(reaction)
                rate_coeffs_components = self.get_rate_coeffs_components(reaction)
                reactant_stoich_coeffs = self.get_reactant_stoich_coeffs(reaction)
                product_stoich_coeffs = self.get_product_stoich_coeffs(reaction)

                if is_reversible == False and rxn_type == "Elementary":
                    rxn = IrrevElemReaction(rxn_type, is_reversible, rxn_equation,
                                            species, rate_coeffs_components,
                                            reactant_stoich_coeffs, product_stoich_coeffs)
                    self.reaction_list.append(rxn)

                elif is_reversible == True and rxn_type == "Elementary":
                    rxn = RevElemReaction(rxn_type, is_reversible, rxn_equation,
                                          species, rate_coeffs_components,
                                          reactant_stoich_coeffs, product_stoich_coeffs)
                    self.reaction_list.append(rxn)

                # Unhandled reaction case
                else:
                    raise NotImplementedError("This type of reaction has not been implemented yet!")
