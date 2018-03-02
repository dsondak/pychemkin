
"""Classes for preprocessing: parsing xml files, identifying reaction type."""

import csv
import numpy
import os
import xml.etree.ElementTree as ET
from pychemkin.config import THIS_DIRECTORY
from pychemkin.reactions.Reactions import *


class XMLParser:
    """Parser for input XML files to retrieve and
    preprocess reaction data."""
    def __init__(self, xml_filename, convert_units=False):
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
        self.convert_units = convert_units
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
        convert_units : boolean, optional (default False)
            converts units if True


        RETURNS:
        --------
        rate_coeffs_components : dict
            dictionary of the form {coefficient component name: coefficient component value}. 
        """
        if self.convert_units:
            # Connect to csv file containing units
            units_file = os.path.join(THIS_DIRECTORY, 'units.csv')
            with open(units_file, 'r') as unit:
                next(unit)
                unit_dict = dict(csv.reader(unit))
            # Create dictionary of units
            dict_ = dict((k, float(v)) for k, v in unit_dict.items())

        # Loop over rateCoeff's and convert units where desired
        rateCoeffs = reaction.find('rateCoeff')

        #output = []
        for rateCoeff in rateCoeffs:
            if rateCoeff.tag == 'Arrhenius':
                # If 'Arrhenius' units are to be converted
                if self.convert_units:
                    try:
                        A_unit = rateCoeff.find('A').attrib['units'].split('/')
                        E_unit = rateCoeff.find('E').attrib['units'].split('/')
                    except:
                        raise ValueError("Input file contains no units. " + 
                                         "Set convert_units to False to continue")
                    A_conv_lis = []
                    for unit in A_unit:
                        try:
                            A_conv_lis.append(dict_[unit])
                        except:
                            raise NotImplementedError(unit +
                                                      " not implemented.")
                    A_conversion = numpy.prod(numpy.array(A_conv_lis))
                    E_conv_lis = []
                    for unit in E_unit:
                        try:
                            E_conv_lis.append(dict_[unit])
                        except:
                            raise NotImplementedError(unit +
                                                      " not implemented.")
                    E_conversion = numpy.prod(numpy.array(E_conv_lis))
                    try:
                        A = float(rateCoeff.find('A').text)*A_conversion
                        E = float(rateCoeff.find('E').text)*E_conversion
                        d = {'A': A, 'E': E}
                    except:
                        print("Conversion failed. " +
                              "Resulting units have not been converted.")
                        A = float(rateCoeff.find('A').text)
                        E = float(rateCoeff.find('E').text)
                        d = {'A': A, 'E': E}
                    if rateCoeff.find('b') is not None:
                        raise ValueError("Cannot use 'b' in Arrhenius type.")
                # If 'Arrhenius' units are not to be converted
                if not self.convert_units:
                    try:
                        A = float(rateCoeff.find('A').text)
                        E = float(rateCoeff.find('E').text)
                        d = {'A': A, 'E': E}
                    except:
                        raise ValueError("Reaction coefficient parameters " +
                                         "not as expected.")
                    if rateCoeff.find('b') is not None:
                        raise ValueError("Cannot use 'b' in Arrhenius type.")

            elif rateCoeff.tag in ('modifiedArrhenius',
                                   'Kooij'):
                kooij_name = ''
                try:
                    kooij_name = rateCoeff.attrib['name']
                except:
                    kooij_name = None
                # If 'modifiedArrhenius' units are to be converted
                if self.convert_units:
                    try:
                        A_unit = rateCoeff.find('A').attrib['units'].split('/')
                        E_unit = rateCoeff.find('E').attrib['units'].split('/')
                    except:
                        raise ValueError("Input file contains no units. " + 
                                         "Set convert_units to False to continue")
                    A_conv_lis = []
                    for unit in A_unit:
                        try:
                            A_conv_lis.append(dict_[unit])
                        except:
                            raise NotImplementedError(unit +
                                                      " not implemented.")
                    A_conversion = numpy.prod(numpy.array(A_conv_lis))
                    E_conv_lis = []
                    for unit in E_unit:
                        try:
                            E_conv_lis.append(dict_[unit])
                        except:
                            raise NotImplementedError(unit +
                                                      " not implemented.")
                    E_conversion = numpy.prod(numpy.array(E_conv_lis))
                    try:
                        A = float(rateCoeff.find('A').text)*A_conversion
                        b = float(rateCoeff.find('b').text)
                        E = float(rateCoeff.find('E').text)*E_conversion
                        d = {'A': A, 'b': b, 'E': E}
                        if rateCoeff.tag == 'Kooij':
                            d['name'] = kooij_name
                    except:
                        print("Conversion failed. " +
                              "Resulting units have not been converted.")
                        A = float(rateCoeff.find('A').text)
                        b = float(rateCoeff.find('b').text)
                        E = float(rateCoeff.find('E').text)
                        d = {'A': A, 'b': b, 'E': E}
                        if rateCoeff.tag == 'Kooij':
                            d['name'] = kooij_name
                # If 'modifiedArrhenius' units are not to be converted
                if not self.convert_units:
                    try:
                        A = float(rateCoeff.find('A').text)
                        b = float(rateCoeff.find('b').text)
                        E = float(rateCoeff.find('E').text)
                        d = {'A': A, 'b': b, 'E': E}
                        if rateCoeff.tag == 'Kooij':
                            d['name'] = kooij_name
                    except:
                        raise ValueError("Reaction coefficient parameters " +
                                         "not as expected.")

            elif rateCoeff.tag == 'Constant':
                try:
                    k = float(rateCoeff.find('k').text)
                    d = {'k': k}
                except:
                    raise ValueError("Non-numeric coefficient parameters.")

            elif rateCoeff.tag == 'efficiencies':
                temp_list = [item.split(':') for item
                             in rateCoeff.text.split(' ')]
                efficiencies = dict()
                efficiencies['Type'] = 'efficiencies'
                for item in temp_list:
                    efficiencies[item[0]] = item[1]
                efficiencies['default'] = rateCoeff.attrib['default']
                d = efficiencies

            elif rateCoeff.tag == 'Troe':
                alpha = float(rateCoeff.find('alpha').text)
                t1 = float(rateCoeff.find('T1').text)
                t2 = float(rateCoeff.find('T2').text)
                t3 = float(rateCoeff.find('T2').text)

                d = {'alpha': alpha,
                     't1': t1,
                     't2': t2,
                     't3': t3}

            else:
                raise NotImplementedError(rateCoeff.tag + " not implemented.")

            #output.append(d)
        #return output
        return d

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
