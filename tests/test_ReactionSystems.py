
"""Test module for """

import numpy
import os
import pytest
import warnings

from pychemkin.config import DB_DIRECTORY
from pychemkin.parsers.XMLParser import XMLParser
from pychemkin.parsers.SQLParser import SQLParser
from pychemkin.reactions.ReactionSystems import ReactionSystem

TEST_DB_PATH = os.path.join(DB_DIRECTORY + "/NASA7_coeffs.sqlite")


# Treat warnings like errors (for testing purposes)
warnings.simplefilter("error")


@pytest.fixture
def test_rxn_sys():
    """Returns a valid reaction system"""
    xml_filename =  "tests/test_xml_files/rxn.xml"
    xml_parser = XMLParser(xml_filename)
    species = xml_parser.get_species()

    sql_parser = SQLParser(TEST_DB_PATH, species)
    thermo_coeffs = sql_parser.get_thermo_coeffs()

    temp = 500 # "low" temperature range in NASA coeffs database
    concentrations = {'H':1, 'O2':1, 'H2O':1}
    rxnsys = ReactionSystem(xml_parser.reaction_list, thermo_coeffs, temp, concentrations)
    return rxnsys

def test_rxn_sys_invalid_temperature():
    """Tests setting up reaction system with invalid temperatures."""
    xml_filename =  "tests/test_xml_files/rxns_mixed.xml"
    xml_parser = XMLParser(xml_filename)
    species = xml_parser.get_species()
    sql_parser = SQLParser(TEST_DB_PATH, species)
    thermo_coeffs = sql_parser.get_thermo_coeffs()
    concentrations = {'H':1, 'O2':2, 'OH':1, 'O':4, 'H2O':0, 'H2':1}
    temp = 0
    with pytest.raises(ValueError):
        rxnsys = ReactionSystem(xml_parser.reaction_list, thermo_coeffs, temp, concentrations)
    temp = -100
    with pytest.raises(ValueError):
        rxnsys = ReactionSystem(xml_parser.reaction_list, thermo_coeffs, temp, concentrations)


def test_rxn_sys_get_reaction_rate_for_1_rxn(test_rxn_sys):
    """Tests function to get reaction rate for a given system of reactions (just 1 reaction)."""
    rates = test_rxn_sys.sort_reaction_rates()
    assert rates['H'] == -30.
    assert rates['O2'] == -15.
    assert rates['H2O'] == 30.


def test_rxn_sys_get_reaction_rate_for_3_rxns():
    """Tests function to get reaction rate for a given system of reactions (more than 1 reaction)."""
    xml_filename =  "tests/test_xml_files/rxnsys.xml"
    xml_parser = XMLParser(xml_filename)
    species = xml_parser.get_species()
    sql_parser = SQLParser(TEST_DB_PATH, species)
    thermo_coeffs = sql_parser.get_thermo_coeffs()
    temp = 10
    concentrations = {'H':1, 'O2':1, 'OH':1, 'O':1, 'H2O':1, 'H2':1}
    rxnsys = ReactionSystem(xml_parser.reaction_list, thermo_coeffs, temp, concentrations)
    rates = rxnsys.sort_reaction_rates()
    assert rates['H'] == -10.
    assert rates['O2'] == -15.
    assert rates['H2O'] == 40.
    assert rates['H2'] == -20.
    assert rates['O'] == -10.
    assert rates['OH'] == 0.


def test_rxn_sys_get_lowT_nasa_matrix(test_rxn_sys):
    """Tests function to fetch NASA coefficients of appropriate T and appropriate species in reaction."""
    # Order of low T range NASA coefficients: H, H2O, O2
    expected_nasa = {'H': numpy.array([2.50000000e+00, 0.00000000e+00,
                                      0.00000000e+00, 0.00000000e+00,
                                      0.00000000e+00, 2.54716270e+04,
                                      -4.60117608e-01]),
                    'H2O': numpy.array([ 3.38684249e+00, 3.47498246e-03,
                                       -6.35469633e-06, 6.96858127e-09,
                                       -2.50658847e-12, -3.02081133e+04,
                                       2.59023285e+00]),
                    'O2': numpy.array([3.21293640e+00, 1.12748635e-03,
                                      -5.75615047e-07, 1.31387723e-09,
                                      -8.76855392e-13, -1.00524902e+03,
                                      6.03473759e+00])}
    assert (numpy.isclose(test_rxn_sys.NASA_matrix['H2O'], expected_nasa['H2O'])).all()
    assert (numpy.isclose(test_rxn_sys.NASA_matrix['O2'], expected_nasa['O2'])).all()
    assert (numpy.isclose(test_rxn_sys.NASA_matrix['H'], expected_nasa['H'])).all()


def test_rxn_sys_get_highT_nasa_matrix():
    """Tests function to fetch NASA coefficients of appropriate T and appropriate species in reaction."""
    xml_filename =  "tests/test_xml_files/rxn.xml"
    xml_parser = XMLParser(xml_filename)
    species = xml_parser.get_species()
    sql_parser = SQLParser(TEST_DB_PATH, species)
    thermo_coeffs = sql_parser.get_thermo_coeffs()
    temp = 5000 # "high" temperature range in NASA coeffs database
    concentrations = {'H':1, 'O2':1, 'H2O':1}
    rxnsys = ReactionSystem(xml_parser.reaction_list, thermo_coeffs, temp, concentrations)
    # Order of high T range NASA coefficients: H, H2O, O2
    expected_nasa = {'O2': numpy.array([3.69757819e+00, 6.13519689e-04,
                                       -1.25884199e-07, 1.77528148e-11,
                                       -1.13643531e-15, -1.23393018e+03,
                                       3.18916559e+00]),
                    'H2O': numpy.array([2.67214561e+00, 3.05629289e-03,
                                       -8.73026011e-07, 1.20099639e-10,
                                       -6.39161787e-15, -2.98992090e+04,
                                       6.86281681e+00]),
                    'H': numpy.array([2.50000000e+00, 0.00000000e+00,
                                     0.00000000e+00, 0.00000000e+00,
                                     0.00000000e+00, 2.54716270e+04,
                                     -4.60117638e-01])}
    assert (numpy.isclose(rxnsys.NASA_matrix['H2O'], expected_nasa['H2O'])).all()
    assert (numpy.isclose(rxnsys.NASA_matrix['O2'], expected_nasa['O2'])).all()
    assert (numpy.isclose(rxnsys.NASA_matrix['H'], expected_nasa['H'])).all()


def test_rxn_sys_rev_reaction():
    """Tests setting up reaction system with reversible reaction."""
    xml_filename =  "tests/test_xml_files/rev_rxn.xml"
    xml_parser = XMLParser(xml_filename)
    species = xml_parser.get_species()
    sql_parser = SQLParser(TEST_DB_PATH, species)
    thermo_coeffs = sql_parser.get_thermo_coeffs()
    temp = 500 # "low" temperature range in NASA coeffs database
    concentrations = {'H':1, 'O2':1, 'H2O':1}
    rxnsys = ReactionSystem(xml_parser.reaction_list, thermo_coeffs, temp, concentrations)
    expected_nasa = {'O2': numpy.array([3.21293640e+00, 1.12748635e-03,
                                       -5.75615047e-07, 1.31387723e-09,
                                       -8.76855392e-13, -1.00524902e+03,
                                       6.03473759e+00]),
                    'H2O': numpy.array([3.38684249e+00, 3.47498246e-03,
                                       -6.35469633e-06, 6.96858127e-09,
                                       -2.50658847e-12, -3.02081133e+04,
                                       2.59023285e+00]),
                    'H': numpy.array([2.50000000e+00, 0.00000000e+00,
                                     0.00000000e+00, 0.00000000e+00,
                                     0.00000000e+00, 2.54716270e+04,
                                     -4.60117608e-01])}
    rev_rxn_obj = xml_parser.reaction_list[0]
    assert (numpy.isclose(rev_rxn_obj.NASA_poly_coefs_dict['H2O'], expected_nasa['H2O'])).all()
    assert (numpy.isclose(rev_rxn_obj.NASA_poly_coefs_dict['O2'], expected_nasa['O2'])).all()
    assert (numpy.isclose(rev_rxn_obj.NASA_poly_coefs_dict['H'], expected_nasa['H'])).all()


def test_rxn_sys_irrev_reaction_antioch1():
    """Test against Antioch irrev rxn results"""
    xml_filename =  "tests/test_xml_files/rxns_irreversible_antioch.xml"
    xml_parser = XMLParser(xml_filename)
    species = xml_parser.get_species()
    sql_parser = SQLParser(TEST_DB_PATH, species)
    thermo_coeffs = sql_parser.get_thermo_coeffs()
    
    # Condition #1
    temp = 2500.0000000000000000
    concentrations = ({'H': 5.0000000000000000e-01,
                      'O': 0.0000000000000000e+00,
                      'OH': 0.0000000000000000e+00,
                      'H2': 2.0000000000000000e+00,
                      'H2O': 0.0000000000000000e+00,
                      'O2': 1.0000000000000000e+00,
                      'HO2': 0.0000000000000000e+00,
                      'H2O2': 0.0000000000000000e+00})
    rxnsys = ReactionSystem(xml_parser.reaction_list, thermo_coeffs, temp, concentrations)
    rates = rxnsys.sort_reaction_rates()

    assert numpy.isclose(rates['H'], -3.3300992971255586e+13)
    assert numpy.isclose(rates['O'], 3.3300992971255586e+13)
    assert numpy.isclose(rates['OH'], 3.3300992971255586e+13)
    assert numpy.isclose(rates['H2'], 0.0000000000000000e+00)
    assert numpy.isclose(rates['H2O'], 0.0000000000000000e+00)
    assert numpy.isclose(rates['O2'], -3.3300992971255586e+13)
    assert numpy.isclose(rates['HO2'], 0.0000000000000000e+00)
    assert numpy.isclose(rates['H2O2'], 0.0000000000000000e+00)

def test_rxn_sys_irrev_reaction_antioch2():
    """Test against Antioch irrev rxn results"""
    xml_filename =  "tests/test_xml_files/rxns_irreversible_antioch.xml"
    xml_parser = XMLParser(xml_filename)
    species = xml_parser.get_species()
    sql_parser = SQLParser(TEST_DB_PATH, species)
    thermo_coeffs = sql_parser.get_thermo_coeffs()

    # Condition #2
    temp = 2500.0000000000000000
    concentrations = ({'H': 5.0000000000000000e-01,
                      'O': 1.0000000000000001e-01,
                      'OH': 1.0000000000000000e-02,
                      'H2': 2.0000000000000000e+00,
                      'H2O': 2.5000000000000000e-01,
                      'O2': 1.0000000000000000e+00,
                      'HO2': 2.9999999999999999e-01,
                      'H2O2': 2.0000000000000000e-02})
    rxnsys = ReactionSystem(xml_parser.reaction_list, thermo_coeffs, temp, concentrations)
    rates = rxnsys.sort_reaction_rates()
    assert numpy.isclose(rates['H'], -3.7324963347340922e+13)
    assert numpy.isclose(rates['O'], 2.3071533925003262e+13)
    assert numpy.isclose(rates['OH'], 6.4368180909475500e+13)
    assert numpy.isclose(rates['H2'], -6.6439941741054521e+12)
    assert numpy.isclose(rates['H2O'], 4.9820020841399396e+11)
    assert numpy.isclose(rates['O2'], -2.9843856969218777e+13)
    assert numpy.isclose(rates['HO2'], -1.3498571473703539e+13)
    assert numpy.isclose(rates['H2O2'], -6.2652907852405969e+11)

def test_rxn_sys_irrev_reaction_antioch3():
    """Test against Antioch irrev rxn results"""
    xml_filename =  "tests/test_xml_files/rxns_irreversible_antioch.xml"
    xml_parser = XMLParser(xml_filename)
    species = xml_parser.get_species()
    sql_parser = SQLParser(TEST_DB_PATH, species)
    thermo_coeffs = sql_parser.get_thermo_coeffs()

    # Condition #3 
    temp = 950.0000000000000000
    concentrations = ({'H': 5.0000000000000000e-01,
                      'O': 0.0000000000000000e+00,
                      'OH':  0.0000000000000000e+00,
                      'H2': 2.0000000000000000e+00,
                      'H2O': 0.0000000000000000e+00,
                      'O2': 1.0000000000000000e+00,
                      'HO2': 0.0000000000000000e+00,
                      'H2O2': 0.0000000000000000e+00})
    rxnsys = ReactionSystem(xml_parser.reaction_list, thermo_coeffs, temp, concentrations)
    rates = rxnsys.sort_reaction_rates()
    assert numpy.isclose(rates['H'], -1.3403448555187156e+13)
    assert numpy.isclose(rates['O'], 1.3403448555187156e+13)
    assert numpy.isclose(rates['OH'], 1.3403448555187156e+13)
    assert numpy.isclose(rates['H2'], 0.0000000000000000e+00)
    assert numpy.isclose(rates['H2O'], 0.0000000000000000e+00)
    assert numpy.isclose(rates['O2'], -1.3403448555187156e+13)
    assert numpy.isclose(rates['HO2'], 0.0000000000000000e+00)
    assert numpy.isclose(rates['H2O2'], 0.0000000000000000e+00)

def test_rxn_sys_irrev_reaction_antioch4():
    """Test against Antioch irrev rxn results"""
    xml_filename =  "tests/test_xml_files/rxns_irreversible_antioch.xml"
    xml_parser = XMLParser(xml_filename)
    species = xml_parser.get_species()
    sql_parser = SQLParser(TEST_DB_PATH, species)
    thermo_coeffs = sql_parser.get_thermo_coeffs()

    # Condition #4
    temp = 950.0000000000000000
    concentrations = ({'H': 5.0000000000000000e-01,
                      'O': 1.0000000000000001e-01,
                      'OH': 1.0000000000000000e-02,
                      'H2': 2.0000000000000000e+00,
                      'H2O': 2.5000000000000000e-01,
                      'O2': 1.0000000000000000e+00,
                      'HO2': 2.9999999999999999e-01,
                      'H2O2': 2.0000000000000000e-02})
    rxnsys = ReactionSystem(xml_parser.reaction_list, thermo_coeffs, temp, concentrations)
    rates = rxnsys.sort_reaction_rates()
    assert numpy.isclose(rates['H'], -2.5701654161377395e+13)
    assert numpy.isclose(rates['O'], 1.1995070447537494e+13)
    assert numpy.isclose(rates['OH'], 3.5250102331863312e+13)
    assert numpy.isclose(rates['H2'], 1.9231756318330654e+12)
    assert numpy.isclose(rates['H2O'], 3.1178427968785107e+11)
    assert numpy.isclose(rates['O2'], -1.0092501951521918e+13)
    assert numpy.isclose(rates['HO2'], -1.3353585162517070e+13)
    assert numpy.isclose(rates['H2O2'], -3.3239141550534125e+11)


def test_rxn_sys_rev_reaction_antioch1():
    """Test against Antioch rev rxn results"""
    xml_filename =  "tests/test_xml_files/rxns_reversible_antioch.xml"
    xml_parser = XMLParser(xml_filename)
    species = xml_parser.get_species()
    sql_parser = SQLParser(TEST_DB_PATH, species)
    thermo_coeffs = sql_parser.get_thermo_coeffs()

    # Condition #1
    temp = 2500.0000000000000000
    concentrations = ({'H': 5.0000000000000000e-01,
                      'O': 0.0000000000000000e+00,
                      'OH': 0.0000000000000000e+00,
                      'H2': 2.0000000000000000e+00,
                      'H2O': 0.0000000000000000e+00,
                      'O2': 1.0000000000000000e+00,
                      'HO2': 0.0000000000000000e+00,
                      'H2O2': 0.0000000000000000e+00})
    rxnsys = ReactionSystem(xml_parser.reaction_list, thermo_coeffs, temp, concentrations)
    rates = rxnsys.sort_reaction_rates()
    assert numpy.isclose(rates['H'], -3.3299453222970820e+13)
    assert numpy.isclose(rates['O'], 3.3300992971255586e+13)
    assert numpy.isclose(rates['OH'], 3.3300992971255586e+13)
    assert numpy.isclose(rates['H2'], -1.5397482847664728e+09)
    assert numpy.isclose(rates['H2O'], 0.0000000000000000e+00)
    assert numpy.isclose(rates['O2'], -3.3302532719540352e+13)
    assert numpy.isclose(rates['HO2'], 1.5397482847664728e+09)
    assert numpy.isclose(rates['H2O2'], 0.0000000000000000e+00)

def test_rxn_sys_rev_reaction_antioch2():
    """Test against Antioch rev rxn results"""
    xml_filename =  "tests/test_xml_files/rxns_reversible_antioch.xml"
    xml_parser = XMLParser(xml_filename)
    species = xml_parser.get_species()
    sql_parser = SQLParser(TEST_DB_PATH, species)
    thermo_coeffs = sql_parser.get_thermo_coeffs()

    # Condition #2
    temp = 2500.0000000000000000
    concentrations = ({'H': 5.0000000000000000e-01,
                      'O': 1.0000000000000001e-01,
                      'OH': 1.0000000000000000e-02,
                      'H2': 2.0000000000000000e+00,
                      'H2O': 2.5000000000000000e-01,
                      'O2': 1.0000000000000000e+00,
                      'HO2': 2.9999999999999999e-01,
                      'H2O2': 2.0000000000000000e-02})
    rxnsys = ReactionSystem(xml_parser.reaction_list, thermo_coeffs, temp, concentrations)
    rates = rxnsys.sort_reaction_rates()
    assert numpy.isclose(rates['H'], -3.7739607416516375e+13)
    assert numpy.isclose(rates['O'], 2.3087395959996922e+13)
    assert numpy.isclose(rates['OH'], 6.4832447404435727e+13)
    assert numpy.isclose(rates['H2'], -6.1109243594152285e+12)
    assert numpy.isclose(rates['H2O'], -2.1877981256944708e+11)
    assert numpy.isclose(rates['O2'], -2.9727063574790047e+13)
    assert numpy.isclose(rates['HO2'], -1.3813504758333098e+13)
    assert numpy.isclose(rates['H2O2'], -3.0996344280845251e+11)

def test_rxn_sys_rev_reaction_antioch3():
    """Test against Antioch rev rxn results"""
    xml_filename =  "tests/test_xml_files/rxns_reversible_antioch.xml"
    xml_parser = XMLParser(xml_filename)
    species = xml_parser.get_species()
    sql_parser = SQLParser(TEST_DB_PATH, species)
    thermo_coeffs = sql_parser.get_thermo_coeffs()

    # Condition #3
    temp = 950.0000000000000000
    concentrations = ({'H': 5.0000000000000000e-01,
                      'O': 0.0000000000000000e+00,
                      'OH':  0.0000000000000000e+00,
                      'H2': 2.0000000000000000e+00,
                      'H2O': 0.0000000000000000e+00,
                      'O2': 1.0000000000000000e+00,
                      'HO2': 0.0000000000000000e+00,
                      'H2O2': 0.0000000000000000e+00})
    rxnsys = ReactionSystem(xml_parser.reaction_list, thermo_coeffs, temp, concentrations)
    rates = rxnsys.sort_reaction_rates()
    assert numpy.isclose(rates['H'], -1.3403448555170947e+13)
    assert numpy.isclose(rates['O'], 1.3403448555187156e+13)
    assert numpy.isclose(rates['OH'], 1.3403448555187156e+13)
    assert numpy.isclose(rates['H2'], -1.6208366095320475e+01)
    assert numpy.isclose(rates['H2O'], 0.0000000000000000e+00)
    assert numpy.isclose(rates['O2'], -1.3403448555203365e+13)
    assert numpy.isclose(rates['HO2'], 1.6208366095320475e+01)
    assert numpy.isclose(rates['H2O2'], 0.0000000000000000e+00)

def test_rxn_sys_rev_reaction_antioch4():
    """Test against Antioch rev rxn results"""
    xml_filename =  "tests/test_xml_files/rxns_reversible_antioch.xml"
    xml_parser = XMLParser(xml_filename)
    species = xml_parser.get_species()
    sql_parser = SQLParser(TEST_DB_PATH, species)
    thermo_coeffs = sql_parser.get_thermo_coeffs()

    # Condition #4
    temp = 950.0000000000000000
    concentrations = ({'H': 5.0000000000000000e-01,
                      'O': 1.0000000000000001e-01,
                      'OH': 1.0000000000000000e-02,
                      'H2': 2.0000000000000000e+00,
                      'H2O': 2.5000000000000000e-01,
                      'O2': 1.0000000000000000e+00,
                      'HO2': 2.9999999999999999e-01,
                      'H2O2': 2.0000000000000000e-02})
    rxnsys = ReactionSystem(xml_parser.reaction_list, thermo_coeffs, temp, concentrations)
    rates = rxnsys.sort_reaction_rates()
    assert numpy.isclose(rates['H'], -1.7750736729894043e+13)
    assert numpy.isclose(rates['O'], 4.0714901298639282e+12)
    assert numpy.isclose(rates['OH'], 2.7224152308533266e+13)
    assert numpy.isclose(rates['H2'], 1.9335776122987725e+12)
    assert numpy.isclose(rates['H2O'], 3.3867579679335474e+11)
    assert numpy.isclose(rates['O2'], -2.1311825395902986e+12)
    assert numpy.isclose(rates['HO2'], -1.3354030759186477e+13)
    assert numpy.isclose(rates['H2O2'], -3.3194581881849854e+11)
