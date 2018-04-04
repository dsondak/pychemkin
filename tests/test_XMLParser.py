
"""Tests for XMLParser"""

import numpy
import pytest
import os

from pychemkin.parsers.XMLParser import XMLParser


def test_XMLParser_file_not_found():
    """Test when XML file is nonexistent."""
    with pytest.raises(OSError):
        parser = XMLParser("no_such_file")

def test_unhandled_xml_components():
    """Test when parts of reaction in the XML file
    are unhandled/not implemented."""

    # Test when type of reaction is unhandled
    xml_filename = "tests/test_xml_files/unhandled_rxn.xml"
    with pytest.raises(NotImplementedError):
        parser = XMLParser(xml_filename)

    # Test when reaction rate coefficient is unhandled
    xml_filename = "tests/test_xml_files/unhandled_k.xml"
    with pytest.raises(NotImplementedError):
        parser = XMLParser(xml_filename)

    # Test when units are undhandled
    with pytest.raises(NotImplementedError):
        xml_filename = "tests/test_xml_files/madeup_units_4_A_arr.xml"
        parser = XMLParser(xml_filename, convert_to_SI_units=True)

    with pytest.raises(NotImplementedError):
        xml_filename = "tests/test_xml_files/madeup_units_4_A_mod_arr.xml"
        parser = XMLParser(xml_filename, convert_to_SI_units=True)

    with pytest.raises(NotImplementedError):
        xml_filename = "tests/test_xml_files/madeup_units_4_E_arr.xml"
        parser = XMLParser(xml_filename, convert_to_SI_units=True)

    with pytest.raises(NotImplementedError):
        xml_filename = "tests/test_xml_files/madeup_units_4_E_mod_arr.xml"
        parser = XMLParser(xml_filename, convert_to_SI_units=True)

def test_XMLParser_functionality():
    """Test various parts of XMLParser."""

    xml_filename = "tests/test_xml_files/rxns.xml"
    parser = XMLParser(xml_filename)

    # get_species() routine
    assert parser.get_species() == ({'H': None,'O': None, 'OH': None,
                                    'H2': None, 'H2O': None, 'O2': None})

    # get_rxn_type() routine
    assert parser.reaction_list[0].rxn_type == 'Elementary'

    # get_rate_coeffs_components() routine
    assert (parser.reaction_list[0].rate_coeffs_components ==
            {'A': 35200000000.0, 'E': 71400.0})

    # get_is_reversible() routine
    assert parser.reaction_list[0].is_reversible == False

    # get_rxn_equation() routine
    assert parser.reaction_list[0].rxn_equation == 'H + O2 =] OH + O'

    # get_reactant_stoich_coeffs() routine
    assert (parser.reaction_list[0].reactant_stoich_coeffs ==
            {'H': 1, 'H2': 0, 'H2O': 0, 'O': 0, 'O2': 1, 'OH': 0})

    # get_product_stoich_coeffs() routine
    assert (parser.reaction_list[0].product_stoich_coeffs ==
            {'H': 0, 'H2': 0, 'H2O': 0, 'O': 1, 'O2': 0, 'OH': 1})

def test_xml_files_with_missing_info():
    """Test various cases in which input XML files
    are missing information."""

    # Test when k is missing from constant type reaction
    with pytest.raises(ValueError):
        xml_filename = "tests/test_xml_files/k_const.xml"
        parser = XMLParser(xml_filename)

    # Test when A is missing from Arrhenius type reaction
    with pytest.raises(ValueError):
        xml_filename = "tests/test_xml_files/A_arr.xml"
        parser = XMLParser(xml_filename)

    # Test when E is missing from Arrhenius type reaction
    with pytest.raises(ValueError):
        xml_filename = "tests/test_xml_files/E_arr.xml"
        parser = XMLParser(xml_filename)

    # Test when A is missing from modified Arrhenius type reaction
    with pytest.raises(ValueError):
        xml_filename = "tests/test_xml_files/A_mod_arr.xml"
        parser = XMLParser(xml_filename)

    # Test when b is missing from modified Arrhenius type reaction
    with pytest.raises(ValueError):
        xml_filename = "tests/test_xml_files/b_mod_arr.xml"
        parser = XMLParser(xml_filename)

    # Test when E is missing from modified Arrhenius type reaction
    with pytest.raises(ValueError):
        xml_filename = "tests/test_xml_files/E_mod_arr.xml"
        parser = XMLParser(xml_filename)

def test_convert_to_SI_units_when_no_units():
    """Test when set convert_to_SI_units to True but no units found in XML."""
    with pytest.raises(ValueError):
        xml_filename = "tests/test_xml_files/A_arr.xml"
        parser = XMLParser(xml_filename, convert_to_SI_units=True)

    with pytest.raises(ValueError):
        xml_filename = "tests/test_xml_files/A_mod_arr.xml"
        parser = XMLParser(xml_filename, convert_to_SI_units=True)

def test_faulty_b_in_arr():
    """Test when b found in Arrhenius (vs. modified Arrhenius)."""
    with pytest.raises(ValueError):
        xml_filename = "tests/test_xml_files/faulty_A_arr.xml"
        parser = XMLParser(xml_filename)

    with pytest.raises(ValueError):
        xml_filename = "tests/test_xml_files/faulty_A_arr.xml"
        parser = XMLParser(xml_filename, convert_to_SI_units=True)

def test_unit_checks():
    """Tests for unit checks."""

    # Arrhenius-type reaction without R
    xml_filename = "tests/test_xml_files/unit_check_arr.xml"
    parser = XMLParser(xml_filename, convert_to_SI_units=True)
    A = parser.reaction_list[0].rate_coeffs_components['A']
    E = parser.reaction_list[0].rate_coeffs_components['E']
    assert numpy.isclose(A, 35200, atol=1e-16)
    assert numpy.isclose(E, 298737.6, atol=1e-16)

    # Arrhenius-type reaction with R
    xml_filename = "tests/test_xml_files/unit_check_arr_with_R.xml"
    parser = XMLParser(xml_filename, convert_to_SI_units=False)
    A = parser.reaction_list[0].rate_coeffs_components['A']
    E = parser.reaction_list[0].rate_coeffs_components['E']
    R = parser.reaction_list[0].rate_coeffs_components['R']
    assert numpy.isclose(A, 3.52e+10, atol=1e-16)
    assert numpy.isclose(E, 7.14e+04, atol=1e-16)
    assert numpy.isclose(R, 8.3144598, atol=1e-16)

    # Modified Arrhenius-type reaction without R
    xml_filename = "tests/test_xml_files/unit_check_modarr.xml"
    parser = XMLParser(xml_filename, convert_to_SI_units=True)
    A = parser.reaction_list[0].rate_coeffs_components['A']
    E = parser.reaction_list[0].rate_coeffs_components['E']
    b = parser.reaction_list[0].rate_coeffs_components['b']
    assert numpy.isclose(A, 35200, atol=1e-16)
    assert numpy.isclose(E, 298737.6, atol=1e-16)
    assert numpy.isclose(b, 2.7, atol=1e-16)

    # Modified Arrhenius-type reaction with R
    xml_filename = "tests/test_xml_files/unit_check_modarr_with_R.xml"
    parser = XMLParser(xml_filename, convert_to_SI_units=False)
    A = parser.reaction_list[0].rate_coeffs_components['A']
    E = parser.reaction_list[0].rate_coeffs_components['E']
    b = parser.reaction_list[0].rate_coeffs_components['b']
    R = parser.reaction_list[0].rate_coeffs_components['R']
    assert numpy.isclose(A, 3.52e+10, atol=1e-16)
    assert numpy.isclose(E, 7.14e+04, atol=1e-16)
    assert numpy.isclose(b, 2.7, atol=1e-16)
    assert numpy.isclose(R, 8.3144598, atol=1e-16)

def test_unit_conversion_fail_arr_onlyE_units():
    """Test when only a few (i.e. not all) components have units."""
    xml_filename = "tests/test_xml_files/unit_conversion_fail_arr_onlyEunits.xml"
    with pytest.raises(ValueError):
        parser = XMLParser(xml_filename, convert_to_SI_units=True)

def test_invalid_component_values():
    """Tests when fields of XML files contain invalid values (non-float)."""

    # Test Arrhenius-type reaction when nonsensical R value (e.g. string)
    xml_filename = "tests/test_xml_files/unit_conversion_fail_arr_invalidR.xml"
    with pytest.raises(ValueError):
        parser = XMLParser(xml_filename, convert_to_SI_units=False)

    # Test Arrhenius-type reaction when nonsensical A value (e.g. string)
    xml_filename = "tests/test_xml_files/unit_conversion_fail_arr_A.xml"
    with pytest.raises(ValueError):
        parser = XMLParser(xml_filename, convert_to_SI_units=True)

    # Test Arrhenius-type reaction when nonsensical E value (e.g. string)
    xml_filename = "tests/test_xml_files/unit_conversion_fail_arr_E.xml"
    with pytest.raises(ValueError):
        parser = XMLParser(xml_filename, convert_to_SI_units=True)

    # Test modified Arrhenius-type reaction when nonsensical R value (e.g. string)
    xml_filename = "tests/test_xml_files/unit_conversion_fail_modarr_invalidR.xml"
    with pytest.raises(ValueError):
        parser = XMLParser(xml_filename, convert_to_SI_units=False)

    # Test modified Arrhenius-type reaction when nonsensical A value (e.g. string)
    xml_filename = "tests/test_xml_files/unit_conversion_fail_modarr_A.xml"
    with pytest.raises(ValueError):
        parser = XMLParser(xml_filename, convert_to_SI_units=True)

    # Test modified Arrhenius-type reaction when nonsensical b value (e.g. string)        
    xml_filename = "tests/test_xml_files/unit_conversion_fail_modarr_b.xml"
    with pytest.raises(ValueError):
        parser = XMLParser(xml_filename, convert_to_SI_units=True)

    # Test modified Arrhenius-type reaction when nonsensical E value (e.g. string)
    xml_filename = "tests/test_xml_files/unit_conversion_fail_modarr_E.xml"
    with pytest.raises(ValueError):
        parser = XMLParser(xml_filename, convert_to_SI_units=True)
