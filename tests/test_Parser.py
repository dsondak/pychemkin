
"""Tests for XmlParser and RxnData classes."""

import pytest

from pychemkin.config import pkg_xml_path
from pychemkin.preprocessing.Parser import XmlParser, RxnType
from pychemkin.pychemkin_errors import PyChemKinError


def test_parse_basic_functionality():
    """Tests attributes of XML file are parsed correctly."""
    xml = XmlParser(pkg_xml_path('rxns_ideal.xml'))
    species, rxns = xml.load()

    # Correct number of reactions returned.
    err_msg = 'Expected 2 reactions but received {}.'.format(len(rxns))
    assert len(rxns) == 2, err_msg

    # Attributes of <reaction> element parsed properly.
    err_msg = 'reversible attribute not parsed properly.'
    assert rxns[0].is_reversible is False, err_msg
    assert rxns[1].is_reversible is False, err_msg

    err_msg = 'rxn_id attribute not parsed properly.'
    assert rxns[0].rxn_id == 'reaction01', err_msg
    assert rxns[1].rxn_id == 'reaction02', err_msg

    err_msg = 'type attribute not parsed properly.'
    assert rxns[0].type == RxnType.Elementary, err_msg
    assert rxns[1].type == RxnType.Elementary, err_msg


def test_parse_reactants_products():
    """Tests parsing of reactants and products."""
    xml = XmlParser(pkg_xml_path('rxns_ideal.xml'))
    species, rxns = xml.load()

    err_msg = 'reactants not parsed correctly.'
    assert rxns[0].reactants == {'H': 1, 'O2': 1}, err_msg
    assert rxns[1].reactants == {'H2': 1, 'O': 1}, err_msg

    err_msg = 'products not parsed correctly.'
    assert rxns[0].products == {'OH': 1, 'O': 1}, err_msg
    assert rxns[1].products == {'OH': 1, 'H': 1}, err_msg


def test_parse_rxn_coeff():
    """Tests parsing of reaction rate coefficients."""
    xml = XmlParser(pkg_xml_path('rxns.xml'))
    species, rxns = xml.load()

    rxn1_coeff = rxns[0].rate_coeff
    rxn2_coeff = rxns[1].rate_coeff
    rxn3_coeff = rxns[2].rate_coeff

    # reaction01
    assert rxn1_coeff[0] == pytest.approx(3.52e+10)
    assert rxn1_coeff[1] == pytest.approx(7.14e+04)

    # reaction02
    assert rxn2_coeff[0] == pytest.approx(5.06e-2)
    assert rxn2_coeff[1] == pytest.approx(2.7)
    assert rxn2_coeff[2] == pytest.approx(2.63e+04)

    # reaction03
    assert rxn3_coeff == pytest.approx(1.0e+03)


def test_rxndata_equation():
    """Tests equation representation of RxnData."""
    xml = XmlParser(pkg_xml_path('rxns_ideal.xml'))
    species, rxns = xml.load()

    err_msg = 'get_equation() method result different than expected.'
    expected_0 = 'H + O2 =] O + OH'
    assert rxns[0].get_equation() == expected_0, err_msg

    err_msg = 'get_equation() method result different than expected.'
    expected_1 = 'H2 + O =] H + OH'
    assert rxns[1].get_equation() == expected_1, err_msg


def test_badparse_negative_A_arr():
    """Tests for case when Arrhenius coefficient component A for a
    reaction is negative.
    """
    xml = XmlParser(pkg_xml_path('rxns_neg_A_2'))
    try:
        rxns = xml.load()
    except PyChemKinError as err:
        assert type(err) == PyChemKinError
        assert str(err).find(
              'A coeff < 0 in reaction with id = reaction02') != -1


def test_badparse_negative_A_modarr():
    """Tests for case when modified Arrhenius coefficient component
    A for a reaction is negative.
    """
    xml = XmlParser(pkg_xml_path('rxns_neg_A'))
    try:
        rxns = xml.load()
    except PyChemKinError as err:
        assert type(err) == PyChemKinError
        assert str(err).find(
              'A coeff < 0 in reaction with id = reaction01') != -1
