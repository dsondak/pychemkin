<<<<<<< HEAD
=======

"""Test module for reactions."""

>>>>>>> origin/develop
import numpy

import warnings
warnings.simplefilter("error")

import pytest
from pychemkin.reactions.Reactions import *
<<<<<<< HEAD


=======


>>>>>>> origin/develop
# ===== Tests for base elementary reactions ===== #

@pytest.fixture
def test_base_reaction():
    """Returns a valid reaction setup (from rxns.xml)"""
    return ElementaryReaction(rxn_type="Elementary",
                              is_reversible=False,
                              rxn_equation="H2 + OH =] H2O + H",
                              species_list=['H', 'O', 'OH', 'H2', 'H2O', 'O2'],
                              rate_coeffs_components={'k': 10},
                              rate_coeffs_type='constant',
                              reactant_stoich_coeffs={'H2' :1, 'OH':1},
                              product_stoich_coeffs={'H2O' :1, 'H':1})

def test_Reaction_functionalities(test_base_reaction):
    """Test functions for Reaction objects."""

    # Special functions
    expected_reaction = "Reaction : H2 + OH =] H2O + H"
    assert str(test_base_reaction) == expected_reaction
    assert len(test_base_reaction) == 4

    # get_unique_species() routine
    assert 'H2O' in test_base_reaction.unique_species
    assert 'H2' in test_base_reaction.unique_species
    assert 'H' in test_base_reaction.unique_species
    assert 'OH' in test_base_reaction.unique_species

    # set_temperature() routine
    test_base_reaction.set_temperature(10)
    assert test_base_reaction.temperature == 10

    # get_ordered_list() routine
    ordering_dict = {'H2': 1, 'OH': 2, 'H2O': 3, 'H': 4}
    expected_ordered_list = [4, 2, 1, 3]
    test_list = test_base_reaction.get_ordered_list(ordering_dict)
    assert test_list == expected_ordered_list

    # set_concentrations() routine
    expected_concentrations = [4, 2, 1, 3]
    test_base_reaction.set_concentrations({'H2':1, 'OH':2, 'H2O':3, 'H':4})
    assert (test_base_reaction.concentrations == expected_concentrations).all()

    # compute_reaction_rate_coeff() routine
    with pytest.raises(NotImplementedError):
        test_base_reaction.compute_reaction_rate_coeff()

    # compute_progress_rate() routine
    with pytest.raises(NotImplementedError):
        test_base_reaction.compute_progress_rate()

    # compute_reaction_rate() routine
    with pytest.raises(ElementaryReactionError):
        test_base_reaction.compute_reaction_rate()


def test_Reaction_set_invalid_temperatures(test_base_reaction):
    """Test when trying to set invalid temperatures (in Kelvin)."""

    # Absolute zero
    with pytest.raises(ValueError):
        test_base_reaction.set_temperature(0)
<<<<<<< HEAD

    # Negative temperature
    with pytest.raises(ValueError):
        test_base_reaction.set_temperature(-100)

def test_Reaction_set_invalid_concentrations(test_base_reaction):
    """Test when trying to set invalid concentrations."""

    # Unrecognized species (CH4)
    with pytest.raises(KeyError):
        test_base_reaction.set_concentrations({'H2':1, 'CH4':2, 'H2O':3, 'H':4})
=======

    # Negative temperature
    with pytest.raises(ValueError):
        test_base_reaction.set_temperature(-100)

def test_Reaction_set_invalid_concentrations(test_base_reaction):
    """Test when trying to set invalid concentrations."""

    # Unrecognized species (CH4)
    with pytest.raises(KeyError):
        test_base_reaction.set_concentrations({'H2':1, 'CH4':2, 'H2O':3, 'H':4})

    # Negative concentration values
    with pytest.raises(ValueError):
        test_base_reaction.set_concentrations({'H2':1, 'OH':2, 'H2O':3, 'H':-4})
>>>>>>> origin/develop

    # Negative concentration values
    with pytest.raises(ValueError):
        test_base_reaction.set_concentrations({'H2':1, 'OH':2, 'H2O':3, 'H':-4})

<<<<<<< HEAD

=======
>>>>>>> origin/develop
# ===== Tests for irreversible elementary reactions ===== #

@pytest.fixture
def test_irrev_reaction():
    """Returns a valid reaction (from rxns.xml)"""
    return IrrevElemReaction(rxn_type="Elementary",
                                is_reversible=False,
                                rxn_equation="A + B =] C",
                                species_list=['A', 'B', 'C'],
                                rate_coeffs_components={'k': 10},
                                rate_coeffs_type='constant',
                                reactant_stoich_coeffs={'A':2, 'B':1},
                                product_stoich_coeffs={'C': 1})

@pytest.fixture
def test_irrev_reaction_arrhenius():
    """Returns a reaction with Arrhenius reaction rate coefficient."""
    return IrrevElemReaction(rxn_type="Elementary",
                             is_reversible=False,
                             rxn_equation="H2 + OH =] H2O + H",
                             species_list=['H', 'O', 'OH', 'H2', 'H2O', 'O2'],
                             rate_coeffs_components={'A': 10, 'E': 100, 'R':8.3144598},
                             rate_coeffs_type='arrhenius',
                             reactant_stoich_coeffs={'H2' :1, 'OH':1},
                             product_stoich_coeffs={'H2O' :1, 'H':1})

@pytest.fixture
def test_irrev_reaction_modified_arr():
    """Returns a reaction with modified Arrhenius reaction rate coefficient."""
    return IrrevElemReaction(rxn_type="Elementary",
                             is_reversible=False,
                             rxn_equation="H2 + OH =] H2O + H",
                             species_list=['H', 'O', 'OH', 'H2', 'H2O', 'O2'],
                             rate_coeffs_components={'A': 10, 'E': 100, 'b':0.5, 'R':8.3144598},
                             rate_coeffs_type='modified arrhenius',
                             reactant_stoich_coeffs={'H2' :1, 'OH':1},
                             product_stoich_coeffs={'H2O' :1, 'H':1})

def test_wrongly_classified_irrev_elem_rxns():
    """Tests for wrongly classified irreversible elementary reactions."""

    # Wrongly classified as irreversible (actually reversible)
    with pytest.raises(IrrevElemReactionError):
        test = IrrevElemReaction(rxn_type="Elementary",
                                 is_reversible=True, # This can't happen for irrev rxns
                                 rxn_equation="H2 + OH =] H2O + H",
                                 species_list=['H', 'O', 'OH', 'H2', 'H2O', 'O2'],
                                 rate_coeffs_components={'k': 10},
                                 rate_coeffs_type='constant',
                                 reactant_stoich_coeffs={'H2' :1, 'OH':1},
                                 product_stoich_coeffs={'H2O' :1, 'H':1})

    # Wrongly classified as elementary (actually non-elementary)
    with pytest.raises(IrrevElemReactionError):
        test = IrrevElemReaction(rxn_type="Non-elementary", # This can't happen for elem rxns
                                 is_reversible=False,
                                 rxn_equation="H2 + OH =] H2O + H",
                                 species_list=['H', 'O', 'OH', 'H2', 'H2O', 'O2'],
                                 rate_coeffs_components={'k': 10},
                                 rate_coeffs_type='constant',
                                 reactant_stoich_coeffs={'H2' :1, 'OH':1},
                                 product_stoich_coeffs={'H2O' :1, 'H':1})

def test_irrev_elem_reaction_functions(test_irrev_reaction):
    """Tests for irreversible elementary reaction functions."""

    # Test computing progress rate without setting concentrations
    with pytest.raises(ValueError):
        test_irrev_reaction.compute_progress_rate()

    # Test computing reaction rate coefficient
    k = test_irrev_reaction.compute_reaction_rate_coeff()
    assert k == 10

    # Test invalid reaction rate coefficient
    test_irrev_reaction.rate_coeffs_components = {'k': -10}
    with pytest.raises(ValueError):
        test_irrev_reaction.compute_reaction_rate_coeff()
    test_irrev_reaction.rate_coeffs_components = {'k': 10} # change back

    # Test computing progress rate
    test_irrev_reaction.set_concentrations({'A': 1, 'B': 2, 'C': 3})
    w = test_irrev_reaction.compute_progress_rate()
    expected = numpy.array([ 20.])
    assert w == expected

    # Test computing reaction rate
    rxnrate = test_irrev_reaction.compute_reaction_rate()
    expected = numpy.array([-40.0, -20.0, 20.0])
    assert (rxnrate == expected).all()

    # Test computing reaction rate when negative reactant stoichiometric coefficients
    test_irrev_reaction.reactant_stoich_coeffs = {'A': 2, 'B': -1}
    with pytest.raises(ValueError):
        rxnrate = test_irrev_reaction.compute_reaction_rate()
    test_irrev_reaction.reactant_stoich_coeffs= {'A': 2, 'B': 1} # change back

    # Test computing reaction rate when negative product stoichiometric coefficients
    test_irrev_reaction.product_stoich_coeffs = {'C': -1}
    with pytest.raises(ValueError):
        rxnrate = test_irrev_reaction.compute_reaction_rate()

def test_irrev_elem_reaction_arrhenius(test_irrev_reaction_arrhenius):
    """Test irrev elem reactions with arrhenius rxn rate coefficient types."""

    # When T is not inputed
    try:
        k = test_irrev_reaction_arrhenius.compute_reaction_rate_coeff() 
    except TypeError:
        test_irrev_reaction_arrhenius.set_temperature(10)
        k = test_irrev_reaction_arrhenius.compute_reaction_rate_coeff()
    assert numpy.isclose(k, 3.0037488791204288, atol=1e-16)

def test_irrev_elem_reaction_modarrhenius(test_irrev_reaction_modified_arr):
    """Test irrev elem reactions with modified arrhenius rxn rate coefficient types."""

    # When T not inputed by user
    try:
        k = test_irrev_reaction_modified_arr.compute_reaction_rate_coeff() 
    except TypeError:
        test_irrev_reaction_modified_arr.set_temperature(10)
        k = test_irrev_reaction_modified_arr.compute_reaction_rate_coeff()
    assert numpy.isclose(k, 9.498687977198342, atol=1e-16)


# ===== Tests for reversible elementary reactions ===== #

@pytest.fixture
def test_rev_reaction():
    """Returns a valid reaction (from rev_rxn.xml)"""
    return RevElemReaction(rxn_type="Elementary",
                    is_reversible=True,
                    rxn_equation="H + O2 [=] H2O ",
                    species_list=['H', 'H2O', 'O2'],
                    rate_coeffs_components={'k': 10},
                    rate_coeffs_type='constant',
                    reactant_stoich_coeffs={'H' :2, 'O2':1},
                    product_stoich_coeffs={'H2O' :2})

def test_wrongly_classified_rev_elem_rxns():
    """Tests for wrongly classified reversible elementary reactions."""

    # Wrongly classified as reversible (actually irreversible)
    with pytest.raises(RevElemReactionError):
        test = RevElemReaction(rxn_type="Elementary",
                                is_reversible=False, # This can't happen for rev rxns
                                rxn_equation="H2 + OH =] H2O + H",
                                species_list=['H', 'O', 'OH', 'H2', 'H2O', 'O2'],
                                rate_coeffs_components={'k': 10},
                                rate_coeffs_type='constant',
                                reactant_stoich_coeffs={'H2' :1, 'OH':1},
                                product_stoich_coeffs={'H2O' :1, 'H':1})

    # Wrongly classified as elementary (actually non-elementary)
    with pytest.raises(RevElemReactionError):
        test = RevElemReaction(rxn_type="Non-elementary", # This can't happen for elem rxns
                                    is_reversible=True,
                                    rxn_equation="H2 + OH =] H2O + H",
                                    species_list=['H', 'O', 'OH', 'H2', 'H2O', 'O2'],
                                    rate_coeffs_components={'k': 10},
                                    rate_coeffs_type='constant',
                                    reactant_stoich_coeffs={'H2' :1, 'OH':1},
                                    product_stoich_coeffs={'H2O' :1, 'H':1})

def test_rev_elem_reaction_functions(test_rev_reaction):
    """Tests for reversible elementary reaction functions.""" 

    # Test computing NASA7 reaction rate coefficient without setting NASA polynomials
    with pytest.raises(ValueError):
        test_rev_reaction.compute_reaction_rate_coeff()

    # Test computing NASA 7 reaction rate coefficient
    T = 400
    test_rev_reaction.set_temperature(T)
    lowT_nasa = {'O2': numpy.array([3.21293640e+00, 1.12748635e-03,
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
    test_rev_reaction.set_NASA_poly_coefs(lowT_nasa)
    kf, kb = test_rev_reaction.compute_reaction_rate_coeff(T=T)
    assert kf == 10
    assert numpy.isclose(kb, 0, atol=1e-16)

    # Test computing invalid type reaction rate coefficient
    test_rev_reaction.bkwd_coeff_type = "unsupported"
    with pytest.raises(NotImplementedError):
        test_rev_reaction.compute_reaction_rate_coeff(T=T)
    test_rev_reaction.bkwd_coeff_type = "NASA7" # change back

    # Test computing progress rate without setting concentrations
    with pytest.raises(ValueError):
        prog_rate = test_rev_reaction.compute_progress_rate(T=T)

    # Test computing progress rate
    test_rev_reaction.set_concentrations(X={'H':1, 'O2':1, 'H2O':1})
    prog_rate = test_rev_reaction.compute_progress_rate(T=T)
    assert prog_rate == 10
<<<<<<< HEAD

=======
>>>>>>> origin/develop
