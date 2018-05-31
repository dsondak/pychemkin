
"""Tests for rxn_rate_coefficients module"""

import numpy
import pytest
import os, sys
#sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import warnings
warnings.simplefilter("error")

from pychemkin.rxn_rate_coefficients.rxn_rate_coefficients import *


# ===== Tests for forward reaction rate coefficients ===== #

def test_unhandled_rxn_rate_coeff_type():
    """Test invalid/unhandled rxn rate coefficient type."""
    rate_coeffs_type = "madeup type"
    with pytest.raises(KeyError):
        test = determine_rxn_rate_coeff_type(rate_coeffs_type)

def test_invalid_components():
    """Test reaction rate coefficients when
    invalid/unhandled component in dictionary."""

    # constant
    k_parameters = {'sdfs': 10}
    with pytest.raises(ValueError):
        k_test = ConstantFwdCoeff(k_parameters)

    # Arrhenius
    k_parameters = {'sdfdds': 10, 'E': 10**3}
    T = 10
    with pytest.raises(ValueError):
        k_test = ArrheniusFwdCoeff(k_parameters, T)

    # modified Arrhenius
    k_parameters = {'sdfdds': 10, 'E': 10**3, 'b': 10}
    T = 10
    with pytest.raises(ValueError):
        k_test = ArrheniusFwdCoeff(k_parameters, T)

def test_invalid_or_extra_components():
    """Test reaction rate coefficients when
    extra component in dictionary."""

    # constant
    k_parameters = {'k': 10, 'q':8}
    with pytest.warns(UserWarning):
        k_test = ConstantFwdCoeff(k_parameters).k

    # Arrhenius
    k_parameters = {'A': 10**7, 'E':10**3, 'R': 8.3144598, 'madeup':120}
    T = 10
    with pytest.warns(UserWarning):
        k_test = ArrheniusFwdCoeff(k_parameters, T).k
<<<<<<< HEAD

    k_parameters = {'A': 10**7, 'E':10**3, 'madeup':120}
    T = 10
    with pytest.warns(UserWarning):
        k_test = ArrheniusFwdCoeff(k_parameters, T).k

    # modified Arrhenius
    k_parameters = {'A': 10**7, 'E':10**3, 'b':0.5, 'madeup':120, 'R': 8.3144598}
    T = 10
    with pytest.warns(UserWarning):
        k_test = ModifiedArrheniusFwdCoeff(k_parameters, T).k

    k_parameters = {'A': 10**7, 'E':10**3, 'b':0.5, 'madeup':120}
    T = 10
    with pytest.warns(UserWarning):
        k_test = ModifiedArrheniusFwdCoeff(k_parameters, T).k

def test_constant_rxn_rate_coefficient():
    """Tests for constant reaction rate coefficient."""

=======

    k_parameters = {'A': 10**7, 'E':10**3, 'madeup':120}
    T = 10
    with pytest.warns(UserWarning):
        k_test = ArrheniusFwdCoeff(k_parameters, T).k

    # modified Arrhenius
    k_parameters = {'A': 10**7, 'E':10**3, 'b':0.5, 'madeup':120, 'R': 8.3144598}
    T = 10
    with pytest.warns(UserWarning):
        k_test = ModifiedArrheniusFwdCoeff(k_parameters, T).k

    k_parameters = {'A': 10**7, 'E':10**3, 'b':0.5, 'madeup':120}
    T = 10
    with pytest.warns(UserWarning):
        k_test = ModifiedArrheniusFwdCoeff(k_parameters, T).k

def test_constant_rxn_rate_coefficient():
    """Tests for constant reaction rate coefficient."""

>>>>>>> origin/develop
    # Test when invalid constant (non-positive k)
    k_parameters = {'k': -10}
    with pytest.raises(ValueError):
        k_test = ConstantFwdCoeff(k_parameters).k
<<<<<<< HEAD

    # Compute with valid input
    k_parameters = {'k': 10}
    k_test = ConstantFwdCoeff(k_parameters).k
    assert k_test == 10

=======

    # Compute with valid input
    k_parameters = {'k': 10}
    k_test = ConstantFwdCoeff(k_parameters).k
    assert k_test == 10

>>>>>>> origin/develop
    # Test when T entered (no effect)
    k_parameters = {'k': 10}
    T = 10
    k_test = ConstantFwdCoeff(k_parameters, T).k
    assert k_test == 10

def test_arrhenius_rxn_rate_coefficient():
    """Tests for Arrhenius reaction rate coefficient."""

    # Test when missing argument T
    k_parameters = {'A': 10, 'E':100, 'R':8.3144598}
    with pytest.raises(TypeError):
        k_test = ArrheniusFwdCoeff(k_parameters).k

    # Compute with valid inputs
    k_parameters = {'A': 10**7, 'E': 10**3, 'R': 8.3144598}
    T = 10**2
    k_test = ArrheniusFwdCoeff(k_parameters, T).k
    assert numpy.isclose(k_test, 3003748.8791204286, atol=1e-16)

    # Test when invalid value for T (non-positive temperature)
    k_parameters = {'A': 10, 'E': 100, 'R': 8.3144598}
    T = -10
    with pytest.raises(ValueError):
        k_test = ArrheniusFwdCoeff(k_parameters, T).k

    # Test when invalid value for A (non-positive prefactor)
    k_parameters = {'A': 0, 'E': 100, 'R': 8.3144598}
    T = 10
    with pytest.raises(ValueError):
        k_test = ArrheniusFwdCoeff(k_parameters, T).k

    # Test when invalid value for R (non-positive gas constant)
    k_parameters = {'A': 10, 'E': 100, 'R': -100}
    T = 10
    with pytest.raises(ValueError):
        k_test = ArrheniusFwdCoeff(k_parameters, T).k

    # Test when changing value of R
    k_parameters = {'A': 10, 'E': 100, 'R': 10.453}
    T = 10
    with pytest.warns(UserWarning):
        k_test = ArrheniusFwdCoeff(k_parameters, T).k

def test_modified_arrhenius_rxn_rate_coefficient():
    """Tests for Modified Arrhenius reaction rate coefficient."""

    # Test when missing argument T
    k_parameters = {'A': 10, 'E':100, 'b':0.5, 'R':8.3144598}
    with pytest.raises(TypeError):
        k_test = ModifiedArrheniusFwdCoeff(k_parameters).k

    # Compute with valid inputs
    k_parameters = {'A': 10**7, 'E':10**3, 'b':0.5, 'R':8.3144598}
    T = 10**2
    k_test = ModifiedArrheniusFwdCoeff(k_parameters, T).k
    assert numpy.isclose(k_test, 30037488.791204285, atol=1e-16)

    # Test when invalid value for A (non-positive prefactor)
    k_parameters = {'A': -10, 'E':100, 'b':0.5, 'R':8.3144598}
    T = 10
    with pytest.raises(ValueError):
        k_test = ModifiedArrheniusFwdCoeff(k_parameters, T).k

    # Test when invalid value for b (non-real constant)
    k_parameters = {'A': 10, 'E':100, 'b':0.5j, 'R':8.3144598}
    T = 10
    with pytest.raises(TypeError):
        k_test = ModifiedArrheniusFwdCoeff(k_parameters, T).k

    # Test when invalid value for T (non-positive temperature)
    k_parameters = {'A': 10, 'E':100, 'b':0.5, 'R':8.3144598}
    T = -10
    with pytest.raises(ValueError):
        k_test = ModifiedArrheniusFwdCoeff(k_parameters, T).k

    # Test when invalid value for R (non-positive gas constant)
    k_parameters = {'A': 10, 'E':100, 'R':-100, 'b':0.5}
    T = 10
    with pytest.raises(ValueError):
        k_test = ModifiedArrheniusFwdCoeff(k_parameters, T).k

    # Test when changing value of R
    k_parameters = {'A': 10, 'E':100, 'R':10.453, 'b':0.5}
    T = 10
    with pytest.warns(UserWarning):
        k_test = ModifiedArrheniusFwdCoeff(k_parameters, T).k


# ===== Tests for backward reaction rate coefficients ===== #

def test_backward_coefficient_base_class():
    """Test BackwardReactionCoeff class"""
    bkwd_coeff = BackwardReactionCoeff()
    with pytest.raises(NotImplementedError):
        bkwd_coeff.compute_bkwd_coefficient()

@pytest.fixture
def NASA7_backward_coeff_setup():
    """Returns a working (but artificial) example of
    backward reaction rate coefficient using NASA
    7 polynomial coefficients."""
    expected_nasa = numpy.array([[1,0,0,0,0,0,0],
                                [1,0,0,0,0,0,0],
                                [1,0,0,0,0,0,0]])
    k_f = 100
    nu_i = numpy.array([-2, -1, 2])
    bkwd_coeff = NASA7BackwardCoeff(nu_i, expected_nasa,
                                     p0=1e5, R=8.3144598)
    return bkwd_coeff

def test_NASA7_backward_coeff(NASA7_backward_coeff_setup):
    """Test backward rxn rate coefficient using NASA
    7 polynomial coefficients."""

    # Test when changing value of R
    expected_nasa = numpy.array([[1,0,0,0,0,0,0],
                                [1,0,0,0,0,0,0],
                                [1,0,0,0,0,0,0]])
    k_f = 100
    nu_i = numpy.array([-2, -1, 2])
    with pytest.warns(UserWarning):
        kwd_coeff = NASA7BackwardCoeff(nu_i, expected_nasa,
                                        p0=1e5, R=43.3423)

    # Test value of gamma
    assert NASA7_backward_coeff_setup.gamma == -1

    # Test computation of H/RT
    T = 100
    expected_H_over_RT = numpy.array([1, 1, 1])
    assert numpy.isclose(NASA7_backward_coeff_setup._H_over_RT(T),
                         expected_H_over_RT, atol=1e-16).all()

    # Test computation of H/RT with invalid T (non-positive temperature)
    T = -100
    with pytest.raises(ValueError):
        NASA7_backward_coeff_setup._H_over_RT(T)

    # Test computation of S/R
    T = 100
    expected_S_over_R = numpy.array([4.60517, 4.60517, 4.60517])
    assert numpy.isclose(NASA7_backward_coeff_setup._S_over_R(T),
                         expected_S_over_R, atol=1e-16).all()

    # Test computation of S/R with invalid T (non-positive temperature)
    T = -100
    with pytest.raises(ValueError):
        NASA7_backward_coeff_setup._S_over_R(T)

    # Test computation of backward rxn rate coefficient k_b
    T = 100
    k_f = 100
    expected_delta_S_over_R = -4.60517
    expected_delta_H_over_RT = -1

    fact =  NASA7_backward_coeff_setup.p0 / NASA7_backward_coeff_setup.R / T
    expected_gamma = -1
    expected_ke = (fact ** expected_gamma) * (numpy.exp(expected_delta_S_over_R - 
                                                        expected_delta_H_over_RT))

    expected_kb_val = 442457 # 100 / 2.260104919e-6
    assert numpy.isclose(NASA7_backward_coeff_setup.compute_bkwd_coefficient(k_f, T),
                         expected_kb_val, atol=1e-16)

    # Test computation of k_b with invalid T (non-positive temperature)
    T = -100
    k_f = 100
    with pytest.raises(ValueError):
        NASA7_backward_coeff_setup.compute_bkwd_coefficient(k_f, T)
