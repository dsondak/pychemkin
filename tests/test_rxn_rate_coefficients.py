
"""Tests for rxn_rate_coefficients module"""

import numpy
import pytest

import warnings
# Treat warnings like errors (for testing purposes)
warnings.simplefilter("error")

from pychemkin.rxn_rate_coefficients.rxn_rate_coefficients import (RxnCoefficientBase,
                                                                   ConstantCoefficient,
                                                                   ArrheniusCoefficient,
                                                                   ModifiedArrheniusCoefficient,
                                                                   BackwardRxnCoefficientBase,
                                                                   NASA7BackwardCoeffs)


def test_instantiating_RxnCoefficientBase():
    """Test creating instance of RxnCoefficientBase."""
    k_base = RxnCoefficientBase()
    assert k_base.k is None
    assert repr(k_base) == "RxnCoefficientBase()"
    with pytest.raises(NotImplementedError):
        k_base.get_coefficient()

def test_instantiating_ConstantCoefficient():
    """Test creating instance of ConstantCoefficient."""
    test_val = 10.0
    k_const = ConstantCoefficient(k=test_val)
    assert k_const.k == 10.0
    assert repr(k_const) == "ConstantCoefficient(k=10.0)"

def test_invalid_ConstantCoefficient():
    """Test non-positive constant values."""
    test_val = 0.0 # zero
    with pytest.raises(ValueError):
        k_const = ConstantCoefficient(k=test_val)

    test_val = -10.0 # negative
    with pytest.raises(ValueError):
        k_const = ConstantCoefficient(k=test_val)

def test_instantiating_ArrheniusCoefficient():
    """Test creating instance of ArrheniusCoefficient."""
    A = 10.0
    E = 1.0
    T = 10.0
    k_arr = ArrheniusCoefficient(A=A, E=E, T=T)
    k_arr.get_coefficient()
    assert numpy.isclose(k_arr.k, 9.8804414136336476)
    assert repr(k_arr) == "ArrheniusCoefficient(A=10.0, E=1.0, T=10.0, R=8.314)"

def test_invalid_A_ArrheniusCoefficient():
    """Test invalid inputs for ArrheniusCoefficient component A."""
    A = -10.0 # negative
    E = 1.0
    T = 10.0
    with pytest.raises(ValueError):
        k_arr = ArrheniusCoefficient(A=A, E=E, T=T)

    A = 0.0 # zero
    with pytest.raises(ValueError):
        k_arr = ArrheniusCoefficient(A=A, E=E, T=T)

def test_invalid_T_ArrheniusCoefficient():
    """Test invalid inputs for ArrheniusCoefficient component T."""
    A = 10.0
    E = 1.0
    T = -10.0 # negative
    with pytest.raises(ValueError):
        k_arr = ArrheniusCoefficient(A=A, E=E, T=T)

    T = 0.0 # zero
    with pytest.raises(ValueError):
        k_arr = ArrheniusCoefficient(A=A, E=E, T=T)

def test_invalid_R_ArrheniusCoefficient():
    """Test invalid inputs for ArrheniusCoefficient component R."""
    A = 10.0
    E = 1.0
    T = 10.0
    R = -10.0 # negative
    with pytest.raises(ValueError):
        k_arr = ArrheniusCoefficient(A=A, E=E, T=T, R=R)

    R = 0.0 # zero
    with pytest.raises(ValueError):
        k_arr = ArrheniusCoefficient(A=A, E=E, T=T, R=R)

def test_overflow_ArrheniusCoefficient():
    """Test too large of an output for ArrheniusCoefficient."""
    A = 10**1000
    E = 1.0
    T = 10.0
    with pytest.raises(OverflowError):
        k_arr = ArrheniusCoefficient(A=A, E=E, T=T)
        k_arr.get_coefficient()

def test_instantiating_ModifiedArrheniusCoefficient():
    """Test creating instance of ModifiedArrheniusCoefficient."""
    A = 10.0
    b = 12.0
    E = 1.0
    T = 10.0
    k_modarr = ModifiedArrheniusCoefficient(A=A, b=b, E=E, T=T)
    k_modarr.get_coefficient()
    assert numpy.isclose(k_modarr.k, 9880441413633.648)
    assert repr(k_modarr) == "ModifiedArrheniusCoefficient(A=10.0, b=12.0, E=1.0, T=10.0, R=8.314)"

def test_invalid_A_ModifiedArrheniusCoefficient():
    """Test invalid inputs for ModifiedArrheniusCoefficient component A."""
    A = -10.0 # negative
    b = 10.0
    E = 1.0
    T = 10.0
    with pytest.raises(ValueError):
        k_modarr = ModifiedArrheniusCoefficient(A=A, b=b, E=E, T=T)

    A = 0.0 # zero
    with pytest.raises(ValueError):
        k_modarr = ModifiedArrheniusCoefficient(A=A, b=b, E=E, T=T)

def test_invalid_T_ModifiedArrheniusCoefficient():
    """Test invalid inputs for ModifiedArrheniusCoefficient component T."""
    A = 10.0
    b = 10.0
    E = 1.0
    T = -10.0 # negative
    with pytest.raises(ValueError):
        k_modarr = ModifiedArrheniusCoefficient(A=A, b=b, E=E, T=T)

    T = 0.0 # zero
    with pytest.raises(ValueError):
        k_modarr = ModifiedArrheniusCoefficient(A=A, b=b, E=E, T=T)

def test_invalid_R_ModifiedArrheniusCoefficient():
    """Test invalid inputs for ModifiedArrheniusCoefficient component R."""
    A = 10.0
    b = 10.0
    E = 1.0
    T = 10.0
    R = -10.0 # negative
    with pytest.raises(ValueError):
        k_modarr = ModifiedArrheniusCoefficient(A=A, b=b, E=E, T=T, R=R)

    R = 0.0 # zero
    with pytest.raises(ValueError):
        k_modarr = ModifiedArrheniusCoefficient(A=A, b=b, E=E, T=T, R=R)

def test_invalid_b_ModifiedArrheniusCoefficient():
    """Test invalid input for ModifiedArrheniusCoefficient component b."""
    A = 10.0
    b = 1.j
    E = 1.0
    T = 10.0
    with pytest.raises(ValueError):
        k_modarr = ModifiedArrheniusCoefficient(A=A, b=b, E=E, T=T)

def test_overflow_ModifiedArrheniusCoefficient():
    """Test too large of an output for ModifiedArrheniusCoefficient."""
    A = 10**1000
    b = 10.0
    E = 1.0
    T = 10.0
    with pytest.raises(OverflowError):
        k_modarr = ModifiedArrheniusCoefficient(A=A, b=b, E=E, T=T)
        k_modarr.get_coefficient()



# ======================= TESTS FOR BACKWARDCOEFF ====================== #

@pytest.fixture
def test_backward_coeff():
    """Returns a working (but artificial) example of
    backward reaction rate coefficient."""
    expected_nasa = numpy.array([[1,0,0,0,0,0,0],
                                [1,0,0,0,0,0,0],
                                [1,0,0,0,0,0,0]])
    k_f = 100
    nu_i = numpy.array([-2, -1, 2])
    bkwd_coeff = NASA7BackwardCoeffs(nu_i, expected_nasa)
    return bkwd_coeff

def test_backwardCoeff_gamma(test_backward_coeff):
    """Tests value of gamma for working example."""
    assert test_backward_coeff.gamma == -1

def test_backwardCoeff_computing_H(test_backward_coeff):
    """Tests computing H/RT for working example."""
    T = 100
    expected_H_over_RT = numpy.array([1, 1, 1])
    assert numpy.isclose(test_backward_coeff.H_over_RT(T),
                         expected_H_over_RT).all()

def test_backwardCoeff_computing_H_neg_T(test_backward_coeff):
    """Tests computing H/RT for working example with neg T."""
    T = -100
    with pytest.raises(ValueError):
        test_backward_coeff.H_over_RT(T)

def test_backwardCoeff_computing_S(test_backward_coeff):
    """Tests computing S/R for working example."""
    T = 100
    expected_S_over_R = numpy.array([4.60517, 4.60517, 4.60517])
    assert numpy.isclose(test_backward_coeff.S_over_R(T),
                         expected_S_over_R).all()

def test_backwardCoeff_computing_S_neg_T(test_backward_coeff):
    """Tests computing S/R for working example with neg T."""
    T = -100
    with pytest.raises(ValueError):
        test_backward_coeff.S_over_R(T)

def test_backwardCoeff_computeCoeff(test_backward_coeff):
    """Tests computing k_b for working example."""
    T = 100
    k_f = 100
    expected_delta_S_over_R = -4.60517
    expected_delta_H_over_RT = -1

    fact =  test_backward_coeff.p0 / test_backward_coeff.R / T
    expected_gamma = -1
    expected_ke = (fact ** expected_gamma) * (numpy.exp(expected_delta_S_over_R - 
                                                        expected_delta_H_over_RT))

    expected_kb_val = 442457 # 100 / 2.260104919e-6

    assert numpy.isclose(test_backward_coeff.compute_backward_coeffs(k_f, T),
                         expected_kb_val)

def test_backwardCoeff_computeCoeff_neg_T(test_backward_coeff):
    """Tests computing k_b for working example with neg T."""
    T = -100
    k_f = 100
    with pytest.raises(ValueError):
        test_backward_coeff.compute_backward_coeffs(k_f, T)
