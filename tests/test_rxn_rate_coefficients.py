
"""Tests for rxn_rate_coefficients module"""

import numpy
import pytest

from pychemkin.rxn_rate_coefficients.rxn_rate_coefficients import (RxnCoefficientBase,
                                                                   ConstantCoefficient,
                                                                   ArrheniusCoefficient,
                                                                   ModifiedArrheniusCoefficient)


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
