"""Classical unsteady-aerodynamics functions against literature values."""
import numpy as np
import pytest
from scipy.special import hankel2

from flaap.aero.strip import sears, theodorsen


@pytest.mark.parametrize("k, F, G", [
    # Theodorsen function C(k) = F + iG, tabulated in Bisplinghoff, Ashley &
    # Halfman, Aeroelasticity (1955) Table 5-1
    (0.10, 0.8319, -0.1723),
    (0.50, 0.5979, -0.1507),
    (1.00, 0.5394, -0.1003),
])
def test_theodorsen_table(k, F, G):
    C = complex(theodorsen(k))
    assert C.real == pytest.approx(F, abs=2e-4)
    assert C.imag == pytest.approx(G, abs=2e-4)


@pytest.mark.parametrize("k, ref", [
    # Sears function referred to the mid-chord, closed form
    # S(k) = 2 / (pi k [H0^(2)(k) - i H1^(2)(k)])  (Sears 1941)
    (0.1, 0.82124 - 0.16348j),
    (1.0, 0.36865 + 0.12594j),
])
def test_sears_closed_form(k, ref):
    S = complex(sears(k))
    closed = 2.0 / (np.pi * k * (hankel2(0, k) - 1j * hankel2(1, k)))
    assert S == pytest.approx(closed, abs=1e-6)
    assert S == pytest.approx(ref, abs=1e-4)


def test_quasi_steady_limits():
    assert complex(theodorsen(0.0)) == 1.0
    assert complex(sears(0.0)) == 1.0
    # C(k) -> 1/2 for k -> infinity
    assert abs(complex(theodorsen(200.0)) - 0.5) < 2e-3
