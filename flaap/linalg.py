"""Linear-algebra helpers ported from the Fortran utility routines.

=====================  =======================================================
Fortran source         Python replacement
=====================  =======================================================
c_inverse.f90          :func:`cinvert_lu` (module CLU / CLU_D, LU inverse)
C_det1.f90             :func:`cdet` (TSCGT/DCGT full-pivot determinant)
c_eigen.f90            :func:`comeig` (COMEIG complex eigen-solver)
c_matrix_sort_imag     :func:`sort_imag`, :func:`sort_real` (PIKSRT)
transpose_vector.f90   :func:`transpose_vector`
Thin_plate_spline.f90  BLAS/LINPACK parts -> numpy; RBF parts -> flaap.structures.spline
TPS_lib.f90            J. Burkardt R8LIB utility library -> numpy/scipy equivalents
=====================  =======================================================

The originals were hand-written numerical recipes; here they are thin,
well-tested wrappers around LAPACK (numpy/scipy) with the same semantics
(including the error code convention of ``CINVERT_LU``).
"""
from __future__ import annotations

import numpy as np
import scipy.linalg as sla


def cinvert_lu(a: np.ndarray) -> tuple[np.ndarray, int]:
    """Inverse via LU decomposition.  Returns ``(inverse, errorcode)``.

    ``errorcode`` is 0 on success and 1 for a singular matrix, as in the
    Fortran ``CINVERT_LU`` / ``CINVERT_LU_D`` routines.
    """
    a = np.asarray(a)
    try:
        lu, piv = sla.lu_factor(a, check_finite=True)
        if np.any(np.abs(np.diag(lu)) == 0.0):
            return np.full_like(a, np.nan), 1
        inv = sla.lu_solve((lu, piv), np.eye(a.shape[0], dtype=lu.dtype))
        return inv, 0
    except (sla.LinAlgError, ValueError):
        return np.full_like(a, np.nan), 1


def cdet(a: np.ndarray, eps: float = 1e-300) -> complex:
    """Determinant of a (complex) square matrix with full pivoting (DCGT)."""
    a = np.asarray(a, dtype=complex)
    if a.size == 0:
        return 1.0 + 0j
    det = np.linalg.det(a)
    return 0j if abs(det) < eps else complex(det)


def comeig(a: np.ndarray, vectors: bool = True):
    """Eigenvalues (and right eigenvectors) of a complex matrix (COMEIG)."""
    a = np.asarray(a, dtype=complex)
    if vectors:
        w, v = np.linalg.eig(a)
        return w, v
    return np.linalg.eigvals(a)


def sort_imag(a: np.ndarray) -> np.ndarray:
    """Sort complex values by ascending imaginary part (PIKSRT_imag)."""
    a = np.asarray(a)
    return a[np.argsort(a.imag, kind="stable")]


def sort_real(a: np.ndarray) -> np.ndarray:
    """Sort complex values by ascending real part (PIKSRT_real)."""
    a = np.asarray(a)
    return a[np.argsort(a.real, kind="stable")]


def transpose_vector(a: np.ndarray) -> np.ndarray:
    """(1,N) row -> (N,1) column."""
    return np.asarray(a).reshape(-1, 1)


def mac(phi1: np.ndarray, phi2: np.ndarray) -> np.ndarray:
    """Modal assurance criterion matrix between column sets (complex ok)."""
    phi1 = np.atleast_2d(phi1.T).T if phi1.ndim == 1 else phi1
    phi2 = np.atleast_2d(phi2.T).T if phi2.ndim == 1 else phi2
    num = np.abs(phi1.conj().T @ phi2) ** 2
    d1 = np.real(np.sum(phi1.conj() * phi1, axis=0))
    d2 = np.real(np.sum(phi2.conj() * phi2, axis=0))
    return num / np.outer(d1, d2)
