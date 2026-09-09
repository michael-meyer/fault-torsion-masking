from sage.all import Matrix, ZZ

# ============================== #
#  Utility functions for masking #
# ============================== #


def canonical_matrix(M):
    """
    We recover M up to a minus sign, we so
    ensure the sign of all matrices are
    canonical
    """
    alpha = M[0][0]
    if alpha < -alpha:
        return -M
    return M


def random_diag_matrix(R, alpha=None):
    r"""
    Compute a random, diagonal, unitary
    invertible matrix M \in Z / (2^b) Z
    """
    if alpha is None:
        alpha = R.random_element()

    # Ensure alpha is an odd unit modulo 2^b
    alpha = alpha + ZZ(alpha % 2) + 1
    M = Matrix(R, 2, 2, [alpha, 0, 0, ~alpha])

    return canonical_matrix(M)


def random_circulant_matrix(R, beta=None):
    r"""
    Compute a random, circulant, unitary
    invertible matrix M \in Z / (2^b) Z
    """
    if beta is None:
        beta = R.random_element()
    b = 4 * beta
    aa = R(b * b + 1)
    a = aa.sqrt()  # TODO: fast sqrt mod 2^k?
    M = Matrix(R, 2, 2, [a, b, b, a])

    return canonical_matrix(M)


def random_matrix(R, ele=None, diag=True):
    """
    Compute a random masking matrix which
    is commutative, invertible and unitary
    in Z / (N Z)
    """
    if diag:
        return random_diag_matrix(R, alpha=ele)
    return random_circulant_matrix(R, beta=ele)


def _bit_flip_scalar(a, bit_index=1):
    """Flip one bit of a masking coefficient."""
    bit_index = int(bit_index)
    if bit_index < 0:
        raise ValueError("bit_index must be non-negative")
    R = a.parent()
    return R(int(ZZ(a)) ^ (1 << bit_index))


def mask_torsion_points(A, P, Q, fault_coeff=None, fault_bit=1):
    """
    Evaluate the masking matrix A on [P,Q].
    If fault_coeff is set, simulate the bit-flip fault in that multiplication.
    """
    a00, a01, a10, a11 = A.list()

    if fault_coeff == "a00":
        a00 = _bit_flip_scalar(a00, fault_bit)
    elif fault_coeff == "a01":
        a01 = _bit_flip_scalar(a01, fault_bit)
    elif fault_coeff == "a10":
        a10 = _bit_flip_scalar(a10, fault_bit)
    elif fault_coeff == "a11":
        a11 = _bit_flip_scalar(a11, fault_bit)
    elif fault_coeff is not None:
        raise ValueError("unknown fault coefficient: %s" % fault_coeff)

    R = a00 * P + a01 * Q
    S = a10 * P + a11 * Q
    return R, S
