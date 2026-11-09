"""Catch absolute probability guards and false acceptance at float ties."""
import math
from decimal import Decimal, localcontext
from fractions import Fraction

import pytest

from panel_mic.numerics import tail_target, _power_ge


def independent_cdf(n, k, q):
    q = Fraction(q)
    return sum((math.comb(n, j) * q**j * (1-q)**(n-j)
                for j in range(k+1)), Fraction(0))


def test_tiny_target_retains_a_positive_safe_inverse():
    target = 1e-100
    q = tail_target(1000, target_power=target)
    assert .2 < q < .3  # Independent bracketing of the k=3 binomial polynomial.
    assert independent_cdf(1000, 3, q) >= Fraction(target)
    with localcontext() as ctx:
        ctx.prec = 100
        # A practical inverse, not a fallback close to zero.
        probability = Decimal(independent_cdf(1000, 3, q).numerator) / Decimal(independent_cdf(1000, 3, q).denominator)
        assert probability / Decimal.from_float(target) < Decimal('1.000000001')


def test_exact_subnormal_power_equality_is_accepted():
    target = math.nextafter(0., 1.)
    # 2**(-1074) is exactly the smallest positive binary64 number.
    assert _power_ge(1074, 0, .5, target)


@pytest.mark.parametrize('n,k,q', [(8, 2, .25), (12, 3, .5), (15, 1, .125)])
def test_power_decision_matches_exact_rational_either_side(n, k, q):
    exact = independent_cdf(n, k, q)
    center = float(exact)
    for target in (math.nextafter(center, 0.), center, math.nextafter(center, 1.)):
        assert _power_ge(n, k, q, target) is (exact >= Fraction(target))
