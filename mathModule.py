from scipy.special import comb, gamma
from math import factorial
from sys import exit


def binomial(term1, term2):
    return float(comb(term1, term2, exact=True))
    # Note: scipy.special.binom is approximate, and scipy.special.comb(a,b,exact=True) is exact


def fact(var):
    if int(var) != var:
        print("Not setup for factoials of decimals in fact function.\n...............\n...............")
        exit()
    elif var < 0:
        return 1.0
    else:
        return float(factorial(var))


def Gamma(var):
    return float(gamma(var))
