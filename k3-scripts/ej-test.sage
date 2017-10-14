def ej_test(pol):
    """
    Test a polynomial for the Elsenhans-Jahnel condition.
    """
    polRing = pol.parent()
    x = polRing.gen()

    pol1 = pol[0].sign() * pol // (1-x)^pol.ord(1-x)
    return(pol1(-1).is_square())
