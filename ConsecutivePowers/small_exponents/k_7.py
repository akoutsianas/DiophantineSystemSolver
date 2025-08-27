from sage.all import ZZ, QQ, polygen, EllipticCurve, NumberField
from sage.misc.functional import sqrt

x = polygen(QQ, 'x')

# Case n = 2
for t in range(2):
    A = 2**6 * 7**t
    E = EllipticCurve([0, A*21, 0, A**2 * 35, A**3 * 7])
    for pt in E.integral_points(both_signs=True):
        if pt[0] % A == 0:
            x0 = pt[0] / A
            if x0.is_square():
                x0 = ZZ(sqrt(x0))
                print(f"Value of x = {x0}")


# Case n = 3
K = NumberField(x^3 + 21*x^2 + 35*x + 7, 'rk')
rk = K.gen()
S = K.primes_above(14)
for u in K.selmer_group_iterator(S, 3):
    Eu = EllipticCurve([0, rk*u^2])
    try:
        r = Eu.rank()
        print(f"r={r} for u={u}")
        gens = Eu.gens()
        print(f"gens={gens}")
    except:
        print(f"Not possible to compute the rank for u={u}")


