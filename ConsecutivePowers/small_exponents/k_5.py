from sage.all import ZZ, QQ, polygen, EllipticCurve
from sage.misc.functional import sqrt

x = polygen(QQ, 'x')

# Case n = 2
# Use Niko's code for integral points on genus one curves or we use Chabauty method for
# y^2 = x(x^4 + 10x^2 + 5) as we compute that x = 0, 1


# Case n = 3
for t in range(1):
    A = 2**4 * 5**t
    E = EllipticCurve([0, 0, 10*A, 0, -5*A**2])
    pts = E.integral_points(both_signs=True)
    for P in pts:
        if (P[0] % A == 0) and (P[1] % A == 0):
            x0 = P[1] / A
            if x0.is_square():
                x0 = sqrt(x0)
                y0 = ((x0-1)/2)^5 + ((x0 + 1)/2)^5
                if y0**(1/3) in ZZ or (-y0)**(1/3) in ZZ:
                    print(f"x0: {x0}, y0: {y0.factor()} for t={t}")

