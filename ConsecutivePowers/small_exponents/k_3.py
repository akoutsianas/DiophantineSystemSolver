from sage.all import ZZ, QQ, polygen, EllipticCurve

x = polygen(QQ, 'x')

# Case n = 2
E = EllipticCurve([3, 0])
pts = E.integral_points()
for P in pts:
    if P[0] % 2 == 1:
        print(f"x0: {(P[0] - 1)/2}")


# Case n = 3
for t in range(3):
    A = 2**2 * 3**t
    E = EllipticCurve([0, -3*A**2])
    pts = E.integral_points(both_signs=True)
    for P in pts:
        if (P[0] % A == 0) and (P[1] % A == 0):
            x0 = P[1] / A
            y0 = ((x0-1)/2)^3 + ((x0 + 1)/2)^3
            if y0**(1/3) in ZZ or (-y0)**(1/3) in ZZ:
                print(f"x0: {x0}, y0: {y0.factor()} for t={t}")
