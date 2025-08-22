from sage.all import ZZ, QQ, polygen, EllipticCurve

x = polygen(QQ, 'x')

# Case n = 2: We can reduce to the computation on integral points on the genus 1 curves
# x^4 + 14x^2 + 1 = 2^4 * 3^t * y_1^2 for t=0, 1


# Case n = 3
for t in range(3):
    A = 2 * 3**t
    E = EllipticCurve([0, -A**2])
    for pt in E.integral_points(both_signs=True):
        if pt[1] % A == 0:
            x0 = pt[1] / A
            y0n = (((x0-1)/2)**6 + ((x0 + 1)/2)**6)
            if (y0n in ZZ) and (y0n.is_perfect_power()):
                print(f"Value of x is {x0}.")

