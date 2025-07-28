import logging

from sage.all import RealField, QQ, log, floor, ZZ, exp, prime_range
from sage.parallel.decorate import parallel


def internal_iteration(n, lower_bound, upper_bound, d1, d2, f, logger):
    import time
    solutions = []
    t0 = time.time()
    for y2 in range(lower_bound, upper_bound + 1):
        x = d2 * y2 ** n
        y1 = (QQ(f(x)) / d1)
        if y1 in ZZ:
            if y1.is_perfect_power():
                y1 = y1 ** (1 / n)
                if y1 in ZZ and x not in solutions:
                    solutions.append(x)
    t1 = time.time()
    logger.info(f"Time for n={n}, lower_bound={lower_bound} and upper_bound={upper_bound}: {t1 - t0}")
    return solutions


class DiophantineSystem:

    def __init__(self, f, d1, d2, prec=100, verbose=0):
        """
        We want to solve the Diophantine system f(x)=d_1y_1^n, x=d_2y_2^n. We give an upper bound n0 of n and
        the solutions when n > n0.

        :param f: a polynomial with integer coefficients
        :param d1: an integer
        :param d2: an integer
        :param prec: the precision for the real number computations
        :param verbose: 0, 1, 2
        """
        self.logger = logging.getLogger(self.__class__.__name__)
        self.set_verbose(verbose)

        self._RR = RealField(prec)
        self.f = f
        self.d1 = d1
        self.d2 = d2
        self.k = f.degree()
        self._ais = self._compute_ais()
        self._cis = self._compute_cis()
        if (self._ais[0] * self.d2**self.k < 0) or (self.d1 < 0):
            raise ValueError(f"The conditions a0*d2^k > 0 and d1 > 0 are not satisfied!")
        self.k0 = self._compute_k0()
        self._thetan = (self._ais[0] * self.d2**self.k) / self.d1
        if self._thetan == 1:
            raise ValueError("The quantity a0*d2^k/d1 is equal to 1.")
        self._theta1, self._theta0 = self._compute_thetas()
        self.m = self._compute_m()
        self._alpha1 = self._compute_alpha1()
        self.x1 = self._compute_x1()
        self.x2 = self._compute_x2()

    def set_verbose(self, verbose):
        # Map integers to logging levels
        level_map = {
            0: logging.WARNING,
            1: logging.INFO,
            2: logging.DEBUG
        }
        if verbose not in level_map:
            raise ValueError(f"Invalid verbose value '{verbose}'. Allowed values are 0 (WARNING), 1 (INFO), 2 (DEBUG).")

        logging.basicConfig(
            format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
            level=level_map[verbose]
        )
        self.logger.setLevel(level_map[verbose])

    def _compute_ais(self):
        ais = [self.f.coefficient(self.k - i) for i in range(self.k + 1)]
        self.logger.debug(f"ais: {ais}")
        return ais

    def _compute_k0(self):
        k0 = min([i for i in range(1, self.k + 1) if self._ais[i] != 0])
        self.logger.debug(f"k0: {k0}")
        return k0

    def _compute_cis(self):
        cis = [self._ais[i] * self.d2**self.k / self.d1 for i in range(self.k + 1)]
        self.logger.debug(f"cis: {cis}")
        return cis

    def _compute_m(self):
        r = self._thetan
        if r in ZZ:
            m = r + 0.5
        else:
            m = r.ceil()
        self.logger.debug(f"m: {m}")
        return m

    def _compute_thetas(self):
        if self._thetan > 1:
            theta1 = (self._RR(self._thetan)) ** (1/QQ(3))
            theta0 = 1
        else:
            theta1 = 1
            theta0 = (self._RR(self._thetan)) ** (1/QQ(3))
        self.logger.debug(f"theta1: {theta1}, theta0: {theta0}")
        return theta1, theta0

    def _compute_c0_c1_c2(self, x0):
        c0 = sum([abs(self._cis[i])/x0**(i - self.k0) for i in range(self.k0, self.k + 1)])
        if self.k == self.k0:
            c1 = abs(self.d2**self.k * self._ais[self.k])
        else:
            c1 = abs(self._ais[self.k0]) - sum([abs(self._ais[i])/x0**(i-self.k0) for i in range(self.k0+1, self.k+1)])
            c1 *= abs(self.d2**self.k)
        c2 = self._RR(exp(c0 / (self._theta0 * x0**self.k0)))
        return c0, c1, c2

    def _compute_alpha1(self):
        a1 = self.m - self._thetan
        self.logger.debug(f"a1: {a1}")
        return a1

    def _compute_alpha(self, x0, c0):
        alpha = self._alpha1 + c0 / x0**self.k0
        return alpha

    def _compute_x1(self):
        x1 = max(1/(1 - self._alpha1), 1/(self._theta1 - self._theta0))
        x1 *= sum([abs(self._cis[i]) for i in range(self.k0, self.k+1)])
        x1 = self._RR(x1)**(1/self.k0)
        self.logger.debug(f"x1: {x1}")
        return x1

    def _compute_x2(self, ):
        x2 = sum([abs(self._ais[i]) for i in range(self.k0+1, self.k+1)]) / abs(self._ais[self.k0])
        self.logger.debug(f"x2: {x2}")
        return x2

    def _compute_lam(self, c0, c1, c2, alpha):
        lam = self._RR(c2 * c0 * self._ais[0] * self.d2**self.k) / (c1 * (self.m + alpha))
        return lam

    def bound_n_for_x0(self, x0, x0_step=100):
        c0, c1, c2 = self._compute_c0_c1_c2(x0)
        alpha = self._compute_alpha(x0, c0)
        lam = self._compute_lam(c0, c1, c2, alpha)
        while lam >= 1:
            x0 += x0_step
            c0, c1, c2 = self._compute_c0_c1_c2(x0)
            alpha = self._compute_alpha(x0, c0)
            lam = self._compute_lam(c0, c1, c2, alpha)
        self.logger.debug(f"c0: {c0}")
        self.logger.debug(f"c1: {c1}")
        self.logger.debug(f"c2: {c2}")
        self.logger.debug(f"alpha: {alpha}")
        self.logger.info(f"lam: {lam}")
        if (x0 > self.x1) and (x0 > self.x2) and (lam < 1):
            num = log(self._RR((self._ais[0] * self.d2**self.k) / (self.d1 * (self.m + alpha))))
            self.logger.debug(f"num: {num}")
            bound = floor(num / log(lam))
            self.logger.debug(f"Bound from computations: {bound}")
            bound = max(bound, 3)
            self.logger.info(f"The bound is {bound}")
        else:
            self.logger.warning(f"It holds x0 <= x1 or x0 <= x2 or lam >= 1")
            bound = None
        return bound, x0

    def sieve(self, x0, method='powers', p_iter='multiprocessing', ncpus=2):
        if method == 'naive':
            return self._sieve_naive(x0)
        elif method == 'powers':
            return self._sieve_powers(x0, p_iter=p_iter, ncpus=ncpus)
        return None

    def _sieve_naive(self, x0):
        sols = []
        for x in range(-x0, x0 + 1):
            w1 = self.f(ZZ(x)) / self.d1
            w2 = x / self.d2
            if w1.is_perfect_power() and w2.is_perfect_power():
                sols.append(x)
        return sols

    def _sieve_powers(self, x0, p_iter='multiprocessing', ncpus=2):
        sols = []
        n0 = max(floor(log(x0 / self.d2) / log(2)), 3)
        for n in prime_range(3, n0 + 1):
            y2_bound = floor((x0/self.d2)**(1/n))
            self.logger.debug(f"y2_bound: {y2_bound}")
            iter_length = floor((2*y2_bound + 1) / ncpus)
            limits = [(n, -y2_bound + i * iter_length, -y2_bound + (i + 1) * iter_length,
                       self.d1, self.d2, self.f, self.logger) for i in range(ncpus-1)]
            limits.append((n, -y2_bound + (ncpus - 1) * iter_length, y2_bound,
                           self.d1, self.d2, self.f, self.logger))
            self.logger.debug(f"limits: {limits}")
            sols_n = list(parallel(p_iter=p_iter, ncpus=ncpus)(internal_iteration)(limits))
            for s in sols_n:
                for sol in s[-1]:
                    if sol not in sols:
                        sols.append(sol)
        return sols

