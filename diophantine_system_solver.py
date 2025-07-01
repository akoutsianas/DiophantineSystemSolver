import logging
from sage.all import RealField, sqrt, log, floor, ZZ, exp

class DiophantineSystem:

    def __init__(self, f, d1, d2, x0=1, prec=100, alpha1=0.5, eps=0.5, verbose=0):
        """
        We want to solve the Diophantine system f(x)=d_1y_1^n, x=d_2y_2^n. We give an upper bound n0 of n and
        the solutions when n > n0.

        :param f: a polynomial with integer coefficients
        :param d1: an integer
        :param d2: an integer
        :param x0: an integer. We assume that |x|>=x0
        :param prec: the precision for the real number computations
        :param alpha1: a value of a' if a0d2^k/d1 is an integer
        :param eps: an upper bound of the difference between a and a'
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
        if self._ais[0] * self.d2**self.k < 0 or self.d1 < 0:
            raise ValueError(f"The conditions a0*d2^k > 0 and d1 > 0 are not satisfied!")
        self.k0 = self._compute_k0()
        self._D = (self._ais[0] * self.d2**self.k) / self.d1
        if self._D == 1:
            raise ValueError("The quantity a0*d2^k/d1 is equal to 1.")
        self._theta1, self._theta0 = self._compute_thetas()
        self.m = self._compute_m()
        self._alpha1 = self._compute_alpha1()

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
        r = self._ais[0] * self.d2**self.k / self.d1
        if r in ZZ:
            m = r + 0.5
        else:
            m = r.ceil()
        self.logger.debug(f"m: {m}")
        return m

    def _compute_thetas(self):
        if self._D > 1:
            theta1 = sqrt(self._D)
            theta0 = 1
        else:
            theta1 = 1
            theta0 = sqrt(self._D)
        self.logger.debug(f"theta1: {theta1}, theta0: {theta0}")
        return theta1, theta0

    def _compute_c0_c1_c2(self, x0):
        c0 = sum([abs(self._cis[self.k0 + i])/x0**i for i in range(self.k - self.k0 + 1)])
        self.logger.debug(f"c0: {c0}")

        if self.k == self.k0:
            c1 = abs(self.d2**self.k * self._ais[self.k])
        else:
            c1 = abs(self._ais[self.k0]) - sum([abs(self._ais[self.k0+i+1])/x0**(i + 1) for i in range(self.k - self.k0)])
            c1 *= abs(self.d2**self.k)
        self.logger.debug(f"c1: {c1}")

        c2 = exp(c1 / (self._theta0 * x0**(self.k0)))
        self.logger.debug(f"c2: {c2}")
        return c0, c1, c2

    def _compute_alpha1(self):
        a1 = self.m - self._ais[0] * self.d2**self.k / self.d1
        self.logger.debug(f"a1: {a1}")
        return a1

    def _compute_alpha(self, x0, c0):
        alpha = self._alpha1 + c0 / x0**self.k0
        self.logger.debug(f"alpha: {alpha}")
        return alpha

    def _compute_x1(self, c1):
        x1 = c1 * max(1/(1 - self._alpha1), 1/(self._theta1 - self._theta0))
        x1 = x1**(1/self.k0)
        self.logger.debug(f"x1: {x1}")
        return x1

    def _compute_x2(self, ):
        x2 = sum([abs(self._ais[i]) for i in range(self.k0+1, self.k+1)]) / abs(self._ais[self.k0])
        self.logger.debug(f"x2: {x2}")
        return x2

    def _compute_lam(self, c0, c1, c2, alpha):
        lam = self._RR(c2 * c0 * self._ais[0] * self.d2**self.k) / (c1 * (self.m + alpha))
        self.logger.info(f"lam: {lam}")
        return lam

    def bound_n_for_x0(self, x0):
        c0, c1, c2 = self._compute_c0_c1_c2(x0)
        alpha = self._compute_alpha(x0, c0)
        x1 = self._compute_x1(c1)
        x2 = self._compute_x2()
        if (x0 > x1) and (x0 > x2):
            lam = self._compute_lam(c0, c1, c2, alpha)
            num = log((self._ais[0] * self.d2**self.k) / (self.d1 * (self.m + alpha)))
            bound = floor(num / log(lam))
            self.logger.info(f"The bound is {bound}")
        else:
            self.logger.warning(f"It holds x0 <= x1 or x0 <= x2")
            bound = None
        return bound


