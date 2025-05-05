import logging
from sage.all import RealField, sqrt, log, floor

class DiophantineSystem:

    def __init__(self, f, d1, d2, prec=100, alpha1=0.5, verbose=0):

        self.logger = logging.getLogger(self.__class__.__name__)
        self.set_verbose(verbose)

        self._RR = RealField(prec=prec)
        self.f = f
        self.d1 = d1
        self.d2 = d2
        self.k = f.degree()
        self._ais = self._compute_ai()
        if self._ais[0] * self.d2**self.k < 0 or self.d1 < 0:
            raise ValueError(f"The conditions a0*d2^k > 0 and d1 > 0 are not satisfied!")
        self.k0 = self._compute_k0()
        self._D = (self._ais[0] * self.d2**self.k) / self.d1
        if self._D == 1:
            raise ValueError("The quantity a0*d2^k/d1 is equal to 1.")
        self._theta1, self._theta0 = self._compute_thetas()
        self.m = self._D.ceil()
        self._alpha1 = self._compute_alpha1(alpha1)
        self.b2 = self._lower_bound_b2()
        self.c0 = self._compute_c0()
        self.c1 = self._compute_c1()
        self.c2 = self._compute_c2()
        self.alpha = self._compute_alpha()
        self.lam = self._compute_lam()

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

    def _compute_ai(self):
        ais = [self.f.coefficient(self.k - i) for i in range(self.k + 1)]
        self.logger.debug(f"ais: {ais}")
        return ais

    def _compute_k0(self):
        k0 = min([i for i in range(1, self.k + 1) if self._ais[i] != 0])
        self.logger.debug(f"k0: {k0}")
        return k0

    def _lower_bound_b2(self):
        low1 = self._theta1**(1/self.k0)
        s = sum([abs(self._ais[i] * self.d2**(self.k - i)) for i in range(self.k0 + 1, self.k + 1)])
        low2 = sqrt(s / abs(self._ais[self.k0] * self.d2**(self.k - self.k0)))
        low1 = self._RR(low1)
        low2 = self._RR(low2)
        self.logger.debug(f"low1: {low1}, low2: {low2}")
        return max([low1, low2]).ceil() + 1

    def _compute_thetas(self):
        if self._D > 1:
            theta1 = sqrt(self._D)
            theta0 = 1
        else:
            theta1 = 1
            theta0 = sqrt(self._D)
        self.logger.debug(f"theta1: {theta1}, theta0: {theta0}")
        return theta1, theta0

    def _compute_c0(self):
        c0 = sum([abs(self._ais[i] * self.d2**(self.k - i) / (self.d1 * self.b2**(i-self.k0)))
                 for i in range(self.k0, self.k + 1)])
        self.logger.debug(f"c0: {c0}")
        return c0

    def _compute_c1(self):
        s = sum([abs(self._ais[i] * self.d2**(self.k - i)) for i in range(self.k0 + 1, self.k + 1)]) / self.b2**2
        s0 = abs(self._ais[self.k0] * self.d2**(self.k - self.k0))
        c1 = s0 - s
        self.logger.debug(f"c1 computation: (s, s0, C1) = ({s}, {s0}, {c1})")
        return c1

    def _compute_c2(self):
        num = abs(self.d1) * self.c0
        denom = self.b2**self.k0 * (log(self.b2**self.k0) - log(self._theta1))
        c2 = 1 + num / denom
        self.logger.debug(f"c2: {c2}")
        return c2

    def _compute_alpha1(self, alpha1):
        a1 = abs(self.m - self._D)
        if a1.is_zero():
            a1 = alpha1
        self.logger.debug(f"a1: {a1}")
        return a1

    def _compute_alpha(self):
        alpha = self._alpha1 + self.c0 / self.b2**self.k0
        self.logger.debug(f"alpha: {alpha}")
        return alpha

    def bound_n(self):
        self._compute_suitable_c0_b2_alpha()
        num = log((self._ais[0] * self.d2**self.k) / (self.d1 * (self.m + self.alpha)))
        bound = floor(num / log(self.lam))
        return bound

    def _compute_suitable_c0_b2_alpha(self):
        while self.lam >= 1 or self.alpha >= 1:
            self.b2 += 1
            self.logger.info(f"The new value of b2: {self.b2}")
            self.c0 = self._compute_c0()
            self.c1 = self._compute_c1()
            self.c2 = self._compute_c2()
            self.alpha = self._compute_alpha()
            self.lam = self._compute_lam()
        self.logger.info(f"The values of parameters are:\n"
                         f"lam: {self._RR(self.lam)}\n"
                         f"c0: {self._RR(self.c0)}\n"
                         f"c1: {self._RR(self.c1)}\n"
                         f"c2: {self._RR(self.c2)}\n"
                         f"b2: {self.b2}\n"
                         f"alpha: {self._RR(self.alpha)}")

    def _compute_lam(self):
        return self._RR(self.c2 * self.c0 * self._ais[0] * self.d2**self.k) / (self.c1 * (self.m + self.alpha))