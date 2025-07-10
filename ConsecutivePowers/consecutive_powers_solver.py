import logging
from sage.all import polygen, QQ, ZZ, ceil, prime_range, legendre_symbol

from diophantine_system_solver import DiophantineSystem

class ConsecutivePowersSolver:

    def __init__(self, k, verbose=0):
        self.logger = logging.getLogger(self.__class__.__name__)
        self.set_verbose(verbose)
        self._verbose = verbose

        if k % 4 == 0:
            raise ValueError(f"The value of k is a multiple of 4, so there are only the trivial solutions. "
                                f"This is known due to the Ellenberg's work.")

        if k == 2:
            raise ValueError(f"The case k=2 is solved with a different approach.")
        self._x = polygen(QQ, 'x')
        self.k = k
        self.fk = self._compute_fk()
        self._lower_bound_n = self._compute_lower_bound_n()
        self.logger.warning(f"For k={self.k} we assume that n > {max(self._lower_bound_n, 3)}!")

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

    def _compute_fk(self):
        g = (self._x - 1)**self.k + (self._x + 1)**self.k
        g /= 2
        if self.k % 2 == 1:
            fk = (g / self._x).numerator()
        else:
            fk = (g / (self._x**2 + 1)).numerator()
        fk = sum([a * self._x ** i for i, a in enumerate(fk.coefficients())])
        self.logger.info(f"The polynomial f{self.k} is equal to {fk}.")
        return fk

    def _compute_lower_bound_n(self):
        return 3*max([self.k.valuation(q) for q in self.k.prime_factors()])

    def _compute_d1_d2(self):
        if self.k % 2 == 1:
            return self._compute_d1_d2_case_odd()
        else:
            return self._compute_d1_d2_case_even()

    def _compute_d1_d2_case_odd(self):
        k_primes_dict = {q: self.k.valuation(q) for q in self.k.prime_factors()}
        self.logger.debug(f"The information about the prime factors of k are {k_primes_dict}.")
        pairs = [(1, 1)]

        for q in k_primes_dict.keys():
            pairs_q = []

            # Case e1q <= e2q
            for e1q in range(k_primes_dict[q] + 1):
                e2q = -e1q
                pairs_q.append((q ** e1q, q ** e2q))

            # Case e2q <= e1q
            for e2q in range(k_primes_dict[q] + 1):
                if 2*e2q <= k_primes_dict[q]:
                    e1q = -e2q
                    pairs_q.append((q ** e1q, q ** e2q))
            self.logger.debug(f"Pairs_q: {pairs_q}")

            pairs_new = []
            for pair in pairs:
                for pair_q in pairs_q:
                    d1 = pair[0] * pair_q[0]
                    d2 = pair[1] * pair_q[1]
                    if (d1, d2) not in pairs_new:
                        pairs_new.append((d1, d2))
            pairs = pairs_new

        return pairs

    def _compute_d1_d2_case_even(self):
        k_primes_dict = {q: self.k.valuation(q) for q in ZZ(self.k/2).prime_factors()
                         if legendre_symbol(-1, q) == 1}
        self.logger.debug(f"The information about the prime factors of k are {k_primes_dict}.")
        pairs = [(1, 1)]

        for q in k_primes_dict.keys():
            pairs_q = []

            # Case e1q <= e2q
            for e1q in range(k_primes_dict[q] + 1):
                e2q = -e1q
                pairs_q.append((q ** e1q, q ** e2q))

            # Case e2q <= e1q
            for e2q in range(k_primes_dict[q] + 1):
                e1q = -e2q
                pairs_q.append((q ** e1q, q ** e2q))
            self.logger.debug(f"Pairs_q: {pairs_q}")

            pairs_new = []
            for pair in pairs:
                for pair_q in pairs_q:
                    d1 = pair[0] * pair_q[0]
                    d2 = pair[1] * pair_q[1]
                    if (d1, d2) not in pairs_new:
                        pairs_new.append((d1, d2))
            pairs = pairs_new

        return pairs

    def solve_equation(self, sieve_method='powers', p_iter='multiprocessing', ncpus=2, prec=100):
        pairs = self._compute_d1_d2()
        f = self._homogenous_equation()
        n0 = self._lower_bound_n
        sols = []
        for (d1, d2) in pairs:
            self.logger.info(f"We solve the system for d1={d1} and d2={d2}.")
            if self.k % 2 == 1:
                ds = DiophantineSystem(f, 2**(self.k - 1) * d1, d2, prec=prec, verbose=self._verbose)
            else:
                ds = DiophantineSystem(f, 2 ** (self.k - 2) * d1, 2*d2, prec=prec, verbose=self._verbose)
            x0 = ceil(max(ds.x1, ds.x2)) + 1
            self.logger.info(f"X0: {x0}.")
            n1 = ds.bound_n_for_x0(x0)
            sols_ds = ds.sieve(x0, method=sieve_method, p_iter=p_iter, ncpus=ncpus)
            self.logger.debug(f"Solutions of the Diophantine system after the sieve: {sols_ds}")
            if self.k % 2 == 0:
                sol_temp = []
                for sol in sols_ds:
                    if ZZ(sol - 1).is_square():
                        for rt, _ in (self._x**2 - (sol - 1)).roots():
                            if rt not in sol_temp:
                                sol_temp.append(rt)
                sols_ds = sol_temp
            self.logger.info(f"We have the small solutions {sols_ds}.")
            for v in sols_ds:
                if v not in sols:
                    sols.append(v)
            self.logger.debug(f"We have n1={n1}.")
            n0 = max(n1, n0)
        n0 = max(prime_range(3, n0+1))
        self.logger.info(f"We have n>={n0}.")
        return n0, sols

    def _homogenous_equation(self):
        if self.k % 2 == 1:
            f = self.fk(self._x**2)
        else:
            f = self.fk(self._x - 1)
        self.logger.debug(f"fk equal to {f}.")
        return f
