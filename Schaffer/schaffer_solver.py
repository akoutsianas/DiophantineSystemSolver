import logging
from sage.all import polygen, QQ

from DiophantineSystemSolver.diophantine_system_solver import DiophantineSystem

from sage.combinat.q_bernoulli import q_bernoulli
from sage.all import bernoulli_polynomial

# Compute d1 and d2 for the Schaffer conjecture

class SchafferSolver:

    def __init__(self, k, verbose=0):
        self.logger = logging.getLogger(self.__class__.__name__)
        self.set_verbose(verbose)
        self.k = k
        self._compute_Tk_Ck()

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

        self._x = polygen(QQ, 'x')

    def _compute_Tk_Ck(self):

        f = bernoulli_polynomial(self._x, self.k + 1)
        B = q_bernoulli(self.k + 1)(1)
        Sk = (f - B)/(self.k + 1)
        Sk = Sk(self._x + 1)
        self.Ck = Sk.denominator()
        self.fk = Sk.numerator()
        R = fk.parent()

        if self.k % 2 == 0:
            self.Tk = self.fk/(self._x * (self._x + 1) * (2*self._x + 1))
        else:
            self.Tk = self.fk/(self._x**2 * (self._x + 1)**2)
        self.Tk = R(self.Tk)
        self.logger.info(f"We have fk: {self.fk.factor()}, Tk:{self.Tk}, Ck: {self.Ck}")

    def compute_d1_d2(self):

        if self.k % 2 == 0:
            return self._compute_d1_d2_k_even()
        else:
            return self._compute_d1_d2_k_odd()

    def _compute_d1_d2_k_even(self):
        self.logger.info(f"Factorization of Tk(0): {self.Tk(0).factor()}.")
        self.logger.info(f"Factorization of Ck: {self.Ck.factor()}.")

        res1 = self.Tk.resultant(self._x)
        res2 = self.Tk.resultant(self._x + 1)
        res3 = self.Tk.resultant(2*self._x + 1)
        self.logger.debug(f"Factorization of res1: {res1.factor()}")
        self.logger.debug(f"Factorization of res2: {res2.factor()}")
        self.logger.debug(f"Factorization of res3: {res3.factor()}")
        ck_primes = {p: self.Ck.valuation(p) for p in self.Ck.prime_factors()}
        res1_primes = {p: res1.valuation(p) for p in res1.prime_factors()}
        res2_primes = {p: res2.valuation(p) for p in res2.prime_factors()}
        res3_primes = {p: res3.valuation(p) for p in res3.prime_factors()}
        prs = list(set(list(ck_primes.keys()) + list(res1_primes.keys()) + list(res2_primes.keys()) + list(res3_primes.keys())))

        pairs = [(1,1)]
        for p in prs:
            dp_exp = []
            if p in ck_primes.keys():
                vpCk = ck_primes[p]

                if p in res1_primes.keys():
                    # The case p | (x, Tk)
                    for v1 in range(1, res1_primes[p] + 1):
                        dp_exp.append((v1, vpCk - v1))
                        dp_exp.append((vpCk - v1, v1))
                elif p in res2_primes.keys():
                    # The case p | (x + 1, Tk)
                    for v1 in range(1, res2_primes[p] + 1):
                        dp_exp.append((vpCk - v1, 0))
                        dp_exp.append((v1, 0))
                elif p in res3_primes.keys():
                    # The case p | (2*x + 1, Tk)
                    for v1 in range(1, res3_primes[p] + 1):
                        dp_exp.append((vpCk - v1, 0))
                        dp_exp.append((v1, 0))
                else:
                    # The case p \nmid (x*(x + 1)*(2*x + 1), Tk)
                    dp_exp.append((vpCk, 0))
                    dp_exp.append((0, vpCk))
            else:
                if p in res1_primes.keys():
                    # The case p | (x, Tk)
                    for v1 in range(1, res1_primes[p] + 1):
                        dp_exp.append((v1, -v1))
                        dp_exp.append((-v1, v1))
                elif p in res2_primes.keys():
                    # The case p | (x + 1, Tk)
                    for v1 in range(1, res2_primes[p] + 1):
                        dp_exp.append((-v1, 0))
                        dp_exp.append((v1, 0))
                elif p in res3_primes.keys():
                    # The case p | (2*x + 1, Tk)
                    for v1 in range(1, res3_primes[p] + 1):
                        dp_exp.append((-v1, 0))
                        dp_exp.append((v1, 0))
                else:
                    # The case p \nmid (x*(x + 1)*(2*x + 1), Tk)
                    dp_exp.append((0, 0))
            dp_exp = list(set(dp_exp))
            dp_pairs = [(p**v[0], p**v[1]) for v in dp_exp]
            pairs_temp = []
            for pair in pairs:
                for dp in dp_pairs:
                    pairs_temp.append((pair[0]*dp[0], pair[1]*dp[1]))
            pairs = pairs_temp
            self.logger.info(f"p: {p}, dp_exp: {dp_exp}, dp_pairs: {dp_pairs}")
        self.logger.info(f"All possible values pairs (d1, d2) are {pairs}.")

    def _compute_d1_d2_k_odd(self):
        pass
