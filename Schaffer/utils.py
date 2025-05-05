from sage.combinat.q_bernoulli import q_bernoulli

# Compute d1 and d2 for the Schaffer conjecture

def compute_Tk_Ck(k):

    x = polygen(QQ, 'x')
    f = bernoulli_polynomial(x, k+ 1)
    B = q_bernoulli(k + 1)(1)
    Sk = (f - B)/(k + 1)
    Sk = Sk(x + 1)
    Ck = Sk.denominator()
    fk = Sk.numerator()
    R = fk.parent()
    # print(f"fk: {fk.factor()}, Ck: {Ck}")
    if k % 2 == 0:
        Tk = fk/(x * (x + 1) * (2*x + 1))
    else:
        Tk = fk/(x**2 * (x + 1)**2)
    Tk = R(Tk)
    # print(f"fk: {fk.factor()}, Tk:{Tk}, Ck: {Ck}")
    return Tk, Ck

def compute_d1_d2(k):

    if k % 2 == 0:
        return compute_d1_d2_k_even(k)
    else:
        return compute_d1_d2_k_odd(k)

def compute_d1_d2_k_even(k):
    Tk, Ck = compute_Tk_Ck(k)
    print(f"Tk(0): {Tk(0).factor()}")
    print(f"Ck: {Ck.factor()}")
    x = Tk.parent().gen()
    res1 = Tk.resultant(x)
    res2 = Tk.resultant(x + 1)
    res3 = Tk.resultant(2*x + 1)
    print(f"res1: {res1.factor()}")
    print(f"res2: {res2.factor()}")
    print(f"res3: {res3.factor()}")
    Ck_primes = {p: Ck.valuation(p) for p in Ck.prime_factors()}
    res1_primes = {p: res1.valuation(p) for p in res1.prime_factors()}
    res2_primes = {p: res2.valuation(p) for p in res2.prime_factors()}
    res3_primes = {p: res3.valuation(p) for p in res3.prime_factors()}
    prs = list(set(list(Ck_primes.keys()) + list(res1_primes.keys()) + list(res2_primes.keys()) + list(res3_primes.keys())))

    Ds = [(1,1)]
    for p in prs:
        dp_exp = []
        if p in Ck_primes.keys():
            vpCk = Ck_primes[p]

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
        pairs = []
        for ds in Ds:
            for dp in dp_pairs:
                pairs.append((ds[0]*dp[0], ds[1]*dp[1]))
        Ds = pairs
        print(f"p: {p}, dp_exp: {dp_exp}, dp_pairs: {dp_pairs}")
    print(f'Ds: {Ds}')

def compute_d1_d2_k_odd(k):
    Tk, Ck = compute_Tk_Ck(k)
    x = Tk.parent().gen()

