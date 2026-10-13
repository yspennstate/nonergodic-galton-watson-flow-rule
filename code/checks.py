"""Checks for "A nonergodic flow rule with Galton-Watson stationary measure".

The Galton-Watson tree has p_2 = 3/5 and p_3 = 2/5. A vertex with two children takes the AND of their colours, a
vertex with three children the OR. q_n(t) is the probability that the root of t gets colour 1 when generation n is
coloured i.i.d. Bernoulli(1/2).

Part S checks the algebra symbolically (SymPy): the fixed point 1/2, the closed family of pair laws, the AND and OR
disagreement polynomials, the recursion d' = d (3 - d)^2 / 10 and the bound d' <= 9d/10 on [0, 1].
Part E checks in exact rational arithmetic: the law of q_n for n <= 3 and the second moment of q_4 against the
recursion; the martingale identity E[q_{n+1} | T|n] = q_n on every tree of depth n <= 2 and on random trees of depth 3;
the mean matrix of the pivot tree; and the stationarity of GW under the flow rule on the sigma-field of the first two
generations, computed from the definition of the flow (uniform choice among the children of the root's colour).
Part M is a numerical sanity check by population dynamics: the law of q_n for n <= 80.
Requires Python 3 with SymPy and NumPy."""
import itertools
import random
from fractions import Fraction as F

import numpy as np
import sympy as sp

P = {2: F(3, 5), 3: F(2, 5)}
HALF = F(1, 2)
FAILS = []
NCHECK = [0]


def check(name, ok, detail=""):
    NCHECK[0] += 1
    print(("PASS " if ok else "FAIL ") + name + ("  " + detail if detail else ""))
    if not ok:
        FAILS.append(name)


def gate(k, a):
    """Probability that a vertex with k children has colour 1 when child i has colour 1 w.p. a[i], independently."""
    if k == 2:
        return a[0] * a[1]
    out = 1
    for x in a:
        out *= 1 - x
    return 1 - out


def gate_bits(k, c):
    return int(all(c)) if k == 2 else int(any(c))


# ---------------------------------------------------------------- Part S: symbolic
def part_s():
    d = sp.symbols("d")
    half = sp.Rational(1, 2)
    p2, p3 = sp.Rational(3, 5), sp.Rational(2, 5)
    check("S1 fixed point: (3/5)(1/2)^2 + (2/5)(1 - (1/2)^3) = 1/2",
          sp.simplify(p2 * half ** 2 + p3 * (1 - (1 - half) ** 3) - half) == 0)
    pair = {(1, 1): (1 - d) / 2, (0, 0): (1 - d) / 2, (1, 0): d / 2, (0, 1): d / 2}
    after = {}
    for k in (2, 3):
        law = {}
        for prs in itertools.product(pair, repeat=k):
            w = 1
            for pr in prs:
                w *= pair[pr]
            r = gate_bits(k, [pr[0] for pr in prs])
            r2 = gate_bits(k, [pr[1] for pr in prs])
            law[(r, r2)] = law.get((r, r2), 0) + w
        after[k] = {key: sp.expand(v) for key, v in law.items()}
    dis2 = sp.expand(after[2][(1, 0)] + after[2][(0, 1)])
    dis3 = sp.expand(after[3][(1, 0)] + after[3][(0, 1)])
    check("S2 AND of two pairs: P(R != R') = d - d^2/2", sp.expand(dis2 - (d - d ** 2 / 2)) == 0, str(dis2))
    check("S3 OR of three pairs: P(R != R') = (3d - 3d^2 + d^3)/4",
          sp.expand(dis3 - (3 * d - 3 * d ** 2 + d ** 3) / 4) == 0, str(dis3))
    mix = {key: sp.expand(p2 * after[2][key] + p3 * after[3][key]) for key in pair}
    dn = sp.expand(mix[(1, 0)] + mix[(0, 1)])
    check("S4 the mixed step keeps the family: P(1,0) = P(0,1) = d'/2, P(1,1) = P(0,0) = (1 - d')/2",
          all(sp.expand(mix[key] - v.subs(d, dn)) == 0 for key, v in pair.items()))
    check("S5 d' = 9d/10 - 3d^2/5 + d^3/10 = d (3 - d)^2 / 10",
          sp.expand(dn - (sp.Rational(9, 10) * d - sp.Rational(3, 5) * d ** 2 + d ** 3 / 10)) == 0
          and sp.expand(dn - d * (3 - d) ** 2 / 10) == 0, str(sp.factor(dn)))
    check("S6 d' - 9d/10 = -d^2 (6 - d)/10, <= 0 on [0, 6]",
          sp.expand(dn - sp.Rational(9, 10) * d + d ** 2 * (6 - d) / 10) == 0)
    marg = sp.expand(mix[(1, 1)] + mix[(1, 0)])
    check("S7 the marginal stays 1/2", sp.simplify(marg - half) == 0)


# ---------------------------------------------------------------- Part E: exact
def law_step(law):
    """Law of q_{n+1} from the law of q_n (dict value -> probability), i.i.d. children."""
    items = list(law.items())
    out = {}
    for k in (2, 3):
        for combo in itertools.product(items, repeat=k):
            w = P[k]
            for _, pw in combo:
                w *= pw
            v = gate(k, [x for x, _ in combo])
            out[v] = out.get(v, 0) + w
    return out


def d_rec(n):
    d = HALF
    for _ in range(n):
        d = d * (3 - d) ** 2 / 10
    return d


def part_e_law():
    law = {HALF: F(1)}
    for n in range(0, 4):
        mean = sum(v * w for v, w in law.items())
        dis = sum(2 * v * (1 - v) * w for v, w in law.items())
        check("E1.%d law of q_%d: total mass 1, E q = 1/2, E 2q(1-q) = d_%d (%d atoms)" % (n, n, n, len(law)),
              sum(law.values()) == 1 and mean == HALF and dis == d_rec(n), "d_%d = %s" % (n, d_rec(n)))
        law = law_step(law) if n < 3 else law
    # second moment of q_4 by enumeration over ordered pairs and triples of q_3 values
    items = list(law_step(law_step(law_step({HALF: F(1)}))).items())
    m1, m2 = F(0), F(0)
    for k in (2, 3):
        for combo in itertools.product(items, repeat=k):
            w = P[k]
            for _, pw in combo:
                w *= pw
            v = gate(k, [x for x, _ in combo])
            m1 += w * v
            m2 += w * v * v
    check("E1.4 E q_4 = 1/2 and E 2 q_4 (1 - q_4) = d_4 (%d atoms of q_3 enumerated)" % len(items),
          m1 == HALF and 2 * (m1 - m2) == d_rec(4), "d_4 = %s = %.6f" % (d_rec(4), float(d_rec(4))))


def shapes(depth):
    """All ordered trees of the given depth with degrees in {2, 3}; a vertex at the bottom is None."""
    if depth == 0:
        return [None]
    sub = shapes(depth - 1)
    out = []
    for k in (2, 3):
        out.extend(itertools.product(sub, repeat=k))
    return out


def weight(t):
    if t is None:
        return F(1)
    w = P[len(t)]
    for c in t:
        w *= weight(c)
    return w


def q(t):
    if t is None:
        return HALF
    return gate(len(t), [q(c) for c in t])


def bottom(t):
    if t is None:
        return 1
    return sum(bottom(c) for c in t)


def extend(t, degs):
    """Give the bottom vertices of t the degrees degs (in order); their children become bottom vertices."""
    it = iter(degs)

    def go(u):
        if u is None:
            return tuple([None] * next(it))
        return tuple(go(c) for c in u)
    return go(t)


def cond_next(t):
    """E[q_{n+1}(T) | T|n = t], by enumerating the degrees of generation n."""
    tot = F(0)
    for degs in itertools.product((2, 3), repeat=bottom(t)):
        w = F(1)
        for k in degs:
            w *= P[k]
        tot += w * q(extend(t, degs))
    return tot


def random_shape(depth, rng):
    if depth == 0:
        return None
    k = 2 if rng.random() < 0.6 else 3
    return tuple(random_shape(depth - 1, rng) for _ in range(k))


def part_e_martingale():
    for n in range(0, 3):
        ts = shapes(n)
        ok = all(cond_next(t) == q(t) for t in ts)
        mass = sum(weight(t) for t in ts)
        check("E2.%d martingale on all %d trees of depth %d (their GW mass sums to %s)" % (n, len(ts), n, mass),
              ok and mass == 1)
    rng = random.Random(20261013)
    tested = 0
    ok = True
    while tested < 40:
        t = random_shape(3, rng)
        if bottom(t) > 12:
            continue
        ok = ok and cond_next(t) == q(t)
        tested += 1
    check("E2.3 martingale on 40 random trees of depth 3 with at most 12 bottom vertices", ok)


def part_e_pivot():
    """Mean matrix of the pivot tree of Johnson, Podder and Skerman for the fixed point Bernoulli(1/2)."""
    M = {(r, s): F(0) for r in (0, 1) for s in (0, 1)}
    for k in (2, 3):
        for c in itertools.product((0, 1), repeat=k):
            w = P[k] / 2 ** k
            r = gate_bits(k, c)
            for i in range(k):
                c2 = list(c)
                c2[i] = 1 - c2[i]
                if gate_bits(k, c2) != r:
                    M[(r, c[i])] += w / F(1, 2)
    check("E3 pivot tree mean matrix = diag(9/10, 9/10), a pivotal child has its parent's colour",
          M == {(0, 0): F(9, 10), (0, 1): F(0), (1, 0): F(0), (1, 1): F(9, 10)}, str(M))


def part_e_stationary():
    """Law of the depth-2 shape of the next state, started from GW, computed from the flow's definition.

    Under GW the colour C(T^i) of child i, given the first three generations of T, is 1 with probability
    q_2(T^i), independently over i (the colours of generation 3 are i.i.d. Bernoulli(1/2) given T|3)."""
    sh = shapes(2)
    qs = {s: q(s) for s in sh}
    ws = {s: weight(s) for s in sh}
    out = {s: F(0) for s in sh}
    colour_one = F(0)
    for k in (2, 3):
        for kids in itertools.product(sh, repeat=k):
            wk = P[k]
            for s in kids:
                wk *= ws[s]
            for c in itertools.product((0, 1), repeat=k):
                w = wk
                for s, ci in zip(kids, c):
                    w *= qs[s] if ci else 1 - qs[s]
                if w == 0:
                    continue
                r = gate_bits(k, c)
                colour_one += w * r
                same = [i for i in range(k) if c[i] == r]
                for i in same:
                    out[kids[i]] += w / len(same)
    check("E4 P(C = 1) = 1/2 computed from the first three generations", colour_one == HALF)
    check("E5 one step of the flow from GW leaves the law of the first two generations unchanged (%d shapes)"
          % len(sh), all(out[s] == ws[s] for s in sh) and sum(out.values()) == 1)


# ---------------------------------------------------------------- Part M: population dynamics
def part_m():
    rng = np.random.default_rng(20261013)
    N = 400000
    x = np.full(N, 0.5)
    d = 0.5
    worst = 0.0
    for n in range(1, 81):
        k3 = rng.random(N) < 0.4
        a, b, c = x[rng.integers(0, N, N)], x[rng.integers(0, N, N)], x[rng.integers(0, N, N)]
        x = np.where(k3, 1 - (1 - a) * (1 - b) * (1 - c), a * b)
        d = d * (3 - d) ** 2 / 10
        worst = max(worst, abs(float(np.mean(2 * x * (1 - x))) - d))
    near = float(np.mean(np.minimum(x, 1 - x) < 1e-6))
    check("M1 population dynamics, 400000 samples: |mean 2q(1-q) - d_n| <= 0.005 for n <= 80",
          worst <= 0.005, "max deviation %.5f" % worst)
    check("M2 at n = 80 the mean of q is 1/2 +- 0.005 and at least 99% of samples are within 1e-6 of {0, 1}",
          abs(float(np.mean(x)) - 0.5) <= 0.005 and near >= 0.99,
          "mean %.5f, fraction %.5f, d_80 = %.3e" % (float(np.mean(x)), near, d))


def main():
    part_s()
    part_e_law()
    part_e_martingale()
    part_e_pivot()
    part_e_stationary()
    part_m()
    print("%d checks, %d failed%s" % (NCHECK[0], len(FAILS), (": " + ", ".join(FAILS)) if FAILS else ""))


if __name__ == "__main__":
    main()
