#!/usr/bin/env python3
r"""greeks.py -- Black-Scholes-Merton prices and 18 analytic greeks (first, second and third order, plus the dual greeks),
with a finite-difference check of every formula.

    python greeks.py        # verifies every analytic greek against a numerical derivative of price()

Conventions: T in years, sig and r and q annualized (continuous compounding). Greeks are per unit: vega per 1.0 of vol,
theta, rho and epsilon per year and per 1.0 of rate. Divide vega and rho by 100 for "per 1%", theta by 365 for "per day".
"""
import math

SQRT2PI = math.sqrt(2.0 * math.pi)

GREEK_NAMES = ("delta", "gamma", "vega", "theta", "rho", "epsilon", "lambda", "vanna", "charm", "vomma", "veta", "vera",
               "speed", "zomma", "color", "ultima", "dual_delta", "dual_gamma")


def ncdf(x): return 0.5 * (1.0 + math.erf(x / math.sqrt(2.0)))
def npdf(x): return math.exp(-0.5 * x * x) / SQRT2PI


def d1d2(S, K, T, r, q, sig):
    v = sig * math.sqrt(T)
    d1 = (math.log(S / K) + (r - q + 0.5 * sig * sig) * T) / v
    return d1, d1 - v


def price(S, K, T, r, q, sig, call):
    """European option price. At or past expiry (T <= 0) or with zero vol it returns the intrinsic value."""
    if T <= 0 or sig <= 0:
        return max(0.0, (S - K) if call else (K - S))
    d1, d2 = d1d2(S, K, T, r, q, sig)
    if call:
        return S * math.exp(-q * T) * ncdf(d1) - K * math.exp(-r * T) * ncdf(d2)
    return K * math.exp(-r * T) * ncdf(-d2) - S * math.exp(-q * T) * ncdf(-d1)


def greeks(S, K, T, r, q, sig, call):
    """All analytic greeks for one option. Returns a dict (plus d1 and d2). Requires T > 0 and sig > 0."""
    sqT = math.sqrt(T)
    v = sig * sqT
    d1, d2 = d1d2(S, K, T, r, q, sig)
    nd1 = npdf(d1)
    eqt, ert = math.exp(-q * T), math.exp(-r * T)
    sign = 1.0 if call else -1.0
    Nsd1 = ncdf(sign * d1)        # N(d1) for a call, N(-d1) for a put
    Nsd2 = ncdf(sign * d2)

    g = {"d1": d1, "d2": d2}
    # first order
    g["delta"] = sign * eqt * Nsd1
    g["gamma"] = eqt * nd1 / (S * v)
    g["vega"] = S * eqt * nd1 * sqT
    g["theta"] = (-eqt * S * nd1 * sig / (2.0 * sqT)
                  - sign * r * K * ert * Nsd2 + sign * q * S * eqt * Nsd1)
    g["rho"] = sign * K * T * ert * Nsd2
    g["epsilon"] = -sign * S * T * eqt * Nsd1            # psi = dPrice/dq
    p = price(S, K, T, r, q, sig, call)
    g["lambda"] = g["delta"] * S / p if p > 0 else float("nan")

    # second order
    g["vanna"] = -eqt * nd1 * d2 / sig                   # dDelta/dSig = dVega/dS
    g["charm"] = (sign * q * eqt * Nsd1                  # -dDelta/dT (delta decay)
                  - eqt * nd1 * (2.0 * (r - q) * T - d2 * v) / (2.0 * T * v))
    g["vomma"] = g["vega"] * d1 * d2 / sig               # dVega/dSig
    g["veta"] = (-S * eqt * nd1 * sqT
                 * (q + (r - q) * d1 / v - (1.0 + d1 * d2) / (2.0 * T)))   # dVega/dT
    g["vera"] = -K * T * ert * npdf(d2) * d1 / sig       # dRho/dSig

    # third order
    g["speed"] = -g["gamma"] / S * (d1 / v + 1.0)        # dGamma/dS
    g["zomma"] = g["gamma"] * (d1 * d2 - 1.0) / sig      # dGamma/dSig
    g["color"] = (-eqt * nd1 / (2.0 * S * T * v)         # dGamma/dT
                  * (2.0 * q * T + 1.0 + (2.0 * (r - q) * T - d2 * v) * d1 / v))
    g["ultima"] = (-g["vega"] / (sig * sig)              # dVomma/dSig
                   * (d1 * d2 * (1.0 - d1 * d2) + d1 * d1 + d2 * d2))

    # dual greeks (with respect to the strike)
    g["dual_delta"] = -sign * ert * Nsd2                 # dPrice/dK
    g["dual_gamma"] = ert * npdf(d2) / (K * v)           # d2Price/dK2
    return g


CASES = [(100, 95, 0.25, 0.045, 0.013, 0.20), (100, 110, 0.05, 0.03, 0.0, 0.45),
         (450, 455, 0.02, 0.05, 0.006, 0.15), (50, 48, 1.0, 0.04, 0.02, 0.30)]


def numeric_greeks(S, K, T, r, q, sig, call, h=1e-4):
    """Central finite differences of price(), and of the analytic lower-order greeks for the higher orders."""
    P = lambda **b: price(b.get("S", S), b.get("K", K), b.get("T", T), b.get("r", r), b.get("q", q), b.get("sig", sig), call)
    G = lambda name, **b: greeks(b.get("S", S), b.get("K", K), b.get("T", T), r, q, b.get("sig", sig), call)[name]
    n = {}
    n["delta"] = (P(S=S + h) - P(S=S - h)) / (2 * h)
    n["gamma"] = (P(S=S + h) - 2 * P() + P(S=S - h)) / (h * h)
    n["vega"] = (P(sig=sig + h) - P(sig=sig - h)) / (2 * h)
    n["theta"] = -(P(T=T + h) - P(T=T - h)) / (2 * h)
    n["rho"] = (P(r=r + h) - P(r=r - h)) / (2 * h)
    n["epsilon"] = (P(q=q + h) - P(q=q - h)) / (2 * h)
    n["dual_delta"] = (P(K=K + h) - P(K=K - h)) / (2 * h)
    n["dual_gamma"] = (P(K=K + h) - 2 * P() + P(K=K - h)) / (h * h)
    n["vanna"] = (G("delta", sig=sig + h) - G("delta", sig=sig - h)) / (2 * h)
    n["charm"] = -(G("delta", T=T + h) - G("delta", T=T - h)) / (2 * h)
    n["vomma"] = (G("vega", sig=sig + h) - G("vega", sig=sig - h)) / (2 * h)
    n["veta"] = (G("vega", T=T + h) - G("vega", T=T - h)) / (2 * h)
    n["vera"] = (G("rho", sig=sig + h) - G("rho", sig=sig - h)) / (2 * h)
    n["speed"] = (G("gamma", S=S + h) - G("gamma", S=S - h)) / (2 * h)
    n["zomma"] = (G("gamma", sig=sig + h) - G("gamma", sig=sig - h)) / (2 * h)
    n["color"] = (G("gamma", T=T + h) - G("gamma", T=T - h)) / (2 * h)
    n["ultima"] = (G("vomma", sig=sig + h) - G("vomma", sig=sig - h)) / (2 * h)
    return n


def fd_check(cases=CASES, verbose=False):
    """Worst relative error between every analytic greek and its numerical counterpart."""
    worst = 0.0
    for (S, K, T, r, q, sig) in cases:
        for call in (True, False):
            a = greeks(S, K, T, r, q, sig, call)
            for name, nv in numeric_greeks(S, K, T, r, q, sig, call).items():
                rel = abs(a[name] - nv) / max(1e-8, abs(nv))
                worst = max(worst, rel)
                if verbose and rel > 1e-3:
                    print(f"  MISMATCH {name:10} {a[name]:>14.6f} {nv:>14.6f} {rel:>10.2e}  ({'C' if call else 'P'} S={S} K={K} T={T})")
    return worst


if __name__ == "__main__":
    w = fd_check(verbose=True)
    print(f"worst relative error across all greeks and cases: {w:.2e}  " + ("ALL PASS (< 1e-3)" if w < 1e-3 else "*** FAILURES ABOVE ***"))
