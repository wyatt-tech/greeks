# greeks

Black-Scholes-Merton prices and 18 analytic greeks in one dependency-free Python file, with every formula checked against a
finite-difference derivative.

| order | greeks |
|---|---|
| first | delta, gamma, vega, theta, rho, epsilon (psi), lambda |
| second | vanna, charm, vomma, veta, vera |
| third | speed, zomma, color, ultima |
| dual (strike) | dual delta, dual gamma |

```python
from greeks import price, greeks

price(42, 40, 0.5, 0.10, 0.0, 0.20, call=True)          # 4.759 (Hull's textbook example)
g = greeks(100, 105, 30 / 365, 0.045, 0.013, 0.22, call=False)
g["delta"], g["gamma"], g["vanna"], g["charm"]
```

Arguments: spot `S`, strike `K`, time `T` in years, rate `r`, dividend yield `q`, volatility `sig` (all annualized,
continuous compounding), and `call` (True or False). Greeks are per unit: divide vega and rho by 100 for "per 1%", and
theta by 365 for "per day".

## Verification

`python greeks.py` compares every analytic greek with a central finite difference (of the price for first-order greeks,
of the analytic lower-order greek for the higher orders) on calls and puts across four cases, from a 7-day option to a
1-year one:

```
worst relative error across all greeks and cases: 2.73e-04  ALL PASS (< 1e-3)
```

Unit tests (put-call parity, the call/put identities, a textbook price, expiry behaviour):
`python -m unittest discover -s tests -v`

Built with AI coding tools under my direction; I designed the checks and verified the results. MIT license.
