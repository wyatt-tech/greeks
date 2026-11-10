import math
import os
import sys
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import greeks as g  # noqa: E402


class TestGreeks(unittest.TestCase):
    def test_every_greek_matches_finite_differences(self):
        self.assertLess(g.fd_check(), 1e-3)

    def test_put_call_parity(self):
        for (S, K, T, r, q, sig) in g.CASES:
            c = g.price(S, K, T, r, q, sig, True)
            p = g.price(S, K, T, r, q, sig, False)
            self.assertAlmostEqual(c - p, S * math.exp(-q * T) - K * math.exp(-r * T), places=9)

    def test_call_and_put_share_gamma_and_vega(self):
        for (S, K, T, r, q, sig) in g.CASES:
            c, p = g.greeks(S, K, T, r, q, sig, True), g.greeks(S, K, T, r, q, sig, False)
            self.assertAlmostEqual(c["gamma"], p["gamma"], places=12)
            self.assertAlmostEqual(c["vega"], p["vega"], places=12)
            self.assertAlmostEqual(c["delta"] - p["delta"], math.exp(-q * T), places=12)

    def test_textbook_value(self):
        # Hull's example: S=42, K=40, T=0.5, r=10%, sig=20%, no dividends -> call 4.76, put 0.81
        self.assertAlmostEqual(g.price(42, 40, 0.5, 0.10, 0.0, 0.20, True), 4.76, places=2)
        self.assertAlmostEqual(g.price(42, 40, 0.5, 0.10, 0.0, 0.20, False), 0.81, places=2)

    def test_expiry_returns_intrinsic(self):
        self.assertEqual(g.price(105, 100, 0, 0.05, 0, 0.2, True), 5)
        self.assertEqual(g.price(95, 100, 0, 0.05, 0, 0.2, True), 0)

    def test_all_names_present(self):
        out = g.greeks(100, 100, 0.5, 0.03, 0.01, 0.25, True)
        for name in g.GREEK_NAMES:
            self.assertIn(name, out)


if __name__ == "__main__":
    unittest.main()
