"""Tests for pipeline/stability_map.py on a synthetic map (no real scores, no patient data)."""
import unittest

from pipeline import stability_map as sm


def row(pos, wt, mut, ddg, plddt=90.0, rsa=0.1, llr=-1.0):
    return {"pos": float(pos), "wt": wt, "mut": mut, "esm_llr": llr, "ddg_destab": ddg, "ddg_raw": ddg,
            "plddt": plddt, "rsa": rsa, "contacts": 10.0}


def synthetic(control_ddg=3.0, no_defect_ddg=0.2):
    m = {}
    for i in range(100):  # reference set: positions 800-899, ddG 0.0 .. 9.9
        m[f"A{800 + i}G"] = row(800 + i, "A", "G", i / 10, rsa=0.5)
    for k in sm.ABUNDANCE:
        m[k] = row(int(k[1:-1]), k[0], k[-1], control_ddg)
    for k in sm.NO_DEFECT:
        m[k] = row(int(k[1:-1]), k[0], k[-1], no_defect_ddg, plddt=50.0)
    m["V700A"] = row(700, "V", "A", 9.0)  # before residue 727: outside the reference set
    return m


class StabilityMapTest(unittest.TestCase):
    def test_supports_high_ddg(self):
        self.assertEqual(sm.call(synthetic(), "A890G")["call"], "supports")

    def test_argues_against_low_ddg_exposed(self):
        self.assertEqual(sm.call(synthetic(), "A805G")["call"], "argues against")

    def test_low_ddg_buried_is_inconclusive(self):
        m = synthetic()
        m["A805G"]["rsa"] = 0.1
        self.assertEqual(sm.call(m, "A805G")["call"], "inconclusive")

    def test_middle_is_inconclusive(self):
        self.assertEqual(sm.call(synthetic(), "A850G")["call"], "inconclusive")

    def test_gate_failure_blocks_the_call(self):
        res = sm.call(synthetic(control_ddg=0.1, no_defect_ddg=5.0), "A890G")
        self.assertFalse(res["gate"]["passes"])
        self.assertTrue(res["call"].startswith("not interpreted"))

    def test_reference_set_excludes_early_residues(self):
        res = sm.call(synthetic(), "A850G")
        self.assertEqual(res["reference_n"], 105)  # 100 reference rows + 5 abundance controls (>= 727, pLDDT 90)

    def test_unknown_substitution(self):
        with self.assertRaises(ValueError):
            sm.call(synthetic(), "W1W")


if __name__ == "__main__":
    unittest.main()
