#!/usr/bin/env python3
"""Offline tests for the R06 drug-triage flags and ranking. Records are synthetic paraphrases of label
language; nothing here is patient data and no network is used.

Run: python -m unittest discover -s pipeline/tests -t .
"""
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from pipeline.drug_triage import flags as F
from pipeline.drug_triage import sources as S
from pipeline.drug_triage import triage as T


def drug(moa="", action="INHIBITOR", targets=(), name="X"):
    mech = [{"moa": moa, "action": action, "targets": list(targets)}] if moa else []
    return {"chembl": "CHEMBL0", "name": name, "stage": "APPROVAL", "mechanisms": mech, "warnings": []}


def label(**sections):
    sections.setdefault("_generic", "DRUGX")
    sections.setdefault("adverse_reactions", "Headache.")
    return sections


class FlagTest(unittest.TestCase):
    def test_checkpoint_kinase_inhibitor_is_excluded(self):
        fl = F.evaluate(drug("Serine/threonine-protein kinase PLK1 inhibitor", targets=["PLK1"]), None, "BUB1B")
        self.assertEqual(fl["checkpoint_weakening"]["status"], "flag")
        self.assertEqual(F.tier(fl), "excluded")

    def test_lowering_the_dosage_deficient_gene_is_flagged(self):
        d = drug("Hypothetical degrader of GENEQ", action="DEGRADER", targets=["GENEQ"])
        self.assertEqual(F.checkpoint(d, "GENEQ")["status"], "flag")
        self.assertEqual(F.checkpoint(d, "GENEQ", dosage_deficient=False)["status"], "no flag found")

    def test_aneugen_classes(self):
        self.assertEqual(F.aneugen(drug("Tubulin beta chain inhibitor", targets=["TUBB"]), None)["status"], "flag")
        self.assertEqual(F.aneugen(drug("Kinesin KIF18A inhibitor", targets=["KIF18A"]), None)["status"], "flag")
        pan = drug("Histone deacetylase inhibitor", targets=["HDAC1", "HDAC2", "HDAC3", "HDAC6"])
        self.assertEqual(F.aneugen(pan, None)["status"], "flag")
        self.assertEqual(F.aneugen(drug("HDAC6 inhibitor", targets=["HDAC6"]), None)["status"], "no flag found")

    def test_label_genotoxicity_respects_negation(self):
        neg = label(nonclinical_toxicology="DrugX was not clastogenic in the in vitro chromosomal aberration assay.")
        pos = label(nonclinical_toxicology="DrugX induced micronuclei and aneuploidy in human lymphocytes.")
        self.assertEqual(F.aneugen(None, neg)["status"], "no flag found")
        self.assertEqual(F.aneugen(None, pos)["status"], "flag")
        # An abbreviation must not split the negation away from the assay name.
        abbrev = label(nonclinical_toxicology="There was no evidence of mutagenicity in the Ames test (S. typhimurium) or the chromosomal aberrations test (human lymphocytes).")
        self.assertEqual(F.aneugen(None, abbrev)["status"], "no flag found")

    def test_paediatric(self):
        yes = label(pediatric_use="8.4 Pediatric Use The safety and effectiveness of DRUGX in children aged 6 years and older have been established.")
        no = label(pediatric_use="Safety and effectiveness in pediatric patients have not been established.")
        self.assertEqual(F.paediatric(yes)["status"], "yes")
        self.assertEqual(F.paediatric(no)["status"], "not established")
        silent = label(pediatric_use="Children under two years are at increased risk of hepatotoxicity.")
        self.assertEqual(F.paediatric(silent)["status"], "not stated")
        narrow = label(pediatric_use="Efficacy in pediatric patients was not established for the treatment of migraine.")
        self.assertEqual(F.paediatric(narrow)["status"], "not stated")
        bracket = label(pediatric_use="The safety and effectiveness of DRUGX in children aged 6 years and older have been established [see Clinical Studies (14) ] . Safety and effectiveness below the age of 6 years have not been established.")
        self.assertEqual(F.paediatric(bracket)["status"], "yes")
        # Any label of the ingredient establishing paediatric use counts.
        multi = label(_pediatric_all=[no["pediatric_use"], yes["pediatric_use"]])
        self.assertEqual(F.paediatric(multi)["status"], "yes")
        self.assertEqual(F.paediatric(None)["status"], "not assessed")

    def test_cyp3a4_perpetrator_vs_victim(self):
        perp = label(drug_interactions="DRUGX is a strong inhibitor of CYP3A4 and increases exposure to CYP3A substrates.")
        victim = label(clinical_pharmacology="DrugX is metabolized primarily by CYP3A4.")
        denies = label(drug_interactions="DRUGX did not inhibit CYP3A4 in vitro.")
        self.assertEqual(F.cyp3a4(perp)["status"], "flag")
        self.assertEqual(F.cyp3a4(victim)["status"], "note")
        self.assertEqual(F.cyp3a4(denies)["status"], "no flag found")
        other_subject = label(_generic="SIROLIMUS", drug_interactions="Strong inhibitors of CYP3A4 increase sirolimus concentrations.")
        self.assertEqual(F.cyp3a4(other_subject)["status"], "no flag found")

    def test_kidney_and_marrow(self):
        renal = label(contraindications="Severe renal impairment (eGFR below 30 mL/min/1.73 m2).")
        self.assertEqual(F.kidney(renal)["status"], "flag")
        plt = label(warnings_and_cautions="Hematological changes: dose-related decreases in platelet counts (thrombocytopenia) occurred.")
        self.assertEqual(F.marrow(plt)["status"], "flag")
        ae_only = label(adverse_reactions="Neutropenia was reported in 1% of patients.")
        self.assertEqual(F.marrow(ae_only)["status"], "note")
        self.assertEqual(F.marrow(label(adverse_reactions="Headache."))["status"], "no flag found")
        self.assertEqual(F.marrow(label(description="A listing with no safety text.", adverse_reactions=""))["status"], "not assessed")
        ttp = label(warnings_and_cautions="Increased risk of thrombotic thrombocytopenic purpura with ciclosporin.")
        self.assertEqual(F.marrow(ttp)["status"], "no flag found")

    def test_tiers(self):
        clean = F.evaluate(drug("AMP-activated kinase activator", action="ACTIVATOR", targets=["PRKAA1"]),
                           label(pediatric_use="Safety and effectiveness have been established in pediatric patients."), "BUB1B")
        self.assertEqual(F.tier(clean), "screen-clear")
        self.assertEqual(F.tier(F.evaluate(None, None, "BUB1B")), "insufficient data")
        caution = F.evaluate(drug("A activator", action="ACTIVATOR", targets=["A"]),
                             label(warnings_and_cautions="Thrombocytopenia has occurred."), "BUB1B")
        self.assertEqual(F.tier(caution), "caution")


class RunTest(unittest.TestCase):
    def test_run_ranks_and_writes(self):
        cands = {"GENEQ": [], "MOD": [{"chembl": "C1", "name": "PANHDACI", "stage": "APPROVAL"}],
                 "PLK1": [{"chembl": "C2", "name": "PLKI", "stage": "PHASE_3"}]}
        drugs = {"PANHDACI": drug("Histone deacetylase inhibitor", targets=["HDAC1", "HDAC2", "HDAC3"], name="PANHDACI"),
                 "PLKI": drug("PLK1 inhibitor", targets=["PLK1"], name="PLKI"),
                 "GOODDRUG": drug("AMPK activator", action="ACTIVATOR", targets=["PRKAA1"], name="GOODDRUG")}
        labels = {"gooddrug": label(pediatric_use="Safety and effectiveness have been established in pediatric patients.")}
        panel = {"query_genes": {"GENEQ": {"disorder": "test", "modifier_targets": ["MOD"]}},
                 "named_drugs": ["gooddrug"], "negative_controls": ["PLK1"]}
        with mock.patch.object(S, "ot_target", lambda h, g: {"symbol": g, "ensembl": "ENSG0", "tractability": [], "candidates": cands.get(g, [])}), \
             mock.patch.object(S, "chembl_potency", lambda h, g: {"target": None, "records": 0, "potent": 0, "best_pchembl": None}), \
             mock.patch.object(S, "ot_drug", lambda h, name=None, chembl=None: drugs.get(name.upper())), \
             mock.patch.object(S, "fda_label", lambda h, n: labels.get(n.lower())):
            genes, rows = T.run(panel, http=None, log=lambda m: None)
        self.assertEqual([r["drug"] for r in rows], ["GOODDRUG", "PANHDACI", "PLKI"])
        self.assertEqual([r["tier"] for r in rows], ["screen-clear", "excluded", "excluded"])
        md = T.render_md(genes, rows, "2000-01-01")
        self.assertIn("| 1 | GOODDRUG |", md)
        with tempfile.TemporaryDirectory() as d:
            T.write_tsv(Path(d) / "x.tsv", rows)
            self.assertEqual(len((Path(d) / "x.tsv").read_text(encoding="utf-8").splitlines()), 4)

    def test_offline_cache_miss_raises(self):
        with tempfile.TemporaryDirectory() as d:
            with self.assertRaises(LookupError):
                S.Http(d, offline=True).request("https://example.invalid/x")


if __name__ == "__main__":
    unittest.main()
