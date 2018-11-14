import unittest
import os
from tiptoft.Gene import Gene


class TestGene(unittest.TestCase):

    def test_full_coverage(self):
        g = Gene('rep5.1_rep(pMW2)_NC005011', 10, 0)
        self.assertTrue(g.is_full_coverage())
        self.assertEqual(str(
            g),
            'rep5.1	Full	100	NC005011	plasmidfinder'
            '	rep5.1_rep(pMW2)_NC005011')	

    def test_no_coverage(self):
        g = Gene('rep5.1_rep(pMW2)_NC005011', 0, 10)
        self.assertFalse(g.is_full_coverage())
        self.assertEqual(str(
            g), 'rep5.1	Partial	0	NC005011	plasmidfinder'
            '	rep5.1_rep(pMW2)_NC005011')

    def test_medium_coverage(self):
        g = Gene('rep5.1_rep(pMW2)_NC005011', 5, 5)
        self.assertFalse(g.is_full_coverage())
        self.assertEqual(str(
            g), 'rep5.1	Partial	50	NC005011	plasmidfinder'
            '	rep5.1_rep(pMW2)_NC005011')

    def test_nc_name(self):
        g = Gene('Col(BS512).1__NC_010656.2', 5, 5)
        self.assertEqual(g.accession(), "NC_010656.2")
        self.assertEqual(g.short_name(), "Col(BS512).1")

    def test_rep_name(self):
        g = Gene('rep5.1_rep(pMW2)_NC005011', 5, 5)
        self.assertEqual(g.short_name(), "rep5.1")
        self.assertEqual(g.accession(), "NC005011")

    def test_rep_with_slash_name(self):
        g = Gene('rep6.1_repA(p703/5)_AF109375', 5, 5)
        self.assertEqual(g.short_name(), "rep6.1")
        self.assertEqual(g.accession(), "AF109375")

    def test_inc_double_dash_name(self):
        g = Gene('IncFII(S).1__CP000851', 5, 5)
        self.assertEqual(g.short_name(), "IncFII(S).1")
        self.assertEqual(g.accession(), "CP000851")
