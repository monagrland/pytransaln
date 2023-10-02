#!/usr/bin/env python3

import unittest

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from pytransaln.translate import translate_striptrailing
from pytransaln.frameshifts import contiguous_runs, report_deletions


class TestTranslate(unittest.TestCase):
    s = SeqRecord(Seq("ATGTTGATAATATTTTGAT"), id="seq1", name="seq1")

    def test_translate_striptrailing(self):
        # Forward frames
        self.assertEqual(
            str(
                translate_striptrailing(
                    TestTranslate.s, frame=0, table=5, id=None, name=None
                ).seq
            ),
            "MLMMFW",
        )
        self.assertEqual(
            str(
                translate_striptrailing(
                    TestTranslate.s, frame=1, table=5, id=None, name=None
                ).seq
            ),
            "CW*YFD",
        )
        # Reverse frames
        self.assertEqual(
            str(
                translate_striptrailing(
                    TestTranslate.s, frame=-1, table=5, id=None, name=None
                ).seq
            ),
            "IKMLST",
        )
        self.assertEqual(
            str(
                translate_striptrailing(
                    TestTranslate.s, frame=-2, table=5, id=None, name=None
                ).seq
            ),
            "SKYYQH",
        )


class TestFrameshifts(unittest.TestCase):
    op1 = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    rp1 = [-1, 1, 3, -1, 6, 7, 8, -1, 9, 10]

    def test_contiguous_runs(self):
        self.assertEqual(
            contiguous_runs([1, 2, 3, 6, 7, 8, 10, 11, 12], report="values"),
            [(1, 3), (6, 8), (10, 12)],
        )

    def test_report_deletions(self):
        self.assertEqual(
            report_deletions(refpos=TestFrameshifts.rp1, origpos=TestFrameshifts.op1),
            [
                {
                    "refstart": -1,
                    "refend": -1,
                    "origstart": 1,
                    "origend": 1,
                    "gaplen": 1,
                },
                {
                    "refstart": -1,
                    "refend": -1,
                    "origstart": 4,
                    "origend": 4,
                    "gaplen": 1,
                },
                {
                    "refstart": -1,
                    "refend": -1,
                    "origstart": 8,
                    "origend": 8,
                    "gaplen": 1,
                },
            ],
        )
