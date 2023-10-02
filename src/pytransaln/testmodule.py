#!/usr/bin/env python3

import unittest

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from pytransaln.translate import translate_striptrailing


class TestTranslate(unittest.TestCase):
    s = SeqRecord(Seq('ATGTTGATAATATTTTGAT'), id='seq1', name='seq1')
    
    def test_translate_striptrailing(self):
        self.assertEqual(
            str(translate_striptrailing(TestTranslate.s, frame=0, table=5, id=None, name=None).seq),
            'MLMMFW'
        )
        self.assertEqual(
            str(translate_striptrailing(TestTranslate.s, frame=1, table=5, id=None, name=None).seq),
            'CW*YFD'
        )
        self.assertEqual(
            str(translate_striptrailing(TestTranslate.s, frame=-1, table=5, id=None, name=None).seq),
            'IKMLST'
        )
        self.assertEqual(
            str(translate_striptrailing(TestTranslate.s, frame=-2, table=5, id=None, name=None).seq),
            'SKYYQH'
        )
